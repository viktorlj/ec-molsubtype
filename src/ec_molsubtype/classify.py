"""Main classification orchestrator implementing the WHO 2020 sequential algorithm."""

from __future__ import annotations

from .evidence import update_result_with_evidence
from .mmr import DEFAULT_MSI_THRESHOLD, assess_mmr
from .models import (
    ClassificationResult,
    ClassificationStep,
    ConfidenceLevel,
    MolecularSubtype,
    MultipleClassifier,
    PoleTier,
    SampleInput,
)
from .pole import assess_pole, get_best_pole_result
from .tp53 import assess_tp53, get_pathogenic_tp53


def classify_sample(
    sample: SampleInput,
    msi_threshold: float = DEFAULT_MSI_THRESHOLD,
    min_tp53_vaf: float = 0.05,
) -> ClassificationResult:
    """Classify a sample into one of four molecular subtypes.

    Implements the WHO 2020 hierarchical algorithm:
    1. POLE pathogenic EDM → POLEmut
    2. MMR deficiency → MMRd
    3. TP53 pathogenic mutation → p53abn
    4. None of the above → NSMP

    Args:
        sample: Complete sample input with variants and metadata.
        msi_threshold: MSI percentage threshold for MSI-H (default 20%).
        min_tp53_vaf: Minimum VAF for TP53 variant consideration.

    Returns:
        ClassificationResult with subtype, confidence, and evidence.
    """
    meta = sample.metadata
    variants = sample.variants
    path: list[ClassificationStep] = []
    clinical_notes: list[str] = []
    flags: list[str] = []

    # Assess all three classifiers upfront (needed for multiple-classifier detection)
    pole_results = assess_pole(variants)
    best_pole = get_best_pole_result(pole_results)

    mmr_result = assess_mmr(
        variants,
        msi_pct=meta.msi_pct,
        msi_status_override=meta.msi_status,
        msi_threshold=msi_threshold,
    )

    tp53_results = assess_tp53(variants, min_vaf=min_tp53_vaf)
    best_tp53 = get_pathogenic_tp53(tp53_results)

    # --- Step 1: POLE EDM ---
    if best_pole and best_pole.tier in (PoleTier.TIER1, PoleTier.TIER2):
        path.append(ClassificationStep(
            step=1,
            test="POLE_EDM",
            result="positive",
            variant=best_pole.variant_str,
            tier=best_pole.tier.value,
        ))

        clinical_notes.append(
            f"{best_pole.variant_str.split(' ')[-1]} is a {best_pole.tier.value} "
            f"{'established pathogenic hotspot' if best_pole.tier == PoleTier.TIER1 else 'confirmed pathogenic variant'}."
        )

        # Check multiple classifier
        mc = _detect_multiple_classifier(
            MolecularSubtype.POLEmut, best_pole, mmr_result, best_tp53
        )
        if mc.is_multiple:
            clinical_notes.extend(_multiple_classifier_notes(MolecularSubtype.POLEmut, mc))

        clinical_notes.append(
            "Swedish guidelines: no adjuvant treatment recommended for stage I-II POLEmut."
        )

        result = ClassificationResult(
            sample_id=meta.sample_id,
            primary_subtype=MolecularSubtype.POLEmut,
            confidence=ConfidenceLevel.HIGH,
            classification_path=path,
            multiple_classifier=mc,
            clinical_notes=clinical_notes,
            flags=flags + mmr_result.flags,
        )

        return _finalize(result, meta, variants)

    elif best_pole and best_pole.tier == PoleTier.VUS:
        # VUS — do NOT auto-classify, but flag
        path.append(ClassificationStep(
            step=1,
            test="POLE_EDM",
            result="negative",
            variant=best_pole.variant_str,
            details="POLE VUS in exonuclease domain — requires secondary evidence for POLEmut classification.",
        ))
        flags.append(
            f"POLE VUS detected: {best_pole.variant_str}. "
            "Secondary evidence (TMB >100, SBS10a/b) needed for POLEmut classification."
        )
    else:
        path.append(ClassificationStep(step=1, test="POLE_EDM", result="negative"))

    # --- Step 2: MMR deficiency ---
    if mmr_result.is_mmrd:
        mmr_step = ClassificationStep(
            step=2,
            test="MMRd",
            result="positive",
            msi_pct=mmr_result.msi_pct,
            msi_threshold=msi_threshold if mmr_result.msi_pct is not None else None,
            mmr_mutations=mmr_result.mmr_mutations if mmr_result.mmr_mutations else None,
            confidence=mmr_result.confidence.value if mmr_result.confidence else None,
        )
        path.append(mmr_step)
        clinical_notes.extend(mmr_result.clinical_notes)
        flags.extend(mmr_result.flags)

        mc = _detect_multiple_classifier(
            MolecularSubtype.MMRd, best_pole, mmr_result, best_tp53
        )
        if mc.is_multiple:
            clinical_notes.extend(_multiple_classifier_notes(MolecularSubtype.MMRd, mc))

        result = ClassificationResult(
            sample_id=meta.sample_id,
            primary_subtype=MolecularSubtype.MMRd,
            confidence=mmr_result.confidence or ConfidenceLevel.MODERATE,
            classification_path=path,
            multiple_classifier=mc,
            clinical_notes=clinical_notes,
            flags=flags,
        )
        return _finalize(result, meta, variants)
    else:
        path.append(ClassificationStep(step=2, test="MMRd", result="negative"))
        flags.extend(mmr_result.flags)
        clinical_notes.extend(mmr_result.clinical_notes)

    # --- Step 3: TP53 pathogenic mutation ---
    if best_tp53 and best_tp53.is_pathogenic:
        path.append(ClassificationStep(
            step=3,
            test="TP53_mutation",
            result="positive",
            variant=best_tp53.variant_str,
        ))
        clinical_notes.extend(best_tp53.clinical_notes)
        clinical_notes.append(
            "Swedish guidelines: postoperative chemotherapy for stage I-II p53abn."
        )
        clinical_notes.append(
            "Note: TP53 sequencing used here — p53 IHC confirmation recommended "
            "(~92% concordance with sequencing)."
        )

        result = ClassificationResult(
            sample_id=meta.sample_id,
            primary_subtype=MolecularSubtype.p53abn,
            confidence=ConfidenceLevel.HIGH,
            classification_path=path,
            clinical_notes=clinical_notes,
            flags=flags,
        )
        return _finalize(result, meta, variants)
    else:
        details = None
        if tp53_results:
            non_pathogenic = [r for r in tp53_results if not r.is_pathogenic and not r.is_benign_polymorphism]
            if non_pathogenic:
                details = f"TP53 variant(s) found but not classified as pathogenic."
        path.append(ClassificationStep(
            step=3, test="TP53_mutation", result="negative", details=details
        ))
        # Add any TP53 clinical notes (e.g., VUS in DBD)
        for r in tp53_results:
            clinical_notes.extend(r.clinical_notes)

    # --- Step 4: NSMP ---
    path.append(ClassificationStep(step=4, test="NSMP", result="positive"))
    clinical_notes.append(
        "No specific molecular profile identified. Classified as NSMP by exclusion."
    )

    result = ClassificationResult(
        sample_id=meta.sample_id,
        primary_subtype=MolecularSubtype.NSMP,
        confidence=ConfidenceLevel.HIGH,
        classification_path=path,
        clinical_notes=clinical_notes,
        flags=flags,
    )
    return _finalize(result, meta)


def _finalize(result: ClassificationResult, meta, variants=None) -> ClassificationResult:
    """Apply secondary evidence and finalize the result."""
    return update_result_with_evidence(
        result,
        tmb=meta.tmb,
        fraction_genome_altered=meta.fraction_genome_altered,
        signature_weights=meta.signature_weights,
        variants=variants,
    )


def _detect_multiple_classifier(
    primary: MolecularSubtype,
    best_pole,
    mmr_result,
    best_tp53,
) -> MultipleClassifier:
    """Detect if sample has multiple classifier features."""
    secondary_features = []
    mc = MultipleClassifier()

    if primary != MolecularSubtype.POLEmut and best_pole and best_pole.is_pathogenic:
        secondary_features.append("POLEmut")
        mc.pole_variant = best_pole.variant_str

    if primary != MolecularSubtype.MMRd and mmr_result.is_mmrd:
        secondary_features.append("MMRd")
        mc.mmr_evidence = (
            f"MSI-H ({mmr_result.msi_pct:.1f}%)" if mmr_result.msi_pct
            else "biallelic MMR mutations"
        )

    if primary != MolecularSubtype.p53abn and best_tp53 and best_tp53.is_pathogenic:
        secondary_features.append("p53abn")
        mc.tp53_variant = best_tp53.variant_str

    if secondary_features:
        mc.is_multiple = True
        mc.secondary_features = secondary_features

    return mc


def _multiple_classifier_notes(
    primary: MolecularSubtype,
    mc: MultipleClassifier,
) -> list[str]:
    """Generate clinical notes for multiple-classifier cases."""
    notes = []
    features_str = ", ".join(mc.secondary_features)

    notes.append(
        f"Multiple classifier: concurrent {features_str} detected. "
        f"Per WHO 2020 hierarchy, classified as {primary.value}."
    )

    # Detail each secondary feature
    if "MMRd" in mc.secondary_features:
        mmr_detail = mc.mmr_evidence or "MMR deficiency"
        notes.append(f"Concurrent MMRd evidence: {mmr_detail}.")
    if "p53abn" in mc.secondary_features:
        tp53_detail = mc.tp53_variant or "pathogenic TP53 variant"
        notes.append(f"Concurrent p53abn evidence: {tp53_detail}.")
    if "POLEmut" in mc.secondary_features:
        pole_detail = mc.pole_variant or "pathogenic POLE EDM"
        notes.append(f"Concurrent POLEmut evidence: {pole_detail}.")

    # Prognostic interpretation based on published evidence
    if primary == MolecularSubtype.POLEmut and "p53abn" in mc.secondary_features:
        notes.append(
            "POLEmut+p53abn: evidence shows prognosis follows POLEmut "
            "(León-Castillo et al. 2020). POLE status dominates."
        )
    if primary == MolecularSubtype.POLEmut and "MMRd" in mc.secondary_features:
        notes.append(
            "POLEmut+MMRd: classify as POLEmut per hierarchy. "
            "Note: patient may still benefit from immunotherapy evaluation."
        )
    if primary == MolecularSubtype.MMRd and "p53abn" in mc.secondary_features:
        notes.append(
            "MMRd+p53abn: evidence shows prognosis follows MMRd "
            "(León-Castillo et al. 2020). Consider immunotherapy."
        )

    return notes
