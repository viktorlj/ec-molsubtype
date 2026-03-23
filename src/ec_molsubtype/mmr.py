"""MMR deficiency assessment for molecular subtyping."""

from __future__ import annotations

import json
from pathlib import Path
from typing import NamedTuple

from .models import ConfidenceLevel, Variant

MMR_GENES = {"MLH1", "MSH2", "MSH6", "PMS2", "EPCAM"}

DEFAULT_MSI_THRESHOLD = 20.0


class MmrResult(NamedTuple):
    """Result of MMR deficiency assessment."""

    is_mmrd: bool
    confidence: ConfidenceLevel | None
    msi_pct: float | None
    msi_status: str | None  # MSI-H, MSI-L, MSS
    mmr_mutations: list[str]
    is_biallelic: bool
    flags: list[str]
    clinical_notes: list[str]


def _load_mmr_data() -> dict:
    data_path = Path(__file__).parent / "data" / "mmr_genes.json"
    with open(data_path) as f:
        return json.load(f)


def classify_msi_status(msi_pct: float, threshold: float = DEFAULT_MSI_THRESHOLD) -> str:
    """Classify MSI status from percentage of unstable markers."""
    if msi_pct >= threshold:
        return "MSI-H"
    elif msi_pct > 0:
        return "MSI-L"
    return "MSS"


def find_mmr_variants(variants: list[Variant]) -> list[Variant]:
    """Filter variants to those in MMR genes."""
    return [v for v in variants if v.hugo_symbol in MMR_GENES]


def is_pathogenic_mmr_variant(variant: Variant) -> bool:
    """Check if an MMR variant is likely pathogenic.

    Accepts:
    - Truncating variants (nonsense, frameshift, splice-site)
    - Variants with ClinVar pathogenic/likely_pathogenic classification
    """
    if variant.is_truncating:
        return True
    clinvar = variant.clinvar_classification.lower()
    if "pathogenic" in clinvar and "benign" not in clinvar:
        return True
    return False


def assess_biallelic(mmr_variants: list[Variant]) -> dict[str, list[Variant]]:
    """Group pathogenic MMR variants by gene to assess biallelic inactivation.

    Returns dict of gene -> list of pathogenic variants.
    A gene with >=2 pathogenic variants suggests biallelic inactivation.
    """
    gene_hits: dict[str, list[Variant]] = {}
    for v in mmr_variants:
        if is_pathogenic_mmr_variant(v):
            gene_hits.setdefault(v.hugo_symbol, []).append(v)
    return gene_hits


def assess_mmr(
    variants: list[Variant],
    msi_pct: float | None = None,
    msi_status_override: str | None = None,
    msi_threshold: float = DEFAULT_MSI_THRESHOLD,
    mmr_ihc: dict[str, str | None] | None = None,
) -> MmrResult:
    """Assess MMR deficiency from IHC, MSI data, and/or MMR gene mutations.

    Priority order:
    0. MMR IHC protein loss → MMRd (high confidence, clinical gold standard)
    1. msi_pct >= threshold → MMRd (high confidence)
    2. msi_pct < threshold BUT biallelic MMR → MMRd (moderate, flag discordance)
    3. msi_pct not provided AND biallelic MMR → MMRd (moderate)
    4. msi_pct not provided AND single het MMR → NOT MMRd, flag
    5. No MSI data and no MMR mutations → MMR proficient
    """
    flags: list[str] = []
    clinical_notes: list[str] = []

    # --- Step 0: MMR IHC (highest priority when available) ---
    if mmr_ihc:
        lost_proteins = [
            gene.upper() for gene, status in mmr_ihc.items()
            if status and status.lower() == "lost"
        ]
        if lost_proteins:
            lost_str = ", ".join(sorted(lost_proteins))
            clinical_notes.append(
                f"MMR IHC: loss of {lost_str} protein expression. "
                "Direct evidence of MMR deficiency."
            )
            # Determine Lynch screening guidance based on which protein is lost
            if "MLH1" in lost_proteins:
                clinical_notes.append(
                    "MLH1 loss detected — check MLH1 promoter hypermethylation. "
                    "If negative, germline MLH1 testing recommended (Lynch syndrome)."
                )
            if "MSH2" in lost_proteins:
                clinical_notes.append(
                    "MSH2 loss detected — germline MSH2/EPCAM testing recommended (Lynch syndrome)."
                )
            if "MSH6" in lost_proteins and "MSH2" not in lost_proteins:
                clinical_notes.append(
                    "Isolated MSH6 loss — germline MSH6 testing recommended."
                )
            if "PMS2" in lost_proteins and "MLH1" not in lost_proteins:
                clinical_notes.append(
                    "Isolated PMS2 loss — germline PMS2 testing recommended."
                )

            # Check for supporting sequencing evidence
            mmr_variants = find_mmr_variants(variants)
            pathogenic_by_gene = assess_biallelic(mmr_variants)
            all_pathogenic = [v for vs in pathogenic_by_gene.values() for v in vs]
            mmr_mutation_strs = [f"{v.hugo_symbol} {v.hgvsp_short}" for v in all_pathogenic]
            if mmr_mutation_strs:
                clinical_notes.append(
                    f"Supporting MMR gene mutations: {', '.join(mmr_mutation_strs)}."
                )

            return MmrResult(
                is_mmrd=True,
                confidence=ConfidenceLevel.HIGH,
                msi_pct=msi_pct,
                msi_status=msi_status_override or (classify_msi_status(msi_pct, msi_threshold) if msi_pct is not None else None),
                mmr_mutations=mmr_mutation_strs,
                is_biallelic=any(len(vs) >= 2 for vs in pathogenic_by_gene.values()),
                flags=flags,
                clinical_notes=clinical_notes,
            )

        # All tested proteins intact
        intact_proteins = [
            gene.upper() for gene, status in mmr_ihc.items()
            if status and status.lower() == "intact"
        ]
        if intact_proteins:
            clinical_notes.append(
                f"MMR IHC: {', '.join(sorted(intact_proteins))} intact."
            )

    # --- MSI / mutation-based assessment (original logic) ---
    computed_msi_status: str | None = None
    if msi_pct is not None:
        computed_msi_status = classify_msi_status(msi_pct, msi_threshold)

        # Validate against override if provided
        if msi_status_override and computed_msi_status != msi_status_override:
            flags.append(
                f"MSI status discordance: computed {computed_msi_status} from "
                f"msi_pct={msi_pct}%, but override says {msi_status_override}"
            )
    elif msi_status_override:
        computed_msi_status = msi_status_override

    # MMR gene mutation assessment
    mmr_variants = find_mmr_variants(variants)
    pathogenic_by_gene = assess_biallelic(mmr_variants)
    all_pathogenic = [v for vs in pathogenic_by_gene.values() for v in vs]
    mmr_mutation_strs = [
        f"{v.hugo_symbol} {v.hgvsp_short}" for v in all_pathogenic
    ]

    has_biallelic = any(len(vs) >= 2 for vs in pathogenic_by_gene.values())
    has_single_het = len(all_pathogenic) > 0 and not has_biallelic

    # Combined logic
    if msi_pct is not None and msi_pct >= msi_threshold:
        # MSI-H from marker panel → MMRd high confidence
        clinical_notes.append(
            f"MSI-H detected: {msi_pct:.1f}% unstable markers (threshold: {msi_threshold}%)."
        )
        if mmr_mutation_strs:
            clinical_notes.append(
                f"Supporting MMR gene mutations: {', '.join(mmr_mutation_strs)}."
            )
        clinical_notes.append(
            "IHC for MLH1/MSH2/MSH6/PMS2 recommended for Lynch syndrome screening."
        )
        return MmrResult(
            is_mmrd=True,
            confidence=ConfidenceLevel.HIGH,
            msi_pct=msi_pct,
            msi_status=computed_msi_status,
            mmr_mutations=mmr_mutation_strs,
            is_biallelic=has_biallelic,
            flags=flags,
            clinical_notes=clinical_notes,
        )

    if msi_pct is not None and msi_pct < msi_threshold and has_biallelic:
        # MSI not high but biallelic MMR → moderate with discordance flag
        flags.append(
            f"Discordance: MSI-L/MSS ({msi_pct:.1f}%) but biallelic MMR gene inactivation detected."
        )
        clinical_notes.append(
            "Biallelic MMR gene inactivation detected despite low MSI marker percentage."
        )
        clinical_notes.append("IHC confirmation strongly recommended.")
        return MmrResult(
            is_mmrd=True,
            confidence=ConfidenceLevel.MODERATE,
            msi_pct=msi_pct,
            msi_status=computed_msi_status,
            mmr_mutations=mmr_mutation_strs,
            is_biallelic=True,
            flags=flags,
            clinical_notes=clinical_notes,
        )

    if msi_pct is None and has_biallelic:
        # No MSI data but biallelic MMR → moderate
        clinical_notes.append(
            "No MSI marker data available. Biallelic MMR gene inactivation supports MMRd."
        )
        clinical_notes.append("IHC confirmation recommended.")
        return MmrResult(
            is_mmrd=True,
            confidence=ConfidenceLevel.MODERATE,
            msi_pct=None,
            msi_status=computed_msi_status,
            mmr_mutations=mmr_mutation_strs,
            is_biallelic=True,
            flags=flags,
            clinical_notes=clinical_notes,
        )

    if msi_pct is None and has_single_het:
        # No MSI data and only single het → NOT MMRd, flag
        flags.append("Possible MMRd — IHC recommended.")
        flags.append("Single heterozygous MMR variant — possible Lynch carrier.")
        clinical_notes.append(
            f"Single heterozygous MMR variant(s): {', '.join(mmr_mutation_strs)}. "
            "Does not confirm somatic MMRd. IHC and germline testing recommended."
        )
        return MmrResult(
            is_mmrd=False,
            confidence=None,
            msi_pct=None,
            msi_status=computed_msi_status,
            mmr_mutations=mmr_mutation_strs,
            is_biallelic=False,
            flags=flags,
            clinical_notes=clinical_notes,
        )

    # MSI-L or MSS without supporting MMR mutations
    if msi_pct is not None and msi_pct < msi_threshold:
        if has_single_het:
            flags.append(
                "Single heterozygous MMR variant with MSS/MSI-L. "
                "Possible germline carrier — IHC and germline testing recommended."
            )
        return MmrResult(
            is_mmrd=False,
            confidence=None,
            msi_pct=msi_pct,
            msi_status=computed_msi_status,
            mmr_mutations=mmr_mutation_strs,
            is_biallelic=False,
            flags=flags,
            clinical_notes=clinical_notes,
        )

    # No MSI data and no MMR mutations → MMR proficient
    return MmrResult(
        is_mmrd=False,
        confidence=None,
        msi_pct=None,
        msi_status=computed_msi_status,
        mmr_mutations=[],
        is_biallelic=False,
        flags=flags,
        clinical_notes=clinical_notes,
    )
