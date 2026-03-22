"""Secondary evidence aggregation and confidence scoring."""

from __future__ import annotations

from .cna import assess_cna
from .models import (
    ClassificationResult,
    ConfidenceLevel,
    EvidenceItem,
    MolecularSubtype,
    SecondaryEvidence,
    Variant,
)
from .signatures import assess_signatures
from .spectrum import assess_spectrum
from .tmb import assess_tmb


def compute_secondary_evidence(
    subtype: MolecularSubtype,
    tmb: float | None = None,
    fraction_genome_altered: float | None = None,
    signature_weights: dict[str, float] | None = None,
    variants: list[Variant] | None = None,
) -> SecondaryEvidence:
    """Compute all secondary evidence for a given subtype."""
    tmb_result = assess_tmb(tmb, subtype)
    cna_result = assess_cna(fraction_genome_altered, subtype)
    sig_result = assess_signatures(signature_weights, subtype)

    # Compute substitution spectrum from variants (always available if variants provided)
    spectrum_item = None
    if variants is not None:
        spec_result = assess_spectrum(variants, subtype)
        spectrum_item = EvidenceItem(
            concordant=spec_result.concordant,
            details={
                "description": spec_result.details,
                "checks": spec_result.checks,
                "n_snvs": spec_result.spectrum.n_snvs if spec_result.spectrum else 0,
                "n_indels": spec_result.spectrum.n_indels if spec_result.spectrum else 0,
                "indel_fraction": (
                    round(spec_result.spectrum.indel_fraction, 4)
                    if spec_result.spectrum else None
                ),
                "fractions": (
                    {k: round(v, 4) for k, v in spec_result.spectrum.fractions.items()}
                    if spec_result.spectrum else {}
                ),
            },
        )

    return SecondaryEvidence(
        tmb=EvidenceItem(
            value=tmb_result.value,
            expected_range=tmb_result.expected_range,
            concordant=tmb_result.concordant,
            details={"description": tmb_result.details},
        ),
        cna_burden=EvidenceItem(
            value=cna_result.value,
            expected=cna_result.expected,
            concordant=cna_result.concordant,
            details={"description": cna_result.details},
        ),
        signatures=EvidenceItem(
            concordant=sig_result.concordant,
            details={
                "expected": sig_result.expected_signatures,
                "present": sig_result.present_signatures,
                "weights": sig_result.signature_weights,
                "description": sig_result.details,
            },
        ),
        substitution_profile=spectrum_item,
    )


def compute_confidence(
    primary_confidence: ConfidenceLevel,
    secondary_evidence: SecondaryEvidence,
) -> ConfidenceLevel:
    """Compute overall confidence based on primary classification and secondary evidence.

    Rules:
    - High: primary classifier + all concordant secondary evidence
    - Moderate: primary present but some discordance
    - Low: borderline cases (VUS, single het MMR)
    - Discordant: primary contradicted by secondary evidence
    """
    if primary_confidence == ConfidenceLevel.LOW:
        return ConfidenceLevel.LOW

    if primary_confidence == ConfidenceLevel.DISCORDANT:
        return ConfidenceLevel.DISCORDANT

    # Count available and concordant evidence
    evidence_items = [
        secondary_evidence.tmb,
        secondary_evidence.cna_burden,
        secondary_evidence.signatures,
        secondary_evidence.substitution_profile,
    ]

    available = [e for e in evidence_items if e is not None and e.concordant is not None]

    if not available:
        # No secondary evidence to evaluate
        return primary_confidence

    concordant_count = sum(1 for e in available if e.concordant)
    discordant_count = sum(1 for e in available if not e.concordant)

    if discordant_count == 0:
        return ConfidenceLevel.HIGH
    elif concordant_count > discordant_count:
        return ConfidenceLevel.MODERATE
    elif discordant_count >= concordant_count and primary_confidence == ConfidenceLevel.HIGH:
        return ConfidenceLevel.MODERATE
    else:
        return ConfidenceLevel.DISCORDANT


def update_result_with_evidence(
    result: ClassificationResult,
    tmb: float | None = None,
    fraction_genome_altered: float | None = None,
    signature_weights: dict[str, float] | None = None,
    variants: list[Variant] | None = None,
) -> ClassificationResult:
    """Add secondary evidence and update confidence for a classification result."""
    secondary = compute_secondary_evidence(
        result.primary_subtype,
        tmb=tmb,
        fraction_genome_altered=fraction_genome_altered,
        signature_weights=signature_weights,
        variants=variants,
    )
    result.secondary_evidence = secondary

    new_confidence = compute_confidence(result.confidence, secondary)
    result.confidence = new_confidence

    return result
