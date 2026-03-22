"""Substitution spectrum analysis from variant lists.

Computes the 6-class SBS spectrum (C>A, C>G, C>T, T>A, T>C, T>G) and indel
fraction from panel sequencing variants. This is a robust proxy for mutational
signatures that works even with limited variant counts from targeted panels.

POLE ultramutated tumors have a distinctive spectrum:
  - C>A > 20% (hallmark of SBS10a/SBS10b)
  - T>G > 4%
  - C>G < 0.6%
  - Indel fraction < 5%

MMRd tumors have:
  - High indel fraction (often > 20%)
  - Enrichment of T>C and C>T

Reference: León-Castillo et al. 2020, J Pathol 250:323-335.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import NamedTuple

from .models import MolecularSubtype, Variant

# Pyrimidine complement for SBS convention (mutations reported on the pyrimidine strand)
_COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# The 6 substitution classes (pyrimidine reference convention)
SBS_CLASSES = ["C>A", "C>G", "C>T", "T>A", "T>C", "T>G"]

# Variant classifications that represent indels
INDEL_CLASSIFICATIONS = {
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "In_Frame_Del",
    "In_Frame_Ins",
}

# SNV variant classifications
SNV_CLASSIFICATIONS = {
    "Missense_Mutation",
    "Nonsense_Mutation",
    "Silent",
    "Splice_Site",
    "Nonstop_Mutation",
    "Translation_Start_Site",
}


class SpectrumResult(NamedTuple):
    """Substitution spectrum analysis result."""

    fractions: dict[str, float]  # 6-class SBS fractions (sum to 1.0)
    counts: dict[str, int]       # Raw counts per substitution class
    n_snvs: int                  # Total SNVs used
    n_indels: int                # Indel count
    indel_fraction: float        # Indels / (SNVs + indels)
    n_total_variants: int        # All variants considered


class SpectrumConcordance(NamedTuple):
    """Concordance assessment of spectrum with expected subtype profile."""

    concordant: bool | None
    spectrum: SpectrumResult | None
    checks: dict[str, dict]   # Per-criterion check results
    details: str


def _to_pyrimidine_context(ref: str, alt: str) -> tuple[str, str] | None:
    """Convert a substitution to pyrimidine reference convention.

    Returns (ref, alt) in pyrimidine context, or None if not a valid SNV.
    """
    if len(ref) != 1 or len(alt) != 1:
        return None
    if ref not in "ACGT" or alt not in "ACGT" or ref == alt:
        return None

    if ref in ("C", "T"):
        return ref, alt
    # Complement both
    return _COMPLEMENT[ref], _COMPLEMENT[alt]


def compute_spectrum(variants: list[Variant]) -> SpectrumResult:
    """Compute 6-class substitution spectrum and indel fraction from variants.

    Only considers coding variants (SNVs and indels with recognized
    Variant_Classification). Ignores non-coding (UTR, intron, IGR, flanking).
    """
    counts = {cls: 0 for cls in SBS_CLASSES}
    n_snvs = 0
    n_indels = 0
    n_total = 0

    for v in variants:
        vc = v.variant_classification

        if vc in INDEL_CLASSIFICATIONS:
            n_indels += 1
            n_total += 1
            continue

        if vc not in SNV_CLASSIFICATIONS:
            continue

        n_total += 1

        ref = v.reference_allele.upper().strip()
        alt = v.tumor_seq_allele2.upper().strip()

        pyrim = _to_pyrimidine_context(ref, alt)
        if pyrim is None:
            continue

        key = f"{pyrim[0]}>{pyrim[1]}"
        if key in counts:
            counts[key] += 1
            n_snvs += 1

    # Compute fractions
    fractions = {}
    for cls in SBS_CLASSES:
        fractions[cls] = counts[cls] / n_snvs if n_snvs > 0 else 0.0

    total_snv_indel = n_snvs + n_indels
    indel_fraction = n_indels / total_snv_indel if total_snv_indel > 0 else 0.0

    return SpectrumResult(
        fractions=fractions,
        counts=counts,
        n_snvs=n_snvs,
        n_indels=n_indels,
        indel_fraction=indel_fraction,
        n_total_variants=n_total,
    )


def _load_spectrum_thresholds() -> dict:
    """Load expected substitution profile thresholds from signature_profiles.json."""
    data_path = Path(__file__).parent / "data" / "signature_profiles.json"
    with open(data_path) as f:
        data = json.load(f)
    return data


def assess_spectrum(
    variants: list[Variant],
    subtype: MolecularSubtype,
    min_snvs: int = 10,
) -> SpectrumConcordance:
    """Assess substitution spectrum concordance with expected subtype profile.

    Args:
        variants: List of variants to analyze.
        subtype: The classified subtype to check concordance against.
        min_snvs: Minimum SNVs required for reliable spectrum analysis.
            Below this threshold, result is indeterminate.

    Returns:
        SpectrumConcordance with concordance assessment and details.
    """
    spectrum = compute_spectrum(variants)

    if spectrum.n_snvs < min_snvs:
        return SpectrumConcordance(
            concordant=None,
            spectrum=spectrum,
            checks={},
            details=f"Insufficient SNVs for spectrum analysis ({spectrum.n_snvs} < {min_snvs}).",
        )

    data = _load_spectrum_thresholds()
    profile = data["profiles"].get(subtype.value, {})
    sub_profile = profile.get("substitution_profile", {})

    if not sub_profile:
        return SpectrumConcordance(
            concordant=None,
            spectrum=spectrum,
            checks={},
            details=f"No substitution profile thresholds defined for {subtype.value}.",
        )

    checks: dict[str, dict] = {}
    all_pass = True

    # Evaluate each criterion in the substitution profile
    for criterion, spec in sub_profile.items():
        if criterion == "C>A_min":
            observed = spectrum.fractions.get("C>A", 0.0)
            passed = observed >= spec
            checks["C>A"] = {
                "observed": round(observed, 4),
                "threshold": spec,
                "direction": ">=",
                "passed": passed,
            }
            if not passed:
                all_pass = False

        elif criterion == "T>G_min":
            observed = spectrum.fractions.get("T>G", 0.0)
            passed = observed >= spec
            checks["T>G"] = {
                "observed": round(observed, 4),
                "threshold": spec,
                "direction": ">=",
                "passed": passed,
            }
            if not passed:
                all_pass = False

        elif criterion == "C>G_max":
            observed = spectrum.fractions.get("C>G", 0.0)
            passed = observed <= spec
            checks["C>G"] = {
                "observed": round(observed, 4),
                "threshold": spec,
                "direction": "<=",
                "passed": passed,
            }
            if not passed:
                all_pass = False

        elif criterion == "indel_fraction_max":
            observed = spectrum.indel_fraction
            passed = observed <= spec
            checks["indel_fraction"] = {
                "observed": round(observed, 4),
                "threshold": spec,
                "direction": "<=",
                "passed": passed,
            }
            if not passed:
                all_pass = False

        elif criterion == "indel_fraction_min":
            observed = spectrum.indel_fraction
            passed = observed >= spec
            checks["indel_fraction"] = {
                "observed": round(observed, 4),
                "threshold": spec,
                "direction": ">=",
                "passed": passed,
            }
            if not passed:
                all_pass = False

    # Build details string
    parts = []
    for name, check in checks.items():
        status = "PASS" if check["passed"] else "FAIL"
        parts.append(
            f"{name}: {check['observed']:.1%} {check['direction']} "
            f"{check['threshold']:.1%} [{status}]"
        )
    details = "; ".join(parts) if parts else "No spectrum criteria to evaluate."

    return SpectrumConcordance(
        concordant=all_pass,
        spectrum=spectrum,
        checks=checks,
        details=details,
    )
