"""Tumor mutational burden (TMB) calculation and assessment."""

from __future__ import annotations

from typing import NamedTuple

from .models import MolecularSubtype, Variant


class TmbResult(NamedTuple):
    """TMB assessment result."""

    value: float | None
    concordant: bool | None
    expected_range: str
    details: str


# Expected TMB ranges per subtype
TMB_EXPECTATIONS = {
    MolecularSubtype.POLEmut: {"min": 100.0, "typical": 200.0, "label": ">100 mut/Mb"},
    MolecularSubtype.MMRd: {"min": 10.0, "max": 100.0, "label": "10-100 mut/Mb"},
    MolecularSubtype.p53abn: {"max": 10.0, "label": "<10 mut/Mb"},
    MolecularSubtype.NSMP: {"max": 10.0, "label": "<10 mut/Mb"},
}


def compute_tmb(
    variants: list[Variant],
    panel_size_mb: float,
    count_silent: bool = False,
) -> float:
    """Compute TMB from variant count and panel size.

    Args:
        variants: List of variants.
        panel_size_mb: Panel size in megabases.
        count_silent: Whether to count silent mutations (default False).
    """
    if panel_size_mb <= 0:
        raise ValueError(f"Panel size must be positive, got {panel_size_mb}")

    if count_silent:
        count = len(variants)
    else:
        count = sum(
            1 for v in variants if v.variant_classification != "Silent"
        )

    return count / panel_size_mb


def assess_tmb(
    tmb: float | None,
    subtype: MolecularSubtype,
) -> TmbResult:
    """Assess TMB concordance with expected range for a given subtype."""
    expectations = TMB_EXPECTATIONS[subtype]

    if tmb is None:
        return TmbResult(
            value=None,
            concordant=None,
            expected_range=expectations["label"],
            details="TMB not available.",
        )

    concordant = True
    if subtype == MolecularSubtype.POLEmut:
        concordant = tmb >= expectations["min"]
        if not concordant:
            details = (
                f"TMB {tmb:.1f} mut/Mb is lower than expected for POLEmut "
                f"(expected {expectations['label']}). Possible passenger/artifact."
            )
        else:
            details = f"TMB {tmb:.1f} mut/Mb concordant with POLEmut phenotype."
    elif subtype == MolecularSubtype.MMRd:
        concordant = expectations["min"] <= tmb <= expectations["max"]
        details = (
            f"TMB {tmb:.1f} mut/Mb "
            f"{'concordant' if concordant else 'outside expected range'} "
            f"for MMRd (expected {expectations['label']})."
        )
    else:
        concordant = tmb <= expectations["max"]
        details = (
            f"TMB {tmb:.1f} mut/Mb "
            f"{'concordant' if concordant else 'higher than expected'} "
            f"for {subtype.value} (expected {expectations['label']})."
        )

    return TmbResult(
        value=tmb,
        concordant=concordant,
        expected_range=expectations["label"],
        details=details,
    )
