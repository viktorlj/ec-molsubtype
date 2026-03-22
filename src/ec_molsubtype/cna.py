"""Copy number alteration (CNA) burden assessment."""

from __future__ import annotations

from typing import NamedTuple

from .models import MolecularSubtype


class CnaResult(NamedTuple):
    """CNA burden assessment result."""

    value: float | None
    concordant: bool | None
    expected: str
    details: str


CNA_EXPECTATIONS = {
    MolecularSubtype.POLEmut: {"max": 0.2, "label": "low"},
    MolecularSubtype.MMRd: {"max": 0.3, "label": "low-intermediate"},
    MolecularSubtype.p53abn: {"min": 0.3, "label": "high"},
    MolecularSubtype.NSMP: {"max": 0.2, "label": "low"},
}


def assess_cna(
    fraction_genome_altered: float | None,
    subtype: MolecularSubtype,
) -> CnaResult:
    """Assess CNA burden concordance with expected pattern for a given subtype."""
    expectations = CNA_EXPECTATIONS[subtype]

    if fraction_genome_altered is None:
        return CnaResult(
            value=None,
            concordant=None,
            expected=expectations["label"],
            details="CNA burden (fraction genome altered) not available.",
        )

    fga = fraction_genome_altered

    if subtype == MolecularSubtype.p53abn:
        concordant = fga >= expectations["min"]
    elif "max" in expectations:
        concordant = fga <= expectations["max"]
    else:
        concordant = True

    label = "concordant" if concordant else "discordant"
    details = (
        f"Fraction genome altered: {fga:.2f} — {label} with "
        f"{subtype.value} (expected {expectations['label']})."
    )

    return CnaResult(
        value=fga,
        concordant=concordant,
        expected=expectations["label"],
        details=details,
    )
