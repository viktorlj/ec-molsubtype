"""FIGO 2023 molecular staging annotation."""

from __future__ import annotations

from .models import MolecularSubtype


# FIGO 2023 introduced molecular classification as part of staging
# These annotations are appended to the traditional stage
FIGO_MOLECULAR_ANNOTATIONS = {
    MolecularSubtype.POLEmut: "IAmPOLEmut",
    MolecularSubtype.MMRd: "IAmMMRd",
    MolecularSubtype.p53abn: "IIICp53abn",
    MolecularSubtype.NSMP: "",  # No specific FIGO annotation for NSMP
}

FIGO_DESCRIPTIONS = {
    MolecularSubtype.POLEmut: (
        "FIGO 2023: Stage IA with POLEmut — favorable prognosis regardless of other features. "
        "Molecularly classified low-risk."
    ),
    MolecularSubtype.MMRd: (
        "FIGO 2023: MMRd molecular classification. "
        "Predictive for checkpoint inhibitor benefit."
    ),
    MolecularSubtype.p53abn: (
        "FIGO 2023: p53abn upgrades staging to at least IIIC2 equivalent. "
        "Aggressive biology — consider intensified treatment."
    ),
    MolecularSubtype.NSMP: (
        "FIGO 2023: No specific molecular annotation. "
        "Standard staging and risk stratification apply."
    ),
}


def get_figo_annotation(subtype: MolecularSubtype) -> str:
    """Get FIGO 2023 molecular staging annotation for a subtype."""
    return FIGO_MOLECULAR_ANNOTATIONS.get(subtype, "")


def get_figo_description(subtype: MolecularSubtype) -> str:
    """Get clinical description of FIGO 2023 molecular staging impact."""
    return FIGO_DESCRIPTIONS.get(subtype, "")
