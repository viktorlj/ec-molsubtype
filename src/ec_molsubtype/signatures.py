"""Mutational signature concordance assessment."""

from __future__ import annotations

import json
from pathlib import Path
from typing import NamedTuple

from .models import MolecularSubtype


class SignatureResult(NamedTuple):
    """Signature concordance assessment result."""

    concordant: bool | None
    expected_signatures: list[str]
    present_signatures: list[str]
    signature_weights: dict[str, float]
    details: str


_SIG_DATA: dict | None = None


def _get_signature_data() -> dict:
    global _SIG_DATA
    if _SIG_DATA is None:
        data_path = Path(__file__).parent / "data" / "signature_profiles.json"
        with open(data_path) as f:
            _SIG_DATA = json.load(f)
    return _SIG_DATA


def assess_signatures(
    signature_weights: dict[str, float] | None,
    subtype: MolecularSubtype,
) -> SignatureResult:
    """Assess mutational signature concordance with expected profile."""
    data = _get_signature_data()
    profile = data["profiles"][subtype.value]
    expected = profile["expected_signatures"]
    threshold = profile["signature_threshold"]

    if signature_weights is None:
        return SignatureResult(
            concordant=None,
            expected_signatures=expected,
            present_signatures=[],
            signature_weights={},
            details="Signature weights not available.",
        )

    # Find expected signatures that are present above threshold
    present = [
        sig for sig in expected
        if signature_weights.get(sig, 0.0) >= threshold
    ]

    # Check concordance
    if subtype == MolecularSubtype.POLEmut:
        # Need at least one of SBS10a/SBS10b
        concordant = len(present) > 0
        sub_profile = profile.get("substitution_profile", {})
        # Additional check for substitution profile if available
        details_parts = [f"Expected signatures: {', '.join(expected)}."]
        if present:
            details_parts.append(f"Present: {', '.join(present)}.")
        else:
            details_parts.append("None of the expected POLE signatures detected.")
        details = " ".join(details_parts)
    elif subtype == MolecularSubtype.MMRd:
        concordant = len(present) > 0
        details = (
            f"Expected signatures: {', '.join(expected)}. "
            f"Present: {', '.join(present) if present else 'none'}."
        )
    else:
        # p53abn and NSMP: clock signatures expected, less strict
        concordant = len(present) > 0
        details = (
            f"Expected signatures: {', '.join(expected)}. "
            f"Present: {', '.join(present) if present else 'none'}."
        )

    return SignatureResult(
        concordant=concordant,
        expected_signatures=expected,
        present_signatures=present,
        signature_weights=signature_weights or {},
        details=details,
    )


def assess_substitution_profile(
    signature_weights: dict[str, float] | None,
    subtype: MolecularSubtype,
) -> dict[str, float | bool | None] | None:
    """Assess substitution class fractions for POLE concordance.

    This is specifically for POLEmut where C>A >20%, T>G >4%, C>G <0.6%.
    Returns None if not applicable or data not available.
    """
    if subtype != MolecularSubtype.POLEmut or signature_weights is None:
        return None

    # These would typically come from the actual mutation spectrum, not signatures.
    # For now, return None — this should be computed from the variant list directly.
    return None
