"""POLE exonuclease domain mutation (EDM) assessment for molecular subtyping."""

from __future__ import annotations

import json
import re
from importlib import resources
from pathlib import Path
from typing import NamedTuple

from .models import PoleTier, Variant


class PoleResult(NamedTuple):
    """Result of POLE EDM assessment."""

    is_pathogenic: bool
    tier: PoleTier
    variant_str: str
    needs_secondary_evidence: bool
    codon: int | None


# Exonuclease domain boundaries
EDM_CODON_START = 268
EDM_CODON_END = 471


def _load_pole_data() -> dict:
    """Load curated POLE pathogenic variant data."""
    data_path = Path(__file__).parent / "data" / "pole_pathogenic.json"
    with open(data_path) as f:
        return json.load(f)


_POLE_DATA: dict | None = None


def _get_pole_data() -> dict:
    global _POLE_DATA
    if _POLE_DATA is None:
        _POLE_DATA = _load_pole_data()
    return _POLE_DATA


def _parse_protein_change(hgvsp: str) -> tuple[str, int, str] | None:
    """Parse a protein change string like 'P286R' into (ref_aa, codon, alt_aa).

    Handles formats: P286R, p.P286R
    """
    s = hgvsp.strip()
    if s.startswith("p."):
        s = s[2:]

    # Standard missense: single letter AA + codon number + single letter AA
    m = re.match(r"^([A-Z])(\d+)([A-Z])$", s)
    if m:
        return m.group(1), int(m.group(2)), m.group(3)

    # Truncating patterns: X123*, X123fs, etc. — extract codon
    m = re.match(r"^([A-Z])(\d+)", s)
    if m:
        return m.group(1), int(m.group(2)), ""

    return None


def check_pole_variant(variant: Variant) -> PoleResult:
    """Check if a POLE variant is a pathogenic exonuclease domain mutation.

    Returns PoleResult with tier classification:
    - TIER1: established pathogenic hotspot
    - TIER2: additional confirmed pathogenic
    - VUS: in exonuclease domain but not in curated list
    - NOT_EDM: outside exonuclease domain
    """
    if variant.hugo_symbol != "POLE":
        return PoleResult(
            is_pathogenic=False,
            tier=PoleTier.NOT_EDM,
            variant_str="",
            needs_secondary_evidence=False,
            codon=None,
        )

    parsed = _parse_protein_change(variant.hgvsp_short)
    if parsed is None:
        return PoleResult(
            is_pathogenic=False,
            tier=PoleTier.NOT_EDM,
            variant_str=f"POLE {variant.hgvsp_short}",
            needs_secondary_evidence=False,
            codon=None,
        )

    ref_aa, codon, alt_aa = parsed
    variant_key = f"{ref_aa}{codon}{alt_aa}" if alt_aa else f"{ref_aa}{codon}"

    # Check if in exonuclease domain
    if not (EDM_CODON_START <= codon <= EDM_CODON_END):
        return PoleResult(
            is_pathogenic=False,
            tier=PoleTier.NOT_EDM,
            variant_str=f"POLE p.{variant_key}",
            needs_secondary_evidence=False,
            codon=codon,
        )

    data = _get_pole_data()

    # Check tier 1
    if variant_key in data["tier1"]["variants"]:
        return PoleResult(
            is_pathogenic=True,
            tier=PoleTier.TIER1,
            variant_str=f"POLE p.{variant_key}",
            needs_secondary_evidence=False,
            codon=codon,
        )

    # Check tier 2
    if variant_key in data["tier2"]["variants"]:
        return PoleResult(
            is_pathogenic=True,
            tier=PoleTier.TIER2,
            variant_str=f"POLE p.{variant_key}",
            needs_secondary_evidence=False,
            codon=codon,
        )

    # VUS in exonuclease domain
    return PoleResult(
        is_pathogenic=False,
        tier=PoleTier.VUS,
        variant_str=f"POLE p.{variant_key}",
        needs_secondary_evidence=True,
        codon=codon,
    )


def assess_pole(variants: list[Variant]) -> list[PoleResult]:
    """Assess all POLE variants in a sample.

    Returns list of PoleResult for all POLE variants found.
    """
    results = []
    for v in variants:
        if v.hugo_symbol == "POLE":
            results.append(check_pole_variant(v))
    return results


def get_best_pole_result(results: list[PoleResult]) -> PoleResult | None:
    """Get the highest-tier POLE result from a list.

    Priority: TIER1 > TIER2 > VUS > NOT_EDM
    """
    if not results:
        return None

    tier_priority = {
        PoleTier.TIER1: 0,
        PoleTier.TIER2: 1,
        PoleTier.VUS: 2,
        PoleTier.NOT_EDM: 3,
    }

    return min(results, key=lambda r: tier_priority[r.tier])
