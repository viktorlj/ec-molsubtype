"""MSI proxy estimation from coding microsatellite frameshift analysis.

When direct MSI marker data (msi_pct) is unavailable — as in AACR GENIE public
data — this module estimates MSI status by counting frameshift indels at coding
microsatellite (cMS) indicator loci covered by the sequencing panel.

The approach:
1. Count frameshift indels at HIGH-SPECIFICITY cMS indicator genes
   (excluding ARID1A/PTEN which have high background rates in MSS EC)
2. Weight by specificity: each high-specificity gene hit contributes
   more to the proxy than non-specific genes
3. Compute the indel fraction (I-index) as supporting evidence
4. Produce a pseudo_msi_pct (0-100) that substitutes for direct MSI markers

Specificity calibration (from GENIE V18 MSK-IMPACT EC data):
  - JAK1:   59.6% in MSI-H, 0.1% in MSS → 448x enrichment
  - RNF43:  44.7% in MSI-H, 0.3% in MSS → 135x
  - MSH3:   25.7% in MSI-H, 0.0% in MSS → perfect specificity
  - MSH6:   17.2% in MSI-H, 0.0% in MSS → perfect specificity
  - KMT2C:  13.1% in MSI-H, 0.3% in MSS → 49x
  - AXIN2:  12.6% in MSI-H, 0.0% in MSS → perfect specificity
  - TCF7L2: 10.1% in MSI-H, 0.1% in MSS → 152x
  - APC:     8.0% in MSI-H, 0.2% in MSS → 40x
  - MRE11A:  6.2% in MSI-H, 0.1% in MSS → 93x
  - B2M:     4.6% in MSI-H, 0.0% in MSS → perfect specificity
  - TGFBR2:  2.3% in MSI-H, 0.1% in MSS → 17x
  - ARID1A: 91.5% in MSI-H, 11.4% in MSS → 8x  (NON-SPECIFIC)
  - PTEN:   62.8% in MSI-H, 8.0% in MSS  → 8x  (NON-SPECIFIC)

References:
- Kim JE et al. 2018 (I-index): DOI:10.1016/j.jmoldx.2018.09.005
- Cortes-Ciriano et al. 2017 (cMS targets): DOI:10.1038/ncomms15180
- Kawaguchi et al. 2009 (EC cMS targets): DOI:10.3892/ijo_00000411
- Ballhausen et al. 2020 (shared frameshift landscape): DOI:10.1038/s41467-020-18514-5
"""

from __future__ import annotations

from typing import NamedTuple

from .models import Variant

# Coding microsatellite indicator genes, stratified by MSI specificity.
#
# HIGH_SPECIFICITY: >99% specificity for MSI-H in EC. Background frameshift
# rate in MSS tumors is <0.3%. Each hit is strong evidence of MSI.
#
# LOW_SPECIFICITY: ARID1A and PTEN have high background frameshift rates
# (11.4% and 8.0%) in MSS endometrioid tumors — they are recurrently
# inactivated tumor suppressors regardless of MSI status. These provide
# supporting but not independent evidence.
#
# Key gap: ACVR2A (the #1 MSI indicator, 50-91% in MSI-H) is NOT on
# MSK-IMPACT. RPL22 (50% in EC) and BAX (23%) are also absent.

CMS_HIGH_SPECIFICITY: dict[str, dict] = {
    "JAK1":    {"repeat": "coding MS", "sensitivity": 0.596, "min_panel": 341},
    "RNF43":   {"repeat": "G(7)",      "sensitivity": 0.447, "min_panel": 341},
    "MSH3":    {"repeat": "A(8)",      "sensitivity": 0.257, "min_panel": 468},
    "MSH6":    {"repeat": "C(8)",      "sensitivity": 0.172, "min_panel": 341},
    "KMT2C":   {"repeat": "coding MS", "sensitivity": 0.131, "min_panel": 341},
    "AXIN2":   {"repeat": "G(7)",      "sensitivity": 0.126, "min_panel": 341},
    "TCF7L2":  {"repeat": "A(9)",      "sensitivity": 0.101, "min_panel": 468},
    "APC":     {"repeat": "multiple",  "sensitivity": 0.080, "min_panel": 341},
    "MRE11A":  {"repeat": "T(11)",     "sensitivity": 0.062, "min_panel": 341},
    "B2M":     {"repeat": "coding MS", "sensitivity": 0.046, "min_panel": 341},
    "TGFBR2":  {"repeat": "A(10)",     "sensitivity": 0.023, "min_panel": 341},
}

CMS_LOW_SPECIFICITY: dict[str, dict] = {
    "ARID1A":  {"repeat": "multiple", "sensitivity": 0.915, "min_panel": 341,
                "note": "11.4% background in MSS EC"},
    "PTEN":    {"repeat": "A(6)",     "sensitivity": 0.628, "min_panel": 341,
                "note": "8.0% background in MSS EC"},
}

# All indicator genes combined
CMS_INDICATOR_GENES = {**CMS_HIGH_SPECIFICITY, **CMS_LOW_SPECIFICITY}

# Frameshift variant classifications
FRAMESHIFT_CLASSES = {"Frame_Shift_Del", "Frame_Shift_Ins"}


class MsiProxyResult(NamedTuple):
    """MSI proxy estimation result."""

    pseudo_msi_pct: float               # Proxy MSI percentage (0-100)
    cms_high_spec_hit: list[str]        # High-specificity genes with frameshifts
    cms_low_spec_hit: list[str]         # Low-specificity genes with frameshifts
    cms_genes_total: int                # Total high-spec indicator genes on panel
    n_cms_frameshifts: int              # Total frameshift count at all indicator genes
    i_index: float                      # Indel fraction (all indels / all mutations)
    n_indels: int                       # Total indel count
    n_mutations: int                    # Total mutation count
    details: str


def get_available_genes(panel_id: str | None = None) -> tuple[set[str], set[str]]:
    """Get high- and low-specificity cMS genes available on a panel.

    Returns:
        Tuple of (high_spec_genes, low_spec_genes) available on the panel.
    """
    min_panel = 341  # default to smallest panel
    if panel_id is not None:
        panel_upper = panel_id.upper()
        for size in [505, 468, 410, 341]:
            if str(size) in panel_upper:
                min_panel = size
                break

    high = {g for g, v in CMS_HIGH_SPECIFICITY.items() if v["min_panel"] <= min_panel}
    low = {g for g, v in CMS_LOW_SPECIFICITY.items() if v["min_panel"] <= min_panel}
    return high, low


def compute_msi_proxy(
    variants: list[Variant],
    panel_id: str | None = None,
) -> MsiProxyResult:
    """Compute MSI proxy score from coding microsatellite frameshift analysis.

    The pseudo_msi_pct is computed from high-specificity genes only:

        base_score = (high_spec_genes_hit / total_high_spec_genes) * 100

    Each high-specificity gene hit contributes equally to the score.
    Low-specificity genes (ARID1A, PTEN) are reported but do NOT
    contribute to pseudo_msi_pct to avoid false positives.

    With 11 high-specificity genes on MSK-IMPACT468/505:
      1 gene hit → 9.1% (below default 20% MSI-H threshold → MSI-L)
      2 gene hits → 18.2% (borderline)
      3 gene hits → 27.3% (above threshold → MSI-H)

    The default msi_threshold of 20% in classify.py means 2-3+ genes
    are needed to call MSI-H, which is appropriately conservative.
    For panels with fewer indicator genes, the threshold may need
    adjustment.

    Args:
        variants: All variants for the sample.
        panel_id: Panel identifier for determining available cMS genes.

    Returns:
        MsiProxyResult with pseudo_msi_pct and supporting metrics.
    """
    high_genes, low_genes = get_available_genes(panel_id)
    all_genes = high_genes | low_genes
    n_high_total = len(high_genes)

    # Scan all variants
    n_mutations = 0
    n_indels = 0
    indel_classes_all = {"Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins"}
    non_coding = {"Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "IGR", "RNA"}

    high_hit: set[str] = set()
    low_hit: set[str] = set()
    n_cms_frameshifts = 0

    for v in variants:
        vc = v.variant_classification
        if vc in non_coding:
            continue
        n_mutations += 1
        if vc in indel_classes_all:
            n_indels += 1

        # Check frameshift at indicator gene
        if vc in FRAMESHIFT_CLASSES and v.hugo_symbol in all_genes:
            n_cms_frameshifts += 1
            if v.hugo_symbol in high_genes:
                high_hit.add(v.hugo_symbol)
            elif v.hugo_symbol in low_genes:
                low_hit.add(v.hugo_symbol)

    # Compute metrics
    high_list = sorted(high_hit)
    low_list = sorted(low_hit)
    n_high_hit = len(high_list)

    # pseudo_msi_pct based on high-specificity genes only
    pseudo_msi_pct = (n_high_hit / n_high_total * 100) if n_high_total > 0 else 0.0
    i_index = (n_indels / n_mutations) if n_mutations > 0 else 0.0

    # Build details
    parts = [
        f"{n_high_hit}/{n_high_total} high-specificity cMS genes",
    ]
    if high_list:
        parts.append(f"({', '.join(high_list)})")
    if low_list:
        parts.append(f"+ {len(low_list)} low-spec ({', '.join(low_list)})")
    parts.append(f"-> pseudo_msi={pseudo_msi_pct:.1f}%.")
    parts.append(f"I-index={i_index:.3f} ({n_indels}/{n_mutations}).")
    details = " ".join(parts)

    return MsiProxyResult(
        pseudo_msi_pct=pseudo_msi_pct,
        cms_high_spec_hit=high_list,
        cms_low_spec_hit=low_list,
        cms_genes_total=n_high_total,
        n_cms_frameshifts=n_cms_frameshifts,
        i_index=i_index,
        n_indels=n_indels,
        n_mutations=n_mutations,
        details=details,
    )
