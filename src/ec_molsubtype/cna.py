"""Copy number alteration (CNA) burden assessment."""

from __future__ import annotations

from typing import NamedTuple

import polars as pl

from .models import MolecularSubtype

# Log2 ratio threshold for calling a segment as altered
SEG_ALTERED_THRESHOLD = 0.2


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


def parse_seg_content(content: str) -> pl.DataFrame:
    """Parse CBS segmentation file content into a polars DataFrame.

    Accepts standard SEG format (tab-separated):
      ID  chrom  loc.start  loc.end  num.mark  seg.mean

    Strips comment lines starting with '#'.
    """
    lines = [line for line in content.splitlines(keepends=True) if not line.startswith("#")]
    if not lines:
        raise ValueError("Empty SEG content")
    cleaned = "".join(lines)
    df = pl.read_csv(cleaned.encode(), separator="\t",
                     schema_overrides={"chrom": pl.Utf8} if "chrom" in cleaned.split("\n")[0] else None)
    return df


def compute_fga_from_seg(content: str, sample_id: str | None = None) -> float:
    """Compute fraction genome altered (FGA) from SEG file content.

    FGA = sum of segment lengths where |seg.mean| > threshold / total covered bases.

    Args:
        content: SEG file content as string.
        sample_id: If the SEG file contains multiple samples, filter to this one.

    Returns:
        FGA as a float between 0 and 1.
    """
    df = parse_seg_content(content)

    # Normalize column names (handle different conventions)
    col_map = {}
    for col in df.columns:
        cl = col.lower().replace(".", "_").replace(" ", "_")
        if cl in ("loc_start", "start", "start_position"):
            col_map[col] = "loc_start"
        elif cl in ("loc_end", "end", "end_position"):
            col_map[col] = "loc_end"
        elif cl in ("seg_mean", "segmean", "log2ratio"):
            col_map[col] = "seg_mean"
        elif cl in ("id", "sample", "sample_id"):
            col_map[col] = "id"
    df = df.rename(col_map)

    # Filter to sample if multi-sample
    if "id" in df.columns:
        unique_ids = df["id"].unique().to_list()
        if len(unique_ids) == 1:
            pass  # single sample, use all rows
        elif sample_id and sample_id in unique_ids:
            df = df.filter(pl.col("id") == sample_id)
        elif sample_id:
            raise ValueError(
                f"Sample '{sample_id}' not found in SEG file. "
                f"Available: {', '.join(str(s) for s in unique_ids[:5])}"
            )
        else:
            # Multiple samples, no sample_id specified — use first
            df = df.filter(pl.col("id") == unique_ids[0])

    seg_length = (pl.col("loc_end") - pl.col("loc_start")).alias("seg_length")
    is_altered = (pl.col("seg_mean").abs() > SEG_ALTERED_THRESHOLD).alias("is_altered")

    result = df.select([
        seg_length,
        is_altered,
    ]).select([
        pl.col("seg_length").sum().alias("total"),
        (pl.col("seg_length") * pl.col("is_altered").cast(pl.Int64)).sum().alias("altered"),
    ])

    total = result["total"][0]
    altered = result["altered"][0]

    if total == 0:
        return 0.0
    return altered / total


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
