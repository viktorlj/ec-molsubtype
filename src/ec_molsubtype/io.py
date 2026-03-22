"""MAF/VCF parsing for ec-molsubtype."""

from __future__ import annotations

import json
from pathlib import Path

import polars as pl

from .models import SampleInput, SampleMetadata, Variant

# Minimum required MAF columns
REQUIRED_MAF_COLUMNS = {
    "Hugo_Symbol",
    "Variant_Classification",
    "HGVSp_Short",
}

OPTIONAL_MAF_COLUMNS = {
    "Chromosome",
    "Start_Position",
    "End_Position",
    "Reference_Allele",
    "Tumor_Seq_Allele2",
    "t_alt_count",
    "t_ref_count",
    "ClinVar_Classification",
    "SIFT",
    "PolyPhen",
    "Tumor_Sample_Barcode",
}


def parse_maf(path: str | Path) -> pl.DataFrame:
    """Parse a MAF file into a polars DataFrame.

    Handles comment lines starting with '#'.
    """
    path = Path(path)

    # Read skipping comment lines
    lines = []
    with open(path) as f:
        for line in f:
            if not line.startswith("#"):
                lines.append(line)

    if not lines:
        raise ValueError(f"Empty MAF file: {path}")

    # Write to temp string and read with polars
    content = "".join(lines)
    df = pl.read_csv(content.encode(), separator="\t", infer_schema_length=0)

    # Validate required columns
    missing = REQUIRED_MAF_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(f"MAF file missing required columns: {missing}")

    return df


def maf_row_to_variant(row: dict) -> Variant:
    """Convert a MAF row (as dict) to a Variant model."""
    return Variant(
        hugo_symbol=row.get("Hugo_Symbol", ""),
        chromosome=row.get("Chromosome", ""),
        start_position=int(row.get("Start_Position", 0) or 0),
        end_position=int(row.get("End_Position", 0) or 0),
        reference_allele=row.get("Reference_Allele", ""),
        tumor_seq_allele2=row.get("Tumor_Seq_Allele2", ""),
        variant_classification=row.get("Variant_Classification", ""),
        hgvsp_short=row.get("HGVSp_Short", ""),
        t_alt_count=int(row.get("t_alt_count", 0) or 0),
        t_ref_count=int(row.get("t_ref_count", 0) or 0),
        clinvar_classification=row.get("ClinVar_Classification", ""),
        sift=row.get("SIFT", ""),
        polyphen=row.get("PolyPhen", ""),
    )


def load_maf_variants(path: str | Path) -> list[Variant]:
    """Load all variants from a MAF file."""
    df = parse_maf(path)
    variants = []
    for row in df.iter_rows(named=True):
        variants.append(maf_row_to_variant(row))
    return variants


def load_sample_metadata(path: str | Path) -> SampleMetadata:
    """Load sample metadata from a JSON sidecar file."""
    path = Path(path)
    with open(path) as f:
        data = json.load(f)

    return SampleMetadata(**data)


def load_sample_metadata_tsv(path: str | Path) -> list[SampleMetadata]:
    """Load sample metadata from a TSV file (for batch processing)."""
    df = pl.read_csv(str(path), separator="\t", infer_schema_length=0)
    results = []

    for row in df.iter_rows(named=True):
        meta = SampleMetadata(
            sample_id=row.get("sample_id", ""),
            tmb=float(row["tmb"]) if row.get("tmb") else None,
            panel_size_mb=float(row["panel_size_mb"]) if row.get("panel_size_mb") else None,
            msi_pct=float(row["msi_pct"]) if row.get("msi_pct") else None,
            msi_status=row.get("msi_status") or None,
            fraction_genome_altered=(
                float(row["fraction_genome_altered"])
                if row.get("fraction_genome_altered")
                else None
            ),
            signature_weights=(
                json.loads(row["signature_weights"])
                if row.get("signature_weights")
                else None
            ),
        )
        results.append(meta)

    return results


def load_sample(
    maf_path: str | Path,
    metadata_path: str | Path | None = None,
    sample_id: str | None = None,
) -> SampleInput:
    """Load a complete sample from MAF + optional metadata.

    If metadata_path is not provided, creates minimal metadata from the MAF.
    """
    variants = load_maf_variants(maf_path)

    if metadata_path:
        metadata = load_sample_metadata(metadata_path)
    else:
        # Infer sample_id from MAF if possible
        if not sample_id:
            # Try to get from first variant's Tumor_Sample_Barcode
            df = parse_maf(maf_path)
            if "Tumor_Sample_Barcode" in df.columns:
                sample_id = df["Tumor_Sample_Barcode"][0]
            else:
                sample_id = Path(maf_path).stem

        metadata = SampleMetadata(sample_id=sample_id)

    return SampleInput(metadata=metadata, variants=variants)
