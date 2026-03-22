#!/usr/bin/env python3
"""Extract endometrial cancer cases from AACR GENIE V18 for MSK-IMPACT panels.

Produces:
  data/genie/ec_msk_samples.tsv     — clinical metadata for EC samples on MSK panels
  data/genie/ec_msk_mutations.maf   — filtered MAF with mutations for those samples
  data/genie/ec_msk_tmb.tsv         — TMB bin + computed TMB for those samples
  data/genie/extraction_summary.txt — summary statistics

Usage:
  python scripts/extract_genie_ec.py /path/to/Genie_V18
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import polars as pl

# All OncoTree codes that map to endometrial cancer
EC_ONCOTREE_CODES = {
    "UCEC",  # Endometrial Carcinoma NOS
    "UEC",   # Endometrioid Carcinoma
    "USC",   # Serous Carcinoma / Papillary
    "UCS",   # Carcinosarcoma / MMMT
    "UCCC",  # Clear Cell Carcinoma
    "UMEC",  # Mixed Endometrial Carcinoma
    "UUC",   # Undifferentiated Uterine Carcinoma
    "UDDC",  # Dedifferentiated Carcinoma
    "UPDC",  # Poorly Differentiated Carcinoma
    "UNEC",  # Neuroendocrine Carcinoma
    "UASC",  # Adenosquamous Carcinoma
    "UMC",   # Mucinous Carcinoma
    "USCC",  # Squamous Cell Carcinoma
    "USMT",  # Smooth Muscle Tumor (borderline)
}

# Only MSK-IMPACT solid tumor panels (not heme, not ACCESS/ctDNA)
MSK_PANELS = {
    "MSK-IMPACT341",
    "MSK-IMPACT410",
    "MSK-IMPACT468",
    "MSK-IMPACT505",
}

# Panel sizes in Mb (computed from genomic_information.txt includeInPanel regions)
PANEL_SIZE_MB = {
    "MSK-IMPACT341": 0.932,
    "MSK-IMPACT410": 1.056,
    "MSK-IMPACT468": 1.171,
    "MSK-IMPACT505": 1.261,
}

# Variant classifications to exclude from TMB computation (silent/non-coding)
SILENT_CLASSIFICATIONS = {
    "Silent",
    "Intron",
    "3'UTR",
    "5'UTR",
    "3'Flank",
    "5'Flank",
    "IGR",
    "RNA",
}


def main(genie_dir: Path, output_dir: Path) -> None:
    t0 = time.time()

    # --- 1. Load clinical sample data ---
    print("Loading clinical sample data...")
    clinical = pl.read_csv(
        genie_dir / "data_clinical_sample.txt",
        separator="\t",
        comment_prefix="#",
        infer_schema_length=0,
    )

    # Filter: endometrial cancer + MSK-IMPACT panels
    ec_samples = clinical.filter(
        pl.col("ONCOTREE_CODE").is_in(EC_ONCOTREE_CODES)
        & pl.col("SEQ_ASSAY_ID").is_in(MSK_PANELS)
    )
    sample_ids = set(ec_samples["SAMPLE_ID"].to_list())
    print(f"  Found {len(sample_ids)} EC samples on MSK-IMPACT panels")
    print("  OncoTree distribution:")
    for row in ec_samples.group_by("ONCOTREE_CODE").len().sort("len", descending=True).iter_rows():
        print(f"    {row[0]}: {row[1]}")
    print("  Panel distribution:")
    for row in ec_samples.group_by("SEQ_ASSAY_ID").len().sort("len", descending=True).iter_rows():
        print(f"    {row[0]}: {row[1]}")

    # Save clinical
    ec_samples.write_csv(output_dir / "ec_msk_samples.tsv", separator="\t")
    print(f"  Wrote {output_dir / 'ec_msk_samples.tsv'}")

    # --- 2. Load TMB bin data ---
    # Note: GENIE public release anonymizes raw TMB values (all near-zero).
    # Only tmb_bin is informative. We compute TMB from mutation counts below.
    print("Loading TMB bin data...")
    tmb = pl.read_csv(
        genie_dir / "tmb_18.0-public.tsv",
        separator="\t",
    )
    ec_tmb = tmb.filter(pl.col("SAMPLE_ID").is_in(sample_ids)).select(["SAMPLE_ID", "tmb_bin"])
    print(f"  Matched {len(ec_tmb)} / {len(sample_ids)} samples with TMB bin data")
    print("  TMB bin distribution:")
    for row in ec_tmb.group_by("tmb_bin").len().sort("len", descending=True).iter_rows():
        print(f"    {row[0]}: {row[1]}")

    # --- 3. Filter mutations (large file — stream with schema overrides) ---
    print("Filtering mutations (this may take a minute)...")
    mutations = pl.scan_csv(
        genie_dir / "data_mutations_extended.txt",
        separator="\t",
        comment_prefix="#",
        infer_schema_length=10000,
        null_values=["", "NA", "nan"],
        schema_overrides={
            "Chromosome": pl.Utf8,
            "FILTER": pl.Utf8,
            "SWISSPROT": pl.Utf8,
            "Codons": pl.Utf8,
        },
    ).filter(
        pl.col("Tumor_Sample_Barcode").is_in(sample_ids)
    ).collect()

    n_samples_with_muts = mutations["Tumor_Sample_Barcode"].n_unique()
    print(f"  Extracted {len(mutations)} mutations for {n_samples_with_muts} samples")

    # Key gene counts
    key_genes = ["POLE", "TP53", "MLH1", "MSH2", "MSH6", "PMS2", "PTEN", "PIK3CA", "ARID1A", "KRAS"]
    gene_counts = (
        mutations.filter(pl.col("Hugo_Symbol").is_in(key_genes))
        .group_by("Hugo_Symbol")
        .len()
        .sort("len", descending=True)
    )
    print("  Key gene mutation counts:")
    for row in gene_counts.iter_rows():
        print(f"    {row[0]}: {row[1]}")

    # Variant classifications
    print("  Variant classification distribution:")
    vc_counts = mutations.group_by("Variant_Classification").len().sort("len", descending=True)
    for row in vc_counts.head(10).iter_rows():
        print(f"    {row[0]}: {row[1]}")

    # Write MAF
    mutations.write_csv(output_dir / "ec_msk_mutations.maf", separator="\t")
    print(f"  Wrote {output_dir / 'ec_msk_mutations.maf'}")

    # --- 4. Compute TMB from mutation counts + panel sizes ---
    print("Computing TMB from mutation counts...")
    # Build sample -> panel mapping
    sample_panel = ec_samples.select(["SAMPLE_ID", "SEQ_ASSAY_ID"])

    # Count non-silent mutations per sample
    non_silent = mutations.filter(
        ~pl.col("Variant_Classification").is_in(SILENT_CLASSIFICATIONS)
    )
    mut_counts = (
        non_silent.group_by("Tumor_Sample_Barcode")
        .len()
        .rename({"Tumor_Sample_Barcode": "SAMPLE_ID", "len": "mutation_count"})
    )

    # Join with panel info and compute TMB
    tmb_computed = (
        sample_panel.join(mut_counts, on="SAMPLE_ID", how="left")
        .with_columns(pl.col("mutation_count").fill_null(0))
        .with_columns(
            pl.col("SEQ_ASSAY_ID")
            .replace_strict(PANEL_SIZE_MB, default=None)
            .alias("panel_size_mb")
        )
        .with_columns(
            (pl.col("mutation_count") / pl.col("panel_size_mb")).alias("tmb_computed")
        )
        .join(ec_tmb, on="SAMPLE_ID", how="left")
    )

    # Save TMB
    tmb_computed.write_csv(output_dir / "ec_msk_tmb.tsv", separator="\t")
    print(f"  Wrote {output_dir / 'ec_msk_tmb.tsv'}")

    # TMB summary
    tmb_vals = tmb_computed["tmb_computed"].drop_nulls()
    print(f"  Computed TMB: median={tmb_vals.median():.1f}, mean={tmb_vals.mean():.1f}, "
          f"max={tmb_vals.max():.1f}")
    high_tmb = tmb_vals.filter(tmb_vals >= 10).len()
    very_high_tmb = tmb_vals.filter(tmb_vals >= 100).len()
    print(f"  TMB >= 10: {high_tmb} samples ({100*high_tmb/len(tmb_computed):.1f}%)")
    print(f"  TMB >= 100: {very_high_tmb} samples ({100*very_high_tmb/len(tmb_computed):.1f}%)")

    # --- 5. Summary ---
    elapsed = time.time() - t0
    summary_lines = [
        "GENIE V18 -- Endometrial Cancer MSK-IMPACT Extraction",
        f"Date: {time.strftime('%Y-%m-%d %H:%M')}",
        f"Extraction time: {elapsed:.0f}s",
        "",
        f"Total EC samples on MSK-IMPACT: {len(sample_ids)}",
        f"Total mutations extracted: {len(mutations)}",
        f"Non-silent mutations: {len(non_silent)}",
        f"Samples with mutations: {n_samples_with_muts}",
        f"Samples with TMB bin: {len(ec_tmb)}",
        "",
        "OncoTree code distribution:",
    ]
    for row in ec_samples.group_by("ONCOTREE_CODE").len().sort("len", descending=True).iter_rows():
        summary_lines.append(f"  {row[0]}: {row[1]}")
    summary_lines.extend(["", "Panel distribution:"])
    for row in ec_samples.group_by("SEQ_ASSAY_ID").len().sort("len", descending=True).iter_rows():
        summary_lines.append(f"  {row[0]}: {row[1]}")
    summary_lines.extend(["", "Key gene mutations:"])
    for row in gene_counts.iter_rows():
        summary_lines.append(f"  {row[0]}: {row[1]}")
    summary_lines.extend([
        "",
        f"Computed TMB (mut/Mb): median={tmb_vals.median():.1f}, "
        f"mean={tmb_vals.mean():.1f}, max={tmb_vals.max():.1f}",
        f"TMB >= 10 mut/Mb: {high_tmb} ({100*high_tmb/len(tmb_computed):.1f}%)",
        f"TMB >= 100 mut/Mb: {very_high_tmb} ({100*very_high_tmb/len(tmb_computed):.1f}%)",
        "",
        "GENIE TMB bin distribution:",
    ])
    for row in ec_tmb.group_by("tmb_bin").len().sort("len", descending=True).iter_rows():
        summary_lines.append(f"  {row[0]}: {row[1]}")
    summary_lines.extend([
        "",
        "Note: Raw TMB values in GENIE public release are anonymized.",
        "TMB was computed from non-silent mutation count / panel size (Mb).",
        "Panel sizes from genomic_information.txt includeInPanel regions.",
    ])

    summary_text = "\n".join(summary_lines) + "\n"
    (output_dir / "extraction_summary.txt").write_text(summary_text)
    print(f"\n{'='*60}")
    print(summary_text)
    print(f"Output files in: {output_dir}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <genie_v18_dir> [output_dir]")
        sys.exit(1)

    genie_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("data/genie")
    output_dir.mkdir(parents=True, exist_ok=True)

    if not (genie_dir / "data_clinical_sample.txt").exists():
        print(f"Error: {genie_dir / 'data_clinical_sample.txt'} not found")
        sys.exit(1)

    main(genie_dir, output_dir)
