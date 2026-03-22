#!/usr/bin/env python3
"""Extract CNA segmentation data for EC MSK-IMPACT samples from GENIE V18.

Computes per-sample CIN (chromosomal instability) metrics:
  - fraction_genome_altered (FGA): fraction of genome with |seg.mean| > threshold
  - weighted_genome_instability_index (wGII): weighted by segment size
  - n_segments: total segments (proxy for breakpoint count)
  - n_amplified / n_deleted segments
  - arm-level CNA counts

Produces:
  data/genie/ec_msk_cna_segments.tsv  — filtered segmentation data
  data/genie/ec_msk_cin_metrics.tsv   — per-sample CIN metrics
  data/genie/ec_msk_tp53_cin.tsv      — merged TP53 status + CIN for analysis

Usage:
  python scripts/extract_genie_cna.py /path/to/Genie_V18
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import polars as pl

# hg19 chromosome sizes (autosomes + X)
CHROM_SIZES_HG19 = {
    "1": 249250621, "2": 243199373, "3": 198022430, "4": 191154276,
    "5": 180915260, "6": 171115067, "7": 159138663, "8": 146364022,
    "9": 141213431, "10": 135534747, "11": 135006516, "12": 133851895,
    "13": 115169878, "14": 107349540, "15": 102531392, "16": 90354753,
    "17": 81195210, "18": 78077248, "19": 59128983, "20": 63025520,
    "21": 48129895, "22": 51304566, "X": 155270560,
}
GENOME_SIZE = sum(CHROM_SIZES_HG19.values())

# hg19 centromere positions (approximate, for arm-level calls)
CENTROMERES_HG19 = {
    "1": 125000000, "2": 93300000, "3": 91000000, "4": 50400000,
    "5": 48400000, "6": 61000000, "7": 59900000, "8": 45600000,
    "9": 49000000, "10": 40200000, "11": 53700000, "12": 35800000,
    "13": 17900000, "14": 17600000, "15": 19000000, "16": 36600000,
    "17": 24000000, "18": 17200000, "19": 26500000, "20": 27500000,
    "21": 13200000, "22": 14700000, "X": 60600000,
}

# Log2 ratio thresholds for calling gains/losses
AMP_THRESHOLD = 0.2    # log2 ratio > 0.2 → gain
DEL_THRESHOLD = -0.2   # log2 ratio < -0.2 → loss

# TP53 hotspot residues (for classification)
TP53_HOTSPOTS = {"R175H", "G245S", "R248W", "R248Q", "R249S", "R273H", "R273C", "R282W"}
TP53_TRUNCATING = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"}


def compute_cin_metrics(segments: pl.DataFrame) -> pl.DataFrame:
    """Compute per-sample CIN metrics from segmentation data."""

    # Add segment length and alteration status
    seg = segments.with_columns([
        (pl.col("loc.end") - pl.col("loc.start")).alias("seg_length"),
        (pl.col("seg.mean") > AMP_THRESHOLD).alias("is_amp"),
        (pl.col("seg.mean") < DEL_THRESHOLD).alias("is_del"),
        (pl.col("seg.mean").abs() > abs(DEL_THRESHOLD)).alias("is_altered"),
    ])

    # Per-sample metrics
    metrics = seg.group_by("ID").agg([
        # Total segments (breakpoint complexity)
        pl.len().alias("n_segments"),

        # Segment counts by type
        pl.col("is_amp").sum().alias("n_amplified_segments"),
        pl.col("is_del").sum().alias("n_deleted_segments"),

        # Total genome covered
        pl.col("seg_length").sum().alias("total_covered_bp"),

        # Altered bases
        (pl.col("seg_length") * pl.col("is_altered").cast(pl.Int64)).sum().alias("altered_bp"),
        (pl.col("seg_length") * pl.col("is_amp").cast(pl.Int64)).sum().alias("amplified_bp"),
        (pl.col("seg_length") * pl.col("is_del").cast(pl.Int64)).sum().alias("deleted_bp"),

        # Weighted mean absolute log2 ratio (wGII proxy)
        (pl.col("seg.mean").abs() * pl.col("seg_length")).sum().alias("weighted_abs_log2r_sum"),
    ]).with_columns([
        # FGA: fraction of covered genome that is altered
        (pl.col("altered_bp") / pl.col("total_covered_bp")).alias("fga"),

        # Fraction amplified / deleted
        (pl.col("amplified_bp") / pl.col("total_covered_bp")).alias("fraction_amplified"),
        (pl.col("deleted_bp") / pl.col("total_covered_bp")).alias("fraction_deleted"),

        # wGII: weighted mean |log2R| across genome
        (pl.col("weighted_abs_log2r_sum") / pl.col("total_covered_bp")).alias("wgii"),
    ]).drop(["weighted_abs_log2r_sum"])

    return metrics


def classify_tp53_status(mutations: pl.DataFrame) -> pl.DataFrame:
    """Classify each sample's TP53 status from mutation data."""

    tp53 = mutations.filter(pl.col("Hugo_Symbol") == "TP53")

    if len(tp53) == 0:
        return pl.DataFrame({
            "SAMPLE_ID": [],
            "tp53_status": [],
            "tp53_variant": [],
            "tp53_variant_class": [],
            "tp53_vaf": [],
        }, schema={
            "SAMPLE_ID": pl.Utf8,
            "tp53_status": pl.Utf8,
            "tp53_variant": pl.Utf8,
            "tp53_variant_class": pl.Utf8,
            "tp53_vaf": pl.Float64,
        })

    # Compute VAF
    tp53 = tp53.with_columns([
        pl.col("HGVSp_Short").fill_null("").alias("protein_change"),
        pl.when(
            (pl.col("t_alt_count").is_not_null()) & (pl.col("t_ref_count").is_not_null())
            & ((pl.col("t_alt_count") + pl.col("t_ref_count")) > 0)
        ).then(
            pl.col("t_alt_count") / (pl.col("t_alt_count") + pl.col("t_ref_count"))
        ).otherwise(None).alias("vaf"),
    ])

    # Classify each TP53 variant
    tp53 = tp53.with_columns([
        pl.when(
            pl.col("protein_change").str.replace("p.", "").is_in(TP53_HOTSPOTS)
        ).then(pl.lit("hotspot"))
        .when(
            pl.col("Variant_Classification").is_in(TP53_TRUNCATING)
        ).then(pl.lit("truncating"))
        .when(
            pl.col("Variant_Classification") == "Missense_Mutation"
        ).then(pl.lit("missense"))
        .otherwise(pl.lit("other"))
        .alias("tp53_category"),
    ])

    # Pick the "worst" TP53 variant per sample (hotspot > truncating > missense)
    priority = {"hotspot": 0, "truncating": 1, "missense": 2, "other": 3}
    tp53 = tp53.with_columns(
        pl.col("tp53_category").replace_strict(priority, default=4).alias("priority")
    )

    best_tp53 = (
        tp53.sort(["Tumor_Sample_Barcode", "priority"])
        .group_by("Tumor_Sample_Barcode")
        .first()
        .select([
            pl.col("Tumor_Sample_Barcode").alias("SAMPLE_ID"),
            pl.lit("mutated").alias("tp53_status"),
            pl.col("protein_change").alias("tp53_variant"),
            pl.col("tp53_category").alias("tp53_variant_class"),
            pl.col("vaf").alias("tp53_vaf"),
        ])
    )

    return best_tp53


def main(genie_dir: Path, output_dir: Path) -> None:
    t0 = time.time()

    # Load sample list
    print("Loading EC sample list...")
    ec_samples = pl.read_csv(output_dir / "ec_msk_samples.tsv", separator="\t",
                             infer_schema_length=0)  # all strings (age has ">89")
    sample_ids = set(ec_samples["SAMPLE_ID"].to_list())
    print(f"  {len(sample_ids)} samples")

    # Extract CNA segments
    print("Filtering CNA segments (424 MB file)...")
    segments = pl.scan_csv(
        genie_dir / "data_cna_hg19.seg",
        separator="\t",
        schema_overrides={"chrom": pl.Utf8},
    ).filter(
        pl.col("ID").is_in(sample_ids)
    ).collect()

    n_samples_cna = segments["ID"].n_unique()
    print(f"  Extracted {len(segments)} segments for {n_samples_cna} / {len(sample_ids)} samples")

    # Save segments
    segments.write_csv(output_dir / "ec_msk_cna_segments.tsv", separator="\t")
    print(f"  Wrote {output_dir / 'ec_msk_cna_segments.tsv'}")

    # Compute CIN metrics
    print("Computing CIN metrics...")
    cin = compute_cin_metrics(segments)

    # Summary stats
    fga = cin["fga"]
    print(f"  FGA: median={fga.median():.3f}, mean={fga.mean():.3f}, max={fga.max():.3f}")
    wgii = cin["wgii"]
    print(f"  wGII: median={wgii.median():.4f}, mean={wgii.mean():.4f}, max={wgii.max():.4f}")
    n_segs = cin["n_segments"]
    print(f"  Segments per sample: median={n_segs.median():.0f}, mean={n_segs.mean():.0f}, max={n_segs.max()}")

    # Classify TP53 status
    print("Classifying TP53 status from mutations...")
    mutations = pl.read_csv(output_dir / "ec_msk_mutations.maf", separator="\t",
                            infer_schema_length=10000, null_values=["", "NA", "nan"])
    tp53_status = classify_tp53_status(mutations)
    print(f"  TP53 mutated: {len(tp53_status)} samples")
    if len(tp53_status) > 0:
        print("  TP53 variant classes:")
        for row in tp53_status.group_by("tp53_variant_class").len().sort("len", descending=True).iter_rows():
            print(f"    {row[0]}: {row[1]}")

    # Merge: all samples with CIN metrics + TP53 status + OncoTree + TMB
    print("Merging TP53 + CIN + clinical data...")
    tmb = pl.read_csv(output_dir / "ec_msk_tmb.tsv", separator="\t")

    merged = (
        cin.rename({"ID": "SAMPLE_ID"})
        .join(tp53_status, on="SAMPLE_ID", how="left")
        .with_columns(
            pl.col("tp53_status").fill_null("wild_type")
        )
        .join(
            ec_samples.select(["SAMPLE_ID", "ONCOTREE_CODE", "SEQ_ASSAY_ID"]),
            on="SAMPLE_ID", how="left"
        )
        .join(
            tmb.select(["SAMPLE_ID", "tmb_computed", "tmb_bin"]),
            on="SAMPLE_ID", how="left"
        )
    )

    # Save CIN metrics
    cin.rename({"ID": "SAMPLE_ID"}).write_csv(output_dir / "ec_msk_cin_metrics.tsv", separator="\t")
    print(f"  Wrote {output_dir / 'ec_msk_cin_metrics.tsv'}")

    # Save merged TP53 + CIN
    merged.write_csv(output_dir / "ec_msk_tp53_cin.tsv", separator="\t")
    print(f"  Wrote {output_dir / 'ec_msk_tp53_cin.tsv'}")

    # Summary statistics: TP53 mutant vs wild-type CIN
    print(f"\n{'='*60}")
    print("TP53 status vs CIN metrics:")
    for status in ["mutated", "wild_type"]:
        subset = merged.filter(pl.col("tp53_status") == status)
        n = len(subset)
        if n == 0:
            continue
        fga_s = subset["fga"]
        wgii_s = subset["wgii"]
        print(f"\n  {status.upper()} (n={n}):")
        print(f"    FGA:  median={fga_s.median():.3f}, mean={fga_s.mean():.3f}, IQR=[{fga_s.quantile(0.25):.3f}-{fga_s.quantile(0.75):.3f}]")
        print(f"    wGII: median={wgii_s.median():.4f}, mean={wgii_s.mean():.4f}")

    # By OncoTree code
    print(f"\n{'='*60}")
    print("CIN by OncoTree code (FGA):")
    for row in (merged.group_by("ONCOTREE_CODE").agg([
        pl.len().alias("n"),
        pl.col("fga").median().alias("median_fga"),
        pl.col("fga").mean().alias("mean_fga"),
        (pl.col("tp53_status") == "mutated").sum().alias("n_tp53_mut"),
    ]).sort("n", descending=True)).iter_rows():
        code, n, med_fga, mean_fga, n_tp53 = row
        tp53_pct = 100 * n_tp53 / n if n > 0 else 0
        print(f"  {code:6s} n={n:4d}  FGA median={med_fga:.3f}  mean={mean_fga:.3f}  TP53mut={n_tp53:3d} ({tp53_pct:.0f}%)")

    # By TP53 variant class
    if len(tp53_status) > 0:
        print(f"\n{'='*60}")
        print("CIN by TP53 variant class:")
        tp53_cin = merged.filter(pl.col("tp53_status") == "mutated")
        for row in (tp53_cin.group_by("tp53_variant_class").agg([
            pl.len().alias("n"),
            pl.col("fga").median().alias("median_fga"),
            pl.col("fga").mean().alias("mean_fga"),
        ]).sort("n", descending=True)).iter_rows():
            cls, n, med, mean = row
            print(f"  {cls:12s} n={n:4d}  FGA median={med:.3f}  mean={mean:.3f}")

    elapsed = time.time() - t0
    print(f"\nDone in {elapsed:.0f}s")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <genie_v18_dir> [output_dir]")
        sys.exit(1)

    genie_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("data/genie")

    if not (output_dir / "ec_msk_samples.tsv").exists():
        print(f"Error: run extract_genie_ec.py first — {output_dir / 'ec_msk_samples.tsv'} not found")
        sys.exit(1)

    main(genie_dir, output_dir)
