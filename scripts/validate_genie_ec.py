#!/usr/bin/env python3
"""Run ec-molsubtype classifier on GENIE V18 MSK-IMPACT EC samples.

Loads the extracted GENIE data, computes MSI proxy from cMS frameshifts,
and runs the hierarchical WHO 2020 classifier on all samples.

Produces:
  data/genie/ec_msk_classifications.tsv  — per-sample classification results
  data/genie/validation_summary.txt      — summary statistics and sanity checks

Usage:
  python scripts/validate_genie_ec.py [data/genie]
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import polars as pl

# Add src to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from ec_molsubtype.classify import classify_sample
from ec_molsubtype.models import SampleInput, SampleMetadata, Variant
from ec_molsubtype.msi import compute_msi_proxy

# Panel sizes in Mb
PANEL_SIZE_MB = {
    "MSK-IMPACT341": 0.932,
    "MSK-IMPACT410": 1.056,
    "MSK-IMPACT468": 1.171,
    "MSK-IMPACT505": 1.261,
}


def build_sample_input(
    sample_id: str,
    panel_id: str,
    variants_df: pl.DataFrame,
    tmb_value: float | None,
    fga_value: float | None,
) -> SampleInput:
    """Convert GENIE data to SampleInput for the classifier."""
    # Convert polars rows to Variant objects
    variant_list: list[Variant] = []
    for row in variants_df.iter_rows(named=True):
        v = Variant(
            hugo_symbol=row.get("Hugo_Symbol", "") or "",
            chromosome=str(row.get("Chromosome", "") or ""),
            start_position=int(row.get("Start_Position", 0) or 0),
            end_position=int(row.get("End_Position", 0) or 0),
            reference_allele=str(row.get("Reference_Allele", "") or ""),
            tumor_seq_allele2=str(row.get("Tumor_Seq_Allele2", "") or ""),
            variant_classification=row.get("Variant_Classification", "") or "",
            hgvsp_short=row.get("HGVSp_Short", "") or "",
            t_alt_count=int(row.get("t_alt_count", 0) or 0),
            t_ref_count=int(row.get("t_ref_count", 0) or 0),
        )
        variant_list.append(v)

    # Compute MSI proxy
    msi_proxy = compute_msi_proxy(variant_list, panel_id=panel_id)

    metadata = SampleMetadata(
        sample_id=sample_id,
        tmb=tmb_value,
        panel_size_mb=PANEL_SIZE_MB.get(panel_id),
        msi_pct=msi_proxy.pseudo_msi_pct,
        fraction_genome_altered=fga_value,
    )

    return SampleInput(metadata=metadata, variants=variant_list), msi_proxy


def main(data_dir: Path, msi_threshold: float = 9.0) -> None:
    t0 = time.time()

    # Load data
    print("Loading extracted GENIE data...")
    samples = pl.read_csv(data_dir / "ec_msk_samples.tsv", separator="\t",
                          infer_schema_length=0)
    maf = pl.read_csv(data_dir / "ec_msk_mutations.maf", separator="\t",
                       infer_schema_length=10000, null_values=["", "NA", "nan"])
    tmb_df = pl.read_csv(data_dir / "ec_msk_tmb.tsv", separator="\t")
    cin_file = data_dir / "ec_msk_cin_metrics.tsv"
    cin_df = pl.read_csv(cin_file, separator="\t") if cin_file.exists() else None

    # Build lookup tables
    tmb_lookup = dict(zip(
        tmb_df["SAMPLE_ID"].to_list(),
        tmb_df["tmb_computed"].to_list(),
    ))
    fga_lookup = {}
    if cin_df is not None:
        fga_lookup = dict(zip(
            cin_df["SAMPLE_ID"].to_list(),
            cin_df["fga"].to_list(),
        ))

    # Group mutations by sample
    print("Grouping mutations by sample...")
    sample_mutations = {}
    for sid in maf["Tumor_Sample_Barcode"].unique().to_list():
        sample_mutations[sid] = maf.filter(pl.col("Tumor_Sample_Barcode") == sid)

    # Classify all samples
    n_total = len(samples)
    print(f"Classifying {n_total} samples...")

    results = []
    errors = []
    for i, row in enumerate(samples.iter_rows(named=True)):
        sid = row["SAMPLE_ID"]
        panel_id = row["SEQ_ASSAY_ID"]
        oncotree = row["ONCOTREE_CODE"]

        if (i + 1) % 500 == 0:
            print(f"  {i+1}/{n_total}...")

        variants_df = sample_mutations.get(sid, pl.DataFrame())
        tmb_val = tmb_lookup.get(sid)
        fga_val = fga_lookup.get(sid)

        try:
            sample_input, msi_proxy = build_sample_input(
                sid, panel_id, variants_df, tmb_val, fga_val,
            )
            result = classify_sample(sample_input, msi_threshold=msi_threshold)

            results.append({
                "SAMPLE_ID": sid,
                "ONCOTREE_CODE": oncotree,
                "SEQ_ASSAY_ID": panel_id,
                "subtype": result.primary_subtype.value,
                "confidence": result.confidence.value,
                "is_multiple_classifier": result.multiple_classifier.is_multiple,
                "secondary_features": ",".join(result.multiple_classifier.secondary_features),
                "n_flags": len(result.flags),
                "flags": "; ".join(result.flags) if result.flags else "",
                "tmb": tmb_val,
                "fga": fga_val,
                "pseudo_msi_pct": msi_proxy.pseudo_msi_pct,
                "cms_genes_hit": len(msi_proxy.cms_high_spec_hit),
                "cms_genes_list": ",".join(msi_proxy.cms_high_spec_hit),
                "cms_low_spec": ",".join(msi_proxy.cms_low_spec_hit) if msi_proxy.cms_low_spec_hit else "",
                "i_index": msi_proxy.i_index,
                "n_mutations": msi_proxy.n_mutations,
                "n_indels": msi_proxy.n_indels,
            })

        except Exception as e:
            errors.append({"SAMPLE_ID": sid, "error": str(e)})

    print(f"  Classified {len(results)} samples ({len(errors)} errors)")

    # Convert to dataframe and save
    results_df = pl.DataFrame(results)
    results_df.write_csv(data_dir / "ec_msk_classifications.tsv", separator="\t")
    print(f"  Wrote {data_dir / 'ec_msk_classifications.tsv'}")

    # === Summary statistics ===
    elapsed = time.time() - t0
    lines = [
        "GENIE V18 EC MSK-IMPACT — Classification Validation",
        f"Date: {time.strftime('%Y-%m-%d %H:%M')}",
        f"Runtime: {elapsed:.0f}s",
        f"Samples classified: {len(results)} / {n_total}",
        f"Errors: {len(errors)}",
        f"MSI threshold: {msi_threshold}%",
        "",
    ]

    # Subtype distribution
    lines.append("=" * 60)
    lines.append("SUBTYPE DISTRIBUTION")
    lines.append("=" * 60)
    for row in results_df.group_by("subtype").len().sort("len", descending=True).iter_rows():
        pct = 100 * row[1] / len(results)
        lines.append(f"  {row[0]:10s}: {row[1]:5d} ({pct:5.1f}%)")

    # Expected from TCGA: POLEmut ~7%, MMRd ~28%, p53abn ~26%, NSMP ~39%
    lines.extend(["", "Expected (TCGA 2013): POLEmut ~7%, MMRd ~28%, p53abn ~26%, NSMP ~39%"])

    # Subtype by OncoTree code
    lines.extend(["", "=" * 60, "SUBTYPE BY HISTOLOGY (ONCOTREE CODE)", "=" * 60])
    for oncotree in ["UEC", "USC", "UCS", "UCEC", "UMEC", "UCCC"]:
        subset = results_df.filter(pl.col("ONCOTREE_CODE") == oncotree)
        if len(subset) == 0:
            continue
        lines.append(f"\n  {oncotree} (n={len(subset)}):")
        for row in subset.group_by("subtype").len().sort("len", descending=True).iter_rows():
            pct = 100 * row[1] / len(subset)
            lines.append(f"    {row[0]:10s}: {row[1]:4d} ({pct:5.1f}%)")

    # Confidence levels
    lines.extend(["", "=" * 60, "CONFIDENCE LEVELS", "=" * 60])
    for row in results_df.group_by("confidence").len().sort("len", descending=True).iter_rows():
        pct = 100 * row[1] / len(results)
        lines.append(f"  {row[0]:12s}: {row[1]:5d} ({pct:5.1f}%)")

    # Confidence by subtype
    lines.extend(["", "Confidence by subtype:"])
    for subtype in ["POLEmut", "MMRd", "p53abn", "NSMP"]:
        subset = results_df.filter(pl.col("subtype") == subtype)
        if len(subset) == 0:
            continue
        lines.append(f"  {subtype}:")
        for row in subset.group_by("confidence").len().sort("len", descending=True).iter_rows():
            pct = 100 * row[1] / len(subset)
            lines.append(f"    {row[0]:12s}: {row[1]:4d} ({pct:5.1f}%)")

    # Multiple classifiers
    mc = results_df.filter(pl.col("is_multiple_classifier"))
    lines.extend(["", "=" * 60, f"MULTIPLE CLASSIFIERS: {len(mc)} ({100*len(mc)/len(results):.1f}%)",
                   "=" * 60])
    if len(mc) > 0:
        for row in mc.group_by("subtype", "secondary_features").len().sort("len", descending=True).iter_rows():
            lines.append(f"  {row[0]} + {row[1]}: {row[2]}")

    # TMB by subtype
    lines.extend(["", "=" * 60, "TMB BY SUBTYPE", "=" * 60])
    for subtype in ["POLEmut", "MMRd", "p53abn", "NSMP"]:
        subset = results_df.filter(pl.col("subtype") == subtype)
        if len(subset) == 0:
            continue
        tmb_v = subset["tmb"].drop_nulls()
        if len(tmb_v) == 0:
            continue
        lines.append(
            f"  {subtype:10s}: median={tmb_v.median():.1f}, mean={tmb_v.mean():.1f}, "
            f"IQR=[{tmb_v.quantile(0.25):.1f}-{tmb_v.quantile(0.75):.1f}]"
        )

    # FGA by subtype
    if fga_lookup:
        lines.extend(["", "=" * 60, "FGA BY SUBTYPE", "=" * 60])
        for subtype in ["POLEmut", "MMRd", "p53abn", "NSMP"]:
            subset = results_df.filter(pl.col("subtype") == subtype)
            fga_v = subset["fga"].drop_nulls()
            if len(fga_v) == 0:
                continue
            lines.append(
                f"  {subtype:10s}: median={fga_v.median():.3f}, mean={fga_v.mean():.3f}, "
                f"IQR=[{fga_v.quantile(0.25):.3f}-{fga_v.quantile(0.75):.3f}]"
            )

    # MSI proxy by subtype
    lines.extend(["", "=" * 60, "MSI PROXY BY SUBTYPE", "=" * 60])
    for subtype in ["POLEmut", "MMRd", "p53abn", "NSMP"]:
        subset = results_df.filter(pl.col("subtype") == subtype)
        if len(subset) == 0:
            continue
        msi_v = subset["pseudo_msi_pct"]
        cms_v = subset["cms_genes_hit"]
        ii_v = subset["i_index"].drop_nulls()
        lines.append(
            f"  {subtype:10s}: pseudo_msi median={msi_v.median():.1f}%, "
            f"cMS genes median={cms_v.median():.0f}, I-index median={ii_v.median():.3f}"
        )

    # Indel fraction (spectrum) by subtype
    lines.extend(["", "=" * 60, "INDEL FRACTION BY SUBTYPE", "=" * 60])
    for subtype in ["POLEmut", "MMRd", "p53abn", "NSMP"]:
        subset = results_df.filter(pl.col("subtype") == subtype)
        if len(subset) == 0:
            continue
        ii_v = subset["i_index"].drop_nulls()
        lines.append(
            f"  {subtype:10s}: median={ii_v.median():.3f}, mean={ii_v.mean():.3f}, "
            f"IQR=[{ii_v.quantile(0.25):.3f}-{ii_v.quantile(0.75):.3f}]"
        )

    # Flags summary
    flagged = results_df.filter(pl.col("n_flags") > 0)
    lines.extend(["", "=" * 60, f"FLAGS: {len(flagged)} samples with flags", "=" * 60])

    summary_text = "\n".join(lines) + "\n"
    (data_dir / "validation_summary.txt").write_text(summary_text)
    print(f"\n{summary_text}")
    print(f"Output: {data_dir / 'ec_msk_classifications.tsv'}")
    print(f"Summary: {data_dir / 'validation_summary.txt'}")


if __name__ == "__main__":
    data_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/genie")
    if not (data_dir / "ec_msk_samples.tsv").exists():
        print(f"Error: run extract_genie_ec.py first — {data_dir} missing data files")
        sys.exit(1)
    main(data_dir)
