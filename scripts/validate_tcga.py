#!/usr/bin/env python3
"""Validate ec-molsubtype against TCGA UCEC PanCancer ground truth subtypes.

Runs the classifier on TCGA WES data and compares against the published
molecular subtypes (UCEC_POLE, UCEC_MSI, UCEC_CN_HIGH, UCEC_CN_LOW).

Since TCGA is WES (not panel), we use the MANTIS MSI score directly
rather than the cMS proxy.

Usage:
  python scripts/validate_tcga.py [data/tcga]
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from ec_molsubtype.classify import classify_sample
from ec_molsubtype.io import load_maf_variants_from_content, parse_maf_content
from ec_molsubtype.models import SampleInput, SampleMetadata, Variant

# TCGA subtype → our subtype mapping
TCGA_TO_SUBTYPE = {
    "UCEC_POLE": "POLEmut",
    "UCEC_MSI": "MMRd",
    "UCEC_CN_HIGH": "p53abn",
    "UCEC_CN_LOW": "NSMP",
}

# MSIsensor threshold for MSI-H (Middha et al. 2017: score >= 10)
MSISENSOR_MSI_H_THRESHOLD = 10.0
# MANTIS threshold for MSI-H (MANTIS paper: >= 0.4)
MANTIS_MSI_H_THRESHOLD = 0.4


def main(data_dir: Path) -> None:
    t0 = time.time()

    # Load clinical data
    print("Loading TCGA clinical data...")
    clin = pl.read_csv(data_dir / "tcga_ucec_clinical.tsv", separator="\t",
                       infer_schema_length=0)

    # Load mutations
    print("Loading mutations...")
    maf = pl.read_csv(data_dir / "tcga_ucec_mutations.maf", separator="\t",
                      infer_schema_length=0)
    print(f"  {len(maf)} mutations, {maf['Tumor_Sample_Barcode'].n_unique()} samples")

    # Build lookups
    clin_dict = {}
    for row in clin.iter_rows(named=True):
        clin_dict[row["SAMPLE_ID"]] = row

    # Group mutations by sample
    muts_by_sample = {}
    for sid in maf["Tumor_Sample_Barcode"].unique().to_list():
        muts_by_sample[sid] = maf.filter(pl.col("Tumor_Sample_Barcode") == sid)

    # Classify all samples with ground truth subtypes
    samples_with_gt = clin.filter(pl.col("SUBTYPE").is_not_null())
    print(f"Classifying {len(samples_with_gt)} samples with ground truth subtypes...")

    results = []
    for i, row in enumerate(samples_with_gt.iter_rows(named=True)):
        sid = row["SAMPLE_ID"]
        gt_tcga = row.get("SUBTYPE", "")
        gt_subtype = TCGA_TO_SUBTYPE.get(gt_tcga, "unknown")

        if (i + 1) % 100 == 0:
            print(f"  {i+1}/{len(samples_with_gt)}...")

        # Build variant list
        variants_df = muts_by_sample.get(sid)
        if variants_df is None:
            continue

        variant_list = []
        for vrow in variants_df.iter_rows(named=True):
            v = Variant(
                hugo_symbol=vrow.get("Hugo_Symbol", "") or "",
                chromosome=str(vrow.get("Chromosome", "") or ""),
                start_position=int(vrow.get("Start_Position", 0) or 0),
                end_position=int(vrow.get("End_Position", 0) or 0),
                reference_allele=str(vrow.get("Reference_Allele", "") or ""),
                tumor_seq_allele2=str(vrow.get("Tumor_Seq_Allele2", "") or ""),
                variant_classification=vrow.get("Variant_Classification", "") or "",
                hgvsp_short=vrow.get("HGVSp_Short", "") or "",
                t_alt_count=int(vrow.get("t_alt_count", 0) or 0),
                t_ref_count=int(vrow.get("t_ref_count", 0) or 0),
            )
            variant_list.append(v)

        # Use MSIsensor score as msi_pct (continuous, well-calibrated)
        # MSIsensor >= 10 = MSI-H (Middha et al.)
        msi_sensor_raw = row.get("MSI_SENSOR_SCORE")
        msi_sensor = float(msi_sensor_raw) if msi_sensor_raw else None

        # Convert MSIsensor score to a percentage-like scale for our classifier
        # MSIsensor scores range 0-40+; threshold at 10 for MSI-H
        # We use the score directly as msi_pct with threshold=10
        msi_pct = msi_sensor

        tmb_raw = row.get("TMB_NONSYNONYMOUS")
        tmb = float(tmb_raw) if tmb_raw else None

        fga_raw = row.get("FRACTION_GENOME_ALTERED")
        fga = float(fga_raw) if fga_raw else None

        mantis_raw = row.get("MSI_SCORE_MANTIS")
        mantis = float(mantis_raw) if mantis_raw else None

        metadata = SampleMetadata(
            sample_id=sid,
            tmb=tmb,
            msi_pct=msi_pct,
            fraction_genome_altered=fga,
        )

        sample = SampleInput(metadata=metadata, variants=variant_list)

        try:
            # Use MSIsensor threshold of 10 (standard)
            result = classify_sample(sample, msi_threshold=10.0)

            results.append({
                "SAMPLE_ID": sid,
                "gt_tcga": gt_tcga,
                "gt_subtype": gt_subtype,
                "pred_subtype": result.primary_subtype.value,
                "confidence": result.confidence.value,
                "is_multiple": result.multiple_classifier.is_multiple,
                "tmb": tmb,
                "fga": fga,
                "msi_sensor": msi_sensor,
                "mantis": mantis,
                "n_variants": len(variant_list),
                "match": gt_subtype == result.primary_subtype.value,
            })
        except Exception as e:
            print(f"  Error on {sid}: {e}")

    res = pl.DataFrame(results)
    res.write_csv(data_dir / "tcga_validation_results.tsv", separator="\t")
    print(f"  Classified {len(res)} samples")

    # === Analysis ===
    n = len(res)
    n_match = res["match"].sum()
    accuracy = 100 * n_match / n

    lines = [
        "TCGA UCEC PanCancer — ec-molsubtype Validation",
        f"Date: {time.strftime('%Y-%m-%d %H:%M')}",
        f"Runtime: {time.time() - t0:.0f}s",
        f"Samples classified: {n}",
        f"",
        f"{'='*60}",
        f"OVERALL ACCURACY: {n_match}/{n} ({accuracy:.1f}%)",
        f"{'='*60}",
    ]

    # Confusion matrix
    lines.extend(["", "CONFUSION MATRIX (rows=ground truth, cols=predicted):", ""])
    subtypes = ["POLEmut", "MMRd", "p53abn", "NSMP"]
    header = f"{'':12s}" + "".join(f"{s:>10s}" for s in subtypes) + f"{'Total':>10s}" + f"{'Recall':>10s}"
    lines.append(header)
    lines.append("-" * len(header))

    for gt in subtypes:
        gt_rows = res.filter(pl.col("gt_subtype") == gt)
        counts = []
        for pred in subtypes:
            c = gt_rows.filter(pl.col("pred_subtype") == pred).height
            counts.append(c)
        total = sum(counts)
        recall = 100 * gt_rows.filter(pl.col("match")).height / total if total > 0 else 0
        row_str = f"{gt:12s}" + "".join(f"{c:10d}" for c in counts) + f"{total:10d}" + f"{recall:9.1f}%"
        lines.append(row_str)

    # Precision row
    lines.append("-" * len(header))
    prec_str = f"{'Precision':12s}"
    for pred in subtypes:
        pred_rows = res.filter(pl.col("pred_subtype") == pred)
        tp = pred_rows.filter(pl.col("match")).height
        total_pred = len(pred_rows)
        prec = 100 * tp / total_pred if total_pred > 0 else 0
        prec_str += f"{prec:9.1f}%"
    lines.append(prec_str)

    # Per-subtype metrics
    lines.extend(["", f"{'='*60}", "PER-SUBTYPE METRICS", f"{'='*60}"])
    for st in subtypes:
        gt_n = res.filter(pl.col("gt_subtype") == st).height
        pred_n = res.filter(pl.col("pred_subtype") == st).height
        tp = res.filter((pl.col("gt_subtype") == st) & (pl.col("pred_subtype") == st)).height
        recall = 100 * tp / gt_n if gt_n > 0 else 0
        precision = 100 * tp / pred_n if pred_n > 0 else 0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0
        lines.append(f"  {st:10s}: n_gt={gt_n:3d}  n_pred={pred_n:3d}  TP={tp:3d}  "
                     f"Recall={recall:5.1f}%  Precision={precision:5.1f}%  F1={f1:5.1f}%")

    # Misclassification analysis
    mismatches = res.filter(~pl.col("match"))
    lines.extend(["", f"{'='*60}", f"MISCLASSIFICATIONS: {len(mismatches)}/{n}", f"{'='*60}"])
    for row in (mismatches.group_by("gt_subtype", "pred_subtype").len()
                .sort("len", descending=True)).iter_rows():
        lines.append(f"  {row[0]:10s} -> {row[1]:10s}: {row[2]:3d}")

    # TMB by ground truth subtype
    lines.extend(["", f"{'='*60}", "TMB BY GROUND TRUTH SUBTYPE", f"{'='*60}"])
    for st in subtypes:
        sub = res.filter(pl.col("gt_subtype") == st)
        tmb_v = sub["tmb"].drop_nulls()
        if len(tmb_v) > 0:
            lines.append(f"  {st:10s}: median={tmb_v.median():.1f}, mean={tmb_v.mean():.1f}, "
                        f"IQR=[{tmb_v.quantile(0.25):.1f}-{tmb_v.quantile(0.75):.1f}]")

    summary = "\n".join(lines) + "\n"
    (data_dir / "tcga_validation_summary.txt").write_text(summary)
    print(f"\n{summary}")


if __name__ == "__main__":
    data_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/tcga")
    main(data_dir)
