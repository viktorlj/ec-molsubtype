#!/usr/bin/env python3
"""Extract and classify EC samples from multiple GENIE V18 panels.

Extracts endometrial cancer cases from PROV, UCSF, DFCI panels
(in addition to existing MSK data) and runs the classifier on each.
Compares subtype distributions across panels.

Usage:
  python scripts/validate_multi_panel.py /path/to/Genie_V18 [data/genie]
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from ec_molsubtype.classify import classify_sample
from ec_molsubtype.models import SampleInput, SampleMetadata, Variant
from ec_molsubtype.msi import compute_msi_proxy

# All endometrial cancer OncoTree codes
EC_ONCOTREE_CODES = {
    "UCEC", "UEC", "USC", "UCS", "UCCC", "UMEC", "UUC", "UDDC",
    "UPDC", "UNEC", "UASC", "UMC", "USCC", "USMT",
}

# Panels to validate (must have POLE + TP53 at minimum)
PANELS = {
    "MSK-IMPACT468":      {"size_mb": 1.171, "msi_threshold": 9.0},
    "MSK-IMPACT505":      {"size_mb": 1.261, "msi_threshold": 9.0},
    "PROV-TSO500HT-V2":   {"size_mb": 1.94,  "msi_threshold": 9.0},
    "UCSF-IDTV5-TO":      {"size_mb": 1.50,  "msi_threshold": 9.0},
    "DFCI-ONCOPANEL-3.1":  {"size_mb": 1.35,  "msi_threshold": 9.0},
    "DFCI-ONCOPANEL-3":    {"size_mb": 1.35,  "msi_threshold": 9.0},
}


def classify_samples_for_panel(
    panel_id: str,
    panel_config: dict,
    sample_ids: set[str],
    mutations_by_sample: dict[str, pl.DataFrame],
    tmb_lookup: dict[str, float],
    fga_lookup: dict[str, float],
    sample_info: dict[str, dict],
) -> list[dict]:
    """Classify all samples for a given panel."""
    results = []
    msi_threshold = panel_config["msi_threshold"]

    for sid in sorted(sample_ids):
        variants_df = mutations_by_sample.get(sid)
        if variants_df is None or len(variants_df) == 0:
            continue

        variant_list = []
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

        msi_proxy = compute_msi_proxy(variant_list, panel_id=panel_id)
        tmb_val = tmb_lookup.get(sid)
        fga_val = fga_lookup.get(sid)

        metadata = SampleMetadata(
            sample_id=sid,
            tmb=tmb_val,
            panel_size_mb=panel_config.get("size_mb"),
            msi_pct=msi_proxy.pseudo_msi_pct,
            fraction_genome_altered=fga_val,
        )
        sample_input = SampleInput(metadata=metadata, variants=variant_list)

        try:
            result = classify_sample(sample_input, msi_threshold=msi_threshold)
            info = sample_info.get(sid, {})
            results.append({
                "SAMPLE_ID": sid,
                "SEQ_ASSAY_ID": panel_id,
                "ONCOTREE_CODE": info.get("ONCOTREE_CODE", ""),
                "subtype": result.primary_subtype.value,
                "confidence": result.confidence.value,
                "is_multiple": result.multiple_classifier.is_multiple,
                "tmb": tmb_val,
                "fga": fga_val,
                "pseudo_msi_pct": msi_proxy.pseudo_msi_pct,
                "cms_high_spec": len(msi_proxy.cms_high_spec_hit),
                "i_index": msi_proxy.i_index,
                "n_mutations": msi_proxy.n_mutations,
            })
        except Exception as e:
            print(f"  Error classifying {sid}: {e}")

    return results


def main(genie_dir: Path, output_dir: Path) -> None:
    t0 = time.time()

    # Load clinical data
    print("Loading clinical data...")
    clinical = pl.read_csv(
        genie_dir / "data_clinical_sample.txt",
        separator="\t", comment_prefix="#", infer_schema_length=0,
    )

    # Filter for EC samples on target panels
    ec_on_panels = clinical.filter(
        pl.col("ONCOTREE_CODE").is_in(EC_ONCOTREE_CODES)
        & pl.col("SEQ_ASSAY_ID").is_in(set(PANELS.keys()))
    )
    print(f"  Total EC samples across {len(PANELS)} panels: {len(ec_on_panels)}")
    for row in ec_on_panels.group_by("SEQ_ASSAY_ID").len().sort("len", descending=True).iter_rows():
        print(f"    {row[0]}: {row[1]}")

    all_sample_ids = set(ec_on_panels["SAMPLE_ID"].to_list())

    # Build sample info lookup
    sample_info = {}
    for row in ec_on_panels.iter_rows(named=True):
        sample_info[row["SAMPLE_ID"]] = row

    # Load TMB
    print("Loading TMB...")
    tmb_df = pl.read_csv(genie_dir / "tmb_18.0-public.tsv", separator="\t")
    # GENIE anonymizes TMB — compute from mutations later if needed

    # Load mutations
    print("Loading mutations (this takes a moment)...")
    mutations = pl.scan_csv(
        genie_dir / "data_mutations_extended.txt",
        separator="\t", comment_prefix="#", infer_schema_length=10000,
        null_values=["", "NA", "nan"],
        schema_overrides={"Chromosome": pl.Utf8, "FILTER": pl.Utf8,
                          "SWISSPROT": pl.Utf8, "Codons": pl.Utf8},
    ).filter(
        pl.col("Tumor_Sample_Barcode").is_in(all_sample_ids)
    ).collect()
    print(f"  {len(mutations)} mutations for {mutations['Tumor_Sample_Barcode'].n_unique()} samples")

    # Group mutations by sample
    mutations_by_sample: dict[str, pl.DataFrame] = {}
    for sid in mutations["Tumor_Sample_Barcode"].unique().to_list():
        mutations_by_sample[sid] = mutations.filter(pl.col("Tumor_Sample_Barcode") == sid)

    # Compute TMB from mutations + panel size
    non_silent = {"Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank", "IGR", "RNA"}
    tmb_lookup: dict[str, float] = {}
    for sid, df in mutations_by_sample.items():
        ns_count = df.filter(~pl.col("Variant_Classification").is_in(non_silent)).height
        panel_id = sample_info.get(sid, {}).get("SEQ_ASSAY_ID", "")
        panel_size = PANELS.get(panel_id, {}).get("size_mb", 1.0)
        tmb_lookup[sid] = ns_count / panel_size

    # Load CNA/FGA if available (from MSK extraction)
    fga_lookup: dict[str, float] = {}
    cin_file = output_dir / "ec_msk_cin_metrics.tsv"
    if cin_file.exists():
        cin_df = pl.read_csv(cin_file, separator="\t")
        fga_lookup = dict(zip(cin_df["SAMPLE_ID"].to_list(), cin_df["fga"].to_list()))

    # Load CNA for non-MSK panels
    print("Loading CNA segments for non-MSK panels...")
    non_msk_ids = {sid for sid in all_sample_ids if not sid.startswith("GENIE-MSK")}
    if non_msk_ids:
        seg = pl.scan_csv(
            genie_dir / "data_cna_hg19.seg",
            separator="\t",
            schema_overrides={"chrom": pl.Utf8},
        ).filter(
            pl.col("ID").is_in(non_msk_ids)
        ).collect()

        # Compute FGA for non-MSK samples
        if len(seg) > 0:
            seg_with_len = seg.with_columns([
                (pl.col("loc.end") - pl.col("loc.start")).alias("seg_length"),
                (pl.col("seg.mean").abs() > 0.2).alias("is_altered"),
            ])
            fga_per_sample = (
                seg_with_len.group_by("ID").agg([
                    pl.col("seg_length").sum().alias("total"),
                    (pl.col("seg_length") * pl.col("is_altered").cast(pl.Int64)).sum().alias("altered"),
                ])
                .with_columns((pl.col("altered") / pl.col("total")).alias("fga"))
            )
            for row in fga_per_sample.iter_rows(named=True):
                fga_lookup[row["ID"]] = row["fga"]
            print(f"  FGA computed for {len(fga_per_sample)} non-MSK samples")

    # Classify per panel
    all_results = []
    for panel_id, config in PANELS.items():
        panel_samples = set(
            ec_on_panels.filter(pl.col("SEQ_ASSAY_ID") == panel_id)["SAMPLE_ID"].to_list()
        )
        if not panel_samples:
            continue
        print(f"\nClassifying {panel_id} ({len(panel_samples)} samples)...")
        panel_results = classify_samples_for_panel(
            panel_id, config, panel_samples,
            mutations_by_sample, tmb_lookup, fga_lookup, sample_info,
        )
        all_results.extend(panel_results)
        print(f"  Classified {len(panel_results)} samples")

    # Save all results
    results_df = pl.DataFrame(all_results)
    results_df.write_csv(output_dir / "ec_multi_panel_classifications.tsv", separator="\t")

    # === Summary ===
    lines = [
        "GENIE V18 EC — Multi-Panel Classification Validation",
        f"Date: {time.strftime('%Y-%m-%d %H:%M')}",
        f"Runtime: {time.time() - t0:.0f}s",
        f"Total samples classified: {len(results_df)}",
        "",
    ]

    # Per-panel subtype distribution
    lines.extend(["=" * 70, "SUBTYPE DISTRIBUTION BY PANEL", "=" * 70])
    for panel_id in PANELS:
        panel_df = results_df.filter(pl.col("SEQ_ASSAY_ID") == panel_id)
        if len(panel_df) == 0:
            continue
        lines.append(f"\n{panel_id} (n={len(panel_df)}):")
        for row in panel_df.group_by("subtype").len().sort("len", descending=True).iter_rows():
            pct = 100 * row[1] / len(panel_df)
            lines.append(f"  {row[0]:10s}: {row[1]:4d} ({pct:5.1f}%)")

        # TMB by subtype
        for subtype in ["POLEmut", "MMRd", "p53abn", "NSMP"]:
            subset = panel_df.filter(pl.col("subtype") == subtype)
            if len(subset) == 0:
                continue
            tmb_v = subset["tmb"].drop_nulls()
            if len(tmb_v) > 0:
                lines.append(
                    f"    {subtype} TMB: median={tmb_v.median():.1f}, "
                    f"I-index={subset['i_index'].drop_nulls().median():.3f}"
                )

    # Cross-panel comparison table
    lines.extend(["", "=" * 70, "CROSS-PANEL COMPARISON (% by subtype)", "=" * 70])
    header = f"{'Panel':25s} {'n':>5s} {'POLEmut':>8s} {'MMRd':>8s} {'p53abn':>8s} {'NSMP':>8s}"
    lines.append(header)
    lines.append("-" * len(header))
    for panel_id in PANELS:
        panel_df = results_df.filter(pl.col("SEQ_ASSAY_ID") == panel_id)
        n = len(panel_df)
        if n == 0:
            continue
        pcts = {}
        for subtype in ["POLEmut", "MMRd", "p53abn", "NSMP"]:
            cnt = panel_df.filter(pl.col("subtype") == subtype).height
            pcts[subtype] = f"{100*cnt/n:.1f}%"
        lines.append(
            f"{panel_id:25s} {n:5d} {pcts['POLEmut']:>8s} {pcts['MMRd']:>8s} "
            f"{pcts['p53abn']:>8s} {pcts['NSMP']:>8s}"
        )
    lines.append(f"{'TCGA reference':25s} {'':>5s} {'~7%':>8s} {'~28%':>8s} {'~26%':>8s} {'~39%':>8s}")

    # UEC-specific (most comparable to TCGA)
    lines.extend(["", "=" * 70, "UEC (ENDOMETRIOID) ONLY — Most comparable to TCGA", "=" * 70])
    uec = results_df.filter(pl.col("ONCOTREE_CODE") == "UEC")
    header2 = f"{'Panel':25s} {'n':>5s} {'POLEmut':>8s} {'MMRd':>8s} {'p53abn':>8s} {'NSMP':>8s}"
    lines.append(header2)
    lines.append("-" * len(header2))
    for panel_id in PANELS:
        panel_uec = uec.filter(pl.col("SEQ_ASSAY_ID") == panel_id)
        n = len(panel_uec)
        if n < 20:
            continue
        pcts = {}
        for subtype in ["POLEmut", "MMRd", "p53abn", "NSMP"]:
            cnt = panel_uec.filter(pl.col("subtype") == subtype).height
            pcts[subtype] = f"{100*cnt/n:.1f}%"
        lines.append(
            f"{panel_id:25s} {n:5d} {pcts['POLEmut']:>8s} {pcts['MMRd']:>8s} "
            f"{pcts['p53abn']:>8s} {pcts['NSMP']:>8s}"
        )

    # Secondary evidence concordance
    lines.extend(["", "=" * 70, "SECONDARY EVIDENCE (ALL PANELS POOLED)", "=" * 70])
    for subtype in ["POLEmut", "MMRd", "p53abn", "NSMP"]:
        subset = results_df.filter(pl.col("subtype") == subtype)
        if len(subset) == 0:
            continue
        tmb_v = subset["tmb"].drop_nulls()
        fga_v = subset["fga"].drop_nulls()
        ii_v = subset["i_index"].drop_nulls()
        lines.append(
            f"  {subtype:10s} (n={len(subset):4d}): "
            f"TMB={tmb_v.median():.1f}, "
            f"FGA={fga_v.median():.3f}, "
            f"I-index={ii_v.median():.3f}"
        )

    summary = "\n".join(lines) + "\n"
    (output_dir / "multi_panel_validation_summary.txt").write_text(summary)
    print(f"\n{summary}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <genie_v18_dir> [output_dir]")
        sys.exit(1)
    genie_dir = Path(sys.argv[1])
    output_dir = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("data/genie")
    output_dir.mkdir(parents=True, exist_ok=True)
    main(genie_dir, output_dir)
