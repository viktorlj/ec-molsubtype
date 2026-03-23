#!/usr/bin/env python3
"""Download CPTAC endometrial cancer data from cBioPortal API.

Downloads clinical data (subtypes, MSI status, TMB, FGA) and mutation data
from the CPTAC UCEC 2020 study for independent validation.

Produces:
  data/cptac/cptac_ucec_clinical.tsv   — clinical + molecular annotations
  data/cptac/cptac_ucec_mutations.maf  — mutation data in MAF-like format
  data/cptac/download_summary.txt      — summary statistics
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import httpx
import polars as pl

BASE_URL = "https://www.cbioportal.org/api"
STUDY_ID = "ucec_cptac_2020"


def fetch_json(path: str, params: dict | None = None) -> list | dict:
    url = f"{BASE_URL}{path}"
    r = httpx.get(url, params=params, timeout=60)
    r.raise_for_status()
    return r.json()


def fetch_post_json(path: str, body: dict) -> list | dict:
    url = f"{BASE_URL}{path}"
    r = httpx.post(url, json=body, timeout=120,
                   headers={"Content-Type": "application/json"})
    r.raise_for_status()
    return r.json()


def main(output_dir: Path) -> None:
    t0 = time.time()
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- 1. Get sample list ---
    print("Fetching sample list...")
    samples = fetch_json(f"/studies/{STUDY_ID}/samples")
    sample_ids = [s["sampleId"] for s in samples]
    patient_ids = list({s["patientId"] for s in samples})
    print(f"  {len(sample_ids)} samples, {len(patient_ids)} patients")

    # --- 2. Fetch clinical data ---
    print("Fetching clinical data...")
    # Fetch all sample-level clinical data
    clinical_data = fetch_json(
        f"/studies/{STUDY_ID}/clinical-data",
        {"clinicalDataType": "SAMPLE", "projection": "SUMMARY"}
    )
    rows: dict[str, dict] = {sid: {"SAMPLE_ID": sid} for sid in sample_ids}
    for item in clinical_data:
        sid = item["sampleId"]
        attr = item["clinicalAttributeId"]
        val = item["value"]
        if sid in rows:
            rows[sid][attr] = val

    # Also fetch patient-level clinical
    patient_clinical = fetch_json(
        f"/studies/{STUDY_ID}/clinical-data",
        {"clinicalDataType": "PATIENT", "projection": "SUMMARY"}
    )
    patient_map: dict[str, dict] = {}
    for item in patient_clinical:
        pid = item["patientId"]
        patient_map.setdefault(pid, {})
        patient_map[pid][item["clinicalAttributeId"]] = item["value"]

    # Merge patient-level into sample rows
    for s in samples:
        sid, pid = s["sampleId"], s["patientId"]
        rows[sid]["PATIENT_ID"] = pid
        for k, v in patient_map.get(pid, {}).items():
            if k not in rows[sid]:
                rows[sid][k] = v

    clinical_df = pl.DataFrame(list(rows.values()))
    clinical_df.write_csv(output_dir / "cptac_ucec_clinical.tsv", separator="\t")
    print(f"  Wrote {output_dir / 'cptac_ucec_clinical.tsv'}")

    # --- 3. Summarize key clinical attributes ---
    print("\n  Available clinical columns:")
    for col in sorted(clinical_df.columns):
        n_non_null = clinical_df.filter(pl.col(col).is_not_null() & (pl.col(col) != "")).height
        if n_non_null > 0:
            print(f"    {col}: {n_non_null}/{clinical_df.height} values")

    # Check for subtype-related columns
    subtype_cols = [c for c in clinical_df.columns if "SUBTYPE" in c.upper()]
    msi_cols = [c for c in clinical_df.columns if "MSI" in c.upper() or "MICRO" in c.upper()]
    tmb_cols = [c for c in clinical_df.columns if "TMB" in c.upper() or "MUTATION" in c.upper()]
    fga_cols = [c for c in clinical_df.columns if "FGA" in c.upper() or "FRACTION_GENOME" in c.upper()]

    print(f"\n  Subtype columns found: {subtype_cols}")
    print(f"  MSI columns found: {msi_cols}")
    print(f"  TMB/mutation columns found: {tmb_cols}")
    print(f"  FGA columns found: {fga_cols}")

    # Print subtype distribution if available
    for col in subtype_cols:
        if col in clinical_df.columns:
            print(f"\n  {col} distribution:")
            for row in clinical_df.group_by(col).len().sort("len", descending=True).iter_rows():
                print(f"    {row[0]}: {row[1]}")

    # Print MSI distribution if available
    for col in msi_cols:
        if col in clinical_df.columns:
            print(f"\n  {col} distribution:")
            for row in clinical_df.group_by(col).len().sort("len", descending=True).iter_rows():
                print(f"    {row[0]}: {row[1]}")

    # --- 4. Fetch mutations ---
    print("\nFetching mutations (this may take a moment)...")
    # Get the molecular profile ID for mutations
    profiles = fetch_json(f"/studies/{STUDY_ID}/molecular-profiles")
    mut_profile = None
    for p in profiles:
        if p["molecularAlterationType"] == "MUTATION_EXTENDED":
            mut_profile = p["molecularProfileId"]
            break

    if not mut_profile:
        print("  ERROR: No mutation profile found!")
        print("  Available profiles:")
        for p in profiles:
            print(f"    {p['molecularProfileId']} ({p['molecularAlterationType']})")
        return

    print(f"  Mutation profile: {mut_profile}")

    # Fetch mutations in batches
    all_mutations = []
    batch_size = 100
    for i in range(0, len(sample_ids), batch_size):
        batch = sample_ids[i:i + batch_size]
        body = {
            "sampleMolecularIdentifiers": [
                {"molecularProfileId": mut_profile, "sampleId": sid}
                for sid in batch
            ],
        }
        muts = fetch_post_json(
            "/mutations/fetch",
            body,
        )
        all_mutations.extend(muts)
        if (i + batch_size) % 500 == 0 or i + batch_size >= len(sample_ids):
            print(f"  Fetched {len(all_mutations)} mutations ({i + len(batch)}/{len(sample_ids)} samples)...")

    print(f"  Total mutations: {len(all_mutations)}")

    # Convert to MAF-like format
    # Gene name is in the keyword field (e.g., "TP53 R248 missense")
    maf_rows = []
    for m in all_mutations:
        keyword = m.get("keyword", "")
        hugo = keyword.split()[0] if keyword else ""
        maf_rows.append({
            "Hugo_Symbol": hugo,
            "Entrez_Gene_Id": m.get("entrezGeneId", ""),
            "Chromosome": m.get("chr", ""),
            "Start_Position": m.get("startPosition", ""),
            "End_Position": m.get("endPosition", ""),
            "Reference_Allele": m.get("referenceAllele", ""),
            "Tumor_Seq_Allele2": m.get("variantAllele", ""),
            "Variant_Classification": m.get("mutationType", ""),
            "HGVSp_Short": f"p.{m['proteinChange']}" if m.get("proteinChange") else "",
            "Variant_Type": m.get("variantType", ""),
            "t_alt_count": m.get("tumorAltCount", ""),
            "t_ref_count": m.get("tumorRefCount", ""),
            "Tumor_Sample_Barcode": m.get("sampleId", ""),
            "Mutation_Status": m.get("mutationStatus", ""),
            "Keyword": m.get("keyword", ""),
        })

    maf_df = pl.DataFrame(maf_rows)
    maf_df.write_csv(output_dir / "cptac_ucec_mutations.maf", separator="\t")
    print(f"  Wrote {output_dir / 'cptac_ucec_mutations.maf'}")

    # --- 5. Verify key gene counts ---
    print("\n  Key gene mutation counts:")
    key_genes = ["POLE", "TP53", "MLH1", "MSH2", "MSH6", "PMS2",
                 "PTEN", "PIK3CA", "ARID1A", "KRAS", "CTNNB1",
                 "PPP2R1A", "FBXW7"]
    for gene in key_genes:
        count = maf_df.filter(pl.col("Hugo_Symbol") == gene).height
        n_samples = maf_df.filter(pl.col("Hugo_Symbol") == gene)["Tumor_Sample_Barcode"].n_unique()
        if count > 0:
            print(f"    {gene}: {count} mutations in {n_samples} samples")
        else:
            print(f"    {gene}: 0 mutations")

    # Check Hugo_Symbol population
    n_empty_hugo = maf_df.filter(
        (pl.col("Hugo_Symbol").is_null()) | (pl.col("Hugo_Symbol") == "")
    ).height
    print(f"\n  Rows with empty Hugo_Symbol: {n_empty_hugo}/{maf_df.height}")

    # --- 6. Summary ---
    elapsed = time.time() - t0

    # Find the best subtype column
    subtype_col = subtype_cols[0] if subtype_cols else None
    n_with_subtype = (
        clinical_df.filter(pl.col(subtype_col).is_not_null() & (pl.col(subtype_col) != "")).height
        if subtype_col else 0
    )
    msi_col = msi_cols[0] if msi_cols else None
    n_with_msi = (
        clinical_df.filter(pl.col(msi_col).is_not_null() & (pl.col(msi_col) != "")).height
        if msi_col else 0
    )

    lines = [
        f"CPTAC UCEC 2020 — Download Summary",
        f"Date: {time.strftime('%Y-%m-%d %H:%M')}",
        f"Runtime: {elapsed:.0f}s",
        f"",
        f"Samples: {len(sample_ids)}",
        f"Patients: {len(patient_ids)}",
        f"Mutations: {len(all_mutations)}",
        f"Samples with subtype ({subtype_col}): {n_with_subtype}",
        f"Samples with MSI data ({msi_col}): {n_with_msi}",
        f"Samples with mutations: {maf_df['Tumor_Sample_Barcode'].n_unique()}",
        f"Rows with empty Hugo_Symbol: {n_empty_hugo}",
        f"",
        f"Key genes: {', '.join(key_genes)}",
    ]
    summary = "\n".join(lines) + "\n"
    (output_dir / "download_summary.txt").write_text(summary)
    print(f"\n{summary}")


if __name__ == "__main__":
    output_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/cptac")
    main(output_dir)
