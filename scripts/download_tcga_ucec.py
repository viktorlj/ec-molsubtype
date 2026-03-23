#!/usr/bin/env python3
"""Download TCGA UCEC PanCancer Atlas data from cBioPortal API.

Downloads clinical data (subtypes, MSI scores, TMB, FGA) and mutation data
for independent validation of the ec-molsubtype classifier.

Produces:
  data/tcga/tcga_ucec_clinical.tsv   — clinical + molecular annotations
  data/tcga/tcga_ucec_mutations.maf  — mutation data
  data/tcga/download_summary.txt     — summary statistics
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import httpx
import polars as pl

BASE_URL = "https://www.cbioportal.org/api"
STUDY_ID = "ucec_tcga_pan_can_atlas_2018"

# Clinical attributes to fetch
CLINICAL_ATTRS = [
    "SUBTYPE", "MSI_SCORE_MANTIS", "MSI_SENSOR_SCORE",
    "TMB_NONSYNONYMOUS", "FRACTION_GENOME_ALTERED", "MUTATION_COUNT",
    "ANEUPLOIDY_SCORE", "ONCOTREE_CODE", "CANCER_TYPE_DETAILED",
    "AJCC_PATHOLOGIC_TUMOR_STAGE", "GRADE",
]


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
    clinical_df.write_csv(output_dir / "tcga_ucec_clinical.tsv", separator="\t")
    print(f"  Wrote {output_dir / 'tcga_ucec_clinical.tsv'}")

    # --- 3. Subtype summary ---
    if "SUBTYPE" in clinical_df.columns:
        print("  Subtype distribution:")
        for row in clinical_df.group_by("SUBTYPE").len().sort("len", descending=True).iter_rows():
            print(f"    {row[0]}: {row[1]}")

    # --- 4. Fetch mutations ---
    print("Fetching mutations (this may take a moment)...")
    # Get the molecular profile ID for mutations
    profiles = fetch_json(f"/studies/{STUDY_ID}/molecular-profiles")
    mut_profile = None
    for p in profiles:
        if p["molecularAlterationType"] == "MUTATION_EXTENDED":
            mut_profile = p["molecularProfileId"]
            break

    if not mut_profile:
        print("  ERROR: No mutation profile found!")
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
            f"/mutations/fetch",
            body,
        )
        all_mutations.extend(muts)
        if (i + batch_size) % 500 == 0 or i + batch_size >= len(sample_ids):
            print(f"  Fetched {len(all_mutations)} mutations ({i + len(batch)}/{len(sample_ids)} samples)...")

    print(f"  Total mutations: {len(all_mutations)}")

    # Convert to MAF-like format
    # Gene name is in the keyword field (e.g., "CHD1 G1129 missense")
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
    maf_df.write_csv(output_dir / "tcga_ucec_mutations.maf", separator="\t")
    print(f"  Wrote {output_dir / 'tcga_ucec_mutations.maf'}")

    # --- 5. Summary ---
    elapsed = time.time() - t0
    n_with_subtype = clinical_df.filter(pl.col("SUBTYPE").is_not_null()).height if "SUBTYPE" in clinical_df.columns else 0
    n_with_msi = clinical_df.filter(pl.col("MSI_SCORE_MANTIS").is_not_null()).height if "MSI_SCORE_MANTIS" in clinical_df.columns else 0

    lines = [
        f"TCGA UCEC PanCancer Atlas — Download Summary",
        f"Date: {time.strftime('%Y-%m-%d %H:%M')}",
        f"Runtime: {elapsed:.0f}s",
        f"",
        f"Samples: {len(sample_ids)}",
        f"Patients: {len(patient_ids)}",
        f"Mutations: {len(all_mutations)}",
        f"Samples with subtype: {n_with_subtype}",
        f"Samples with MANTIS score: {n_with_msi}",
        f"Samples with mutations: {maf_df['Tumor_Sample_Barcode'].n_unique()}",
    ]
    summary = "\n".join(lines) + "\n"
    (output_dir / "download_summary.txt").write_text(summary)
    print(f"\n{summary}")


if __name__ == "__main__":
    output_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/tcga")
    main(output_dir)
