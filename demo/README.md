# Demo Data

Synthetic endometrial cancer samples for testing the ec-molsubtype classifier. Each sample demonstrates a specific subtype or edge case in the WHO 2020 molecular classification.

**These are fabricated data for demonstration purposes only. They do not represent real patients.**

## Samples

### Clear-cut subtypes

| File | Expected subtype | Description |
|------|-----------------|-------------|
| `polemut_clear` | **POLEmut** | POLE P286R (tier 1 hotspot) + concurrent TP53 R248W. TMB 232, low CNA, SBS10a/10b dominant. Multiple classifier: reports p53abn as secondary feature but classifies as POLEmut per WHO hierarchy. |
| `mmrd_msih` | **MMRd** | MSI-H (48% unstable markers) with MSH6 frameshift + frameshifts in RNF43, JAK1, AXIN2, B2M, KMT2C. TMB 42, high indel fraction. Classic MSI-H endometrioid pattern. |
| `p53abn_serous` | **p53abn** | TP53 R175H (DNA-binding domain hotspot) in serous-like context. TMB 3.8, high CNA (FGA 0.58). Co-mutations: PIK3CA, PPP2R1A, FBXW7 — typical for serous histology. |
| `nsmp_endometrioid` | **NSMP** | No POLE EDM, MSS, no pathogenic TP53. Classic low-grade endometrioid with PTEN, PIK3CA, ARID1A, CTNNB1, KRAS. TMB 2.4, low CNA. |

### Edge cases

| File | Expected subtype | Description |
|------|-----------------|-------------|
| `polemut_vus` | **NSMP** (flagged) | POLE A300T — a VUS in the exonuclease domain (codon 268-471) but not on the curated pathogenic list. Despite TMB 145 and SBS10a/10b signatures, the VUS is NOT auto-classified as POLEmut. Flagged for manual review with secondary evidence. Demonstrates conservative VUS handling. |
| `mmrd_biallelic` | **MMRd** | Biallelic MLH1 inactivation (nonsense + frameshift) with low MSI marker percentage (5%). Demonstrates the biallelic MMR pathway: classified as MMRd despite discordant MSI markers, with a discordance flag. |
| `multiple_pole_p53` | **POLEmut** (triple) | POLE V411L + MSI-H (32%) + TP53 R273H + biallelic MSH6. A rare triple-classifier case. Classified as POLEmut per hierarchy, with MMRd and p53abn reported as secondary features. Clinical note highlights that immunotherapy may still benefit despite POLEmut primary classification. |

## File structure

Each sample has two files:
- **`<name>.maf`** — Somatic variants in MAF format (tab-separated)
- **`<name>.json`** — Sample metadata sidecar (TMB, MSI%, CNA, signature weights)

Plus:
- **`metadata.tsv`** — Batch metadata for all samples (for use with `classify-batch`)
- **`results.tsv`** — Pre-computed batch classification results

## Usage

### Single sample

```bash
ec-molsubtype classify demo/polemut_clear.maf \
  -m demo/polemut_clear.json \
  --human
```

### Batch mode

```bash
ec-molsubtype classify-batch \
  --input-dir demo/ \
  --metadata-tsv demo/metadata.tsv \
  --output demo/results.tsv
```

### JSON output

```bash
ec-molsubtype classify demo/multiple_pole_p53.maf \
  -m demo/multiple_pole_p53.json \
  -o result.json
```

## What to look for

- **POLEmut cases**: Check that POLE tier 1/2 hotspots are identified, TMB >100 is concordant, and concurrent TP53/MMRd features are reported as secondary.
- **MMRd cases**: Both MSI-H driven (`mmrd_msih`) and biallelic MMR gene driven (`mmrd_biallelic`) pathways should classify correctly. Note the discordance flag on the biallelic case.
- **p53abn**: TP53 hotspot recognition, high CNA concordance, serous-associated co-mutations.
- **NSMP**: Diagnosis of exclusion — no POLE EDM, no MSI-H, no pathogenic TP53.
- **POLE VUS**: The `polemut_vus` sample demonstrates that VUS in the exonuclease domain are NOT auto-classified, even with suggestive secondary evidence. This is intentional — VUS require pathologist review.
- **Multiple classifiers**: The `polemut_clear` and `multiple_pole_p53` samples show how concurrent subtype features are reported alongside the hierarchical primary classification.

## MAF column reference

| Column | Required | Description |
|--------|----------|-------------|
| Hugo_Symbol | Yes | Gene symbol (e.g., POLE, TP53) |
| Variant_Classification | Yes | Effect type (Missense_Mutation, Frame_Shift_Del, etc.) |
| HGVSp_Short | Yes | Protein change (e.g., p.P286R) |
| Reference_Allele | No | Reference allele (used for spectrum analysis) |
| Tumor_Seq_Allele2 | No | Alternate allele (used for spectrum analysis) |
| t_alt_count | No | Alt allele read count (used for VAF) |
| t_ref_count | No | Ref allele read count (used for VAF) |
| Chromosome | No | Chromosome |
| Start_Position | No | Genomic position |
| Tumor_Sample_Barcode | No | Sample ID (fallback if not in metadata) |

## Metadata sidecar fields

| Field | Type | Description |
|-------|------|-------------|
| sample_id | string | Sample identifier |
| tmb | float | Tumor mutational burden (mutations/Mb) |
| panel_size_mb | float | Sequencing panel size in Mb |
| msi_pct | float | % unstable MSI markers (0-100). Primary MMRd input. |
| fraction_genome_altered | float | CNA burden (0-1). Secondary evidence for p53abn. |
| signature_weights | object | SBS signature weights from decomposition tools |
