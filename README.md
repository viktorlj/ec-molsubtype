# ec-molsubtype

Molecular subtyping of endometrial cancer from panel sequencing data, following the **WHO 2020** classification framework, **FIGO 2023** molecular staging, and **Swedish national guidelines** (Nationellt vardprogram livmoderkroppscancer v3.2, 2026-03-10).

Classifies tumors into four prognostic molecular subtypes using a hierarchical sequential algorithm with secondary genomic evidence for confidence scoring.

## The four molecular subtypes

| Priority | Subtype     | Marker                              | Prognosis     | Clinical implication (Swedish guidelines)                |
|----------|-------------|--------------------------------------|---------------|----------------------------------------------------------|
| 1        | **POLEmut** | Pathogenic POLE exonuclease domain mutation | Best          | No adjuvant treatment for stage I--II                   |
| 2        | **MMRd**    | Mismatch repair deficiency           | Intermediate  | Checkpoint inhibitor benefit; Lynch syndrome screening   |
| 3        | **p53abn**  | Pathogenic TP53 mutation             | Worst         | Postoperative chemotherapy for stage I--II               |
| 4        | **NSMP**    | No specific molecular profile        | Intermediate  | Diagnosis of exclusion; ER status sub-stratifies         |

## Classification algorithm

The WHO 2020 hierarchical sequential classifier:

```
Step 1: POLE pathogenic EDM (exons 9-14)?  --> YES --> POLEmut
Step 2: MMR deficiency (MSI-H or biallelic MMR)?  --> YES --> MMRd
Step 3: TP53 pathogenic mutation?  --> YES --> p53abn
Step 4: None of the above  --> NSMP
```

Priority is absolute: a tumor with both a POLE hotspot and a TP53 mutation is classified as **POLEmut** (multiple-classifier cases are reported but do not change the primary assignment).

### Secondary evidence layer

Each classification is scored against concordant secondary genomic features:

| Subtype | Expected TMB | Expected CNA | Expected signatures |
|---------|-------------|--------------|---------------------|
| POLEmut | >100 mut/Mb | Low          | SBS10a, SBS10b      |
| MMRd    | 10--100 mut/Mb | Low--intermediate | SBS6, SBS15, SBS21, SBS26 |
| p53abn  | <10 mut/Mb  | High         | SBS1, SBS5          |
| NSMP    | <10 mut/Mb  | Low          | SBS1, SBS5          |

Confidence levels: **High** (all concordant), **Moderate** (some discordance), **Low** (borderline/VUS), **Discordant** (primary contradicted by secondary evidence).

## Installation

Requires Python >= 3.11 and [uv](https://docs.astral.sh/uv/).

```bash
# Clone and install
git clone https://github.com/<your-org>/ec-molsubtype.git
cd ec-molsubtype
uv venv
source .venv/bin/activate
uv pip install -e ".[dev]"
```

## Usage

### Classify a single sample

```bash
ec-molsubtype classify sample.maf --metadata sample_meta.json --output result.json
```

### Classify a batch

```bash
ec-molsubtype classify-batch --input-dir ./samples/ --output results.tsv
```

### Inspect a POLE variant

```bash
ec-molsubtype check-pole "p.P286R" --tmb 250.0
```

### Input format

**MAF file** (minimum required columns):

| Column | Description |
|--------|-------------|
| `Hugo_Symbol` | Gene symbol |
| `Variant_Classification` | Missense_Mutation, Nonsense_Mutation, Frame_Shift_Del, etc. |
| `HGVSp_Short` | Protein change (e.g., p.P286R) |
| `t_alt_count` | Tumor alt allele count (optional, for VAF) |
| `t_ref_count` | Tumor ref allele count (optional, for VAF) |

**Sample metadata** (JSON sidecar):

```json
{
  "sample_id": "SAMPLE_001",
  "tmb": 285.3,
  "panel_size_mb": 1.2,
  "msi_pct": 0.0,
  "fraction_genome_altered": 0.05,
  "signature_weights": {"SBS10a": 0.45, "SBS10b": 0.22}
}
```

`msi_pct` is the percentage of unstable markers from the panel's MSI marker set (0.0--100.0). This is the primary input for MMRd assessment. If unavailable, biallelic MMR gene mutations in the MAF are used as an alternative.

### Output

JSON per sample with classification result, evidence, and clinical notes:

```json
{
  "sample_id": "SAMPLE_001",
  "primary_subtype": "POLEmut",
  "confidence": "high",
  "classification_path": [
    {"step": 1, "test": "POLE_EDM", "result": "positive", "variant": "POLE p.P286R", "tier": 1}
  ],
  "secondary_evidence": {
    "tmb": {"value": 285.3, "expected_range": ">100", "concordant": true}
  },
  "clinical_notes": [
    "POLE P286R is a tier 1 established pathogenic hotspot.",
    "Swedish guidelines: no adjuvant treatment recommended for stage I-II POLEmut."
  ]
}
```

Batch output is TSV with one row per sample.

## POLE variant tiers

| Tier | Variants | Action |
|------|----------|--------|
| **Tier 1** (established hotspots) | P286R, V411L, S297F, A456P, S459F | Auto-classify as POLEmut |
| **Tier 2** (confirmed pathogenic) | F367S, L424I, L424V, M444K, D368N, P286H | Auto-classify as POLEmut |
| **VUS in EDM** | Other missense in codons 268--471 | Flag only; require secondary evidence (TMB >100, SBS10a/b) |
| **Outside EDM** | Any POLE variant outside exons 9--14 | Not classified as POLEmut |

Based on Leon-Castillo et al. 2020 scoring system.

## Development

```bash
# Run tests
source .venv/bin/activate
python -m pytest tests/ -v

# Run with coverage
python -m pytest tests/ --cov=ec_molsubtype --cov-report=term-missing
```

### Project structure

```
src/ec_molsubtype/
  models.py       # Pydantic data models
  classify.py     # WHO 2020 hierarchical classifier
  pole.py         # POLE EDM assessment (tiers, VUS)
  mmr.py          # MMR deficiency (MSI markers + gene mutations)
  tp53.py         # TP53 variant interpretation
  evidence.py     # Secondary evidence aggregation + confidence
  tmb.py          # Tumor mutational burden
  cna.py          # Copy number alteration burden
  signatures.py   # Mutational signature concordance
  io.py           # MAF/VCF parsing
  report.py       # JSON/TSV/human-readable output
  figo.py         # FIGO 2023 molecular staging
  cli.py          # Typer CLI
  data/           # Curated variant lists and signature profiles
```

## Important caveats

1. **Research and clinical decision support only.** Final classification must be confirmed by a molecular pathologist.
2. **MMR status from sequencing is less reliable than IHC.** Always recommend IHC confirmation for MLH1/MSH2/MSH6/PMS2.
3. **POLE VUS require secondary evidence** before classification as POLEmut.
4. **TP53 sequencing is not identical to p53 IHC** (~92% concordance). This tool uses sequencing data as a surrogate.
5. **Signature decomposition from small panels may be unreliable** due to limited variant counts.

## References

- Swedish national guidelines: Nationellt vardprogram livmoderkroppscancer v3.2 (2026-03-10)
- ESGO-ESTRO-ESP guidelines update 2025 (Lancet Oncology)
- FIGO 2023 staging system for endometrial cancer
- Leon-Castillo A, et al. Interpretation of somatic POLE mutations in endometrial carcinoma. *J Pathol*. 2020;250(3):323-335.
- Leon-Castillo A, et al. Clinicopathological and molecular characterisation of multiple-classifier endometrial carcinomas. *J Pathol*. 2020;250(4):413-422.
- Cancer Genome Atlas Research Network. Integrated genomic characterization of endometrial carcinoma. *Nature*. 2013;497(7447):67-73.
- WHO Classification of Female Genital Tumours, 5th edition (2020)

## License

MIT
