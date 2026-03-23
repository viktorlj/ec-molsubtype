# ec-molsubtype

Molecular subtyping of endometrial cancer from panel sequencing data, following the **WHO 2020** classification framework, **FIGO 2023** molecular staging, and **Swedish national guidelines** (Nationellt vardprogram livmoderkroppscancer v3.2, 2026-03-10).

Classifies tumors into four prognostic molecular subtypes using a hierarchical sequential algorithm with secondary genomic evidence for confidence scoring. Includes a web interface for clinical use and PDF report generation.

## Validation

Independently validated on two cohorts with ground truth molecular subtypes:

| Cohort | Samples | Accuracy | POLEmut F1 | MMRd F1 | p53abn F1 | NSMP F1 |
|--------|---------|----------|-----------|---------|-----------|---------|
| **TCGA PanCancer** | 507 | **86.8%** | 96.8% | 90.5% | 86.7% | 81.0% |
| **CPTAC** | 95 | **90.5%** | -- | 100% recall | -- | -- |

Additionally validated on **5,141 GENIE V18 samples** across 5 sequencing panels (MSK-IMPACT, UCSF, DFCI), with consistent subtype distributions matching TCGA reference rates.

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

Priority is absolute: a tumor with both a POLE hotspot and a TP53 mutation is classified as **POLEmut**. Multiple-classifier cases (~4% of samples) are reported with specific evidence for each secondary feature but do not change the primary assignment.

### POLE assessment

| Tier | Variants | Action |
|------|----------|--------|
| **Tier 1** (established hotspots) | P286R, V411L, S297F, A456P, S459F | Auto-classify as POLEmut |
| **Tier 2** (confirmed pathogenic) | F367S, L424I, L424V, M444K, D368N, P286H | Auto-classify as POLEmut |
| **VUS in EDM** | Other missense in codons 268--471 | Flag only; require secondary evidence (TMB >100, SBS10a/b) |

Based on Leon-Castillo et al. 2020.

### MMR / MSI assessment

Three independent pathways for MMRd detection, in priority order:

1. **Direct MSI markers** (`msi_pct`): percentage of unstable microsatellite sites from the panel. MSI-H >= threshold (default 20%) -- high confidence.
2. **Coding microsatellite (cMS) proxy**: when direct MSI markers are unavailable, counts frameshift indels at 11 high-specificity indicator genes (JAK1, RNF43, MSH3, MSH6, AXIN2, KMT2C, TCF7L2, APC, MRE11A, B2M, TGFBR2). Boosted by the indel fraction (I-index) to compensate for missing ACVR2A.
3. **Biallelic MMR gene inactivation**: two pathogenic hits in MLH1/MSH2/MSH6/PMS2/EPCAM.

### TP53 assessment

Pathogenic TP53 variants include:
- **Curated hotspots**: R175H, G245S, R248W/Q, R249S, R273H/C, R282W
- **All missense in the DNA-binding domain** (codons 102--292): classified as pathogenic based on TCGA validation showing 82% recall (vs 42% with hotspots alone)
- **Truncating variants**: nonsense, frameshift, splice-site
- **Exclusions**: P72R/R72P benign polymorphism; subclonal variants below 5% VAF

### Secondary evidence

Each classification is scored against concordant secondary features:

| Subtype | Expected TMB | Expected CNA | Substitution spectrum |
|---------|-------------|--------------|----------------------|
| POLEmut | >100 mut/Mb | Low          | C>A >20%, T>G >4%, C>G <0.6%, indels <5% |
| MMRd    | 10--100 mut/Mb | Low--intermediate | Indel fraction >15% |
| p53abn  | <10 mut/Mb  | High (FGA >0.3) | Clock-like (SBS1, SBS5) |
| NSMP    | <10 mut/Mb  | Low          | Clock-like (SBS1, SBS5) |

The substitution spectrum is computed directly from the variant list -- no external signature decomposition tools required.

## Installation

Requires Python >= 3.11 and [uv](https://docs.astral.sh/uv/).

```bash
git clone https://github.com/<your-org>/ec-molsubtype.git
cd ec-molsubtype
uv venv
source .venv/bin/activate
uv pip install -e .
```

## Usage

### Web interface

```bash
ec-molsubtype serve
# Open http://127.0.0.1:8000
```

Upload a MAF file, optionally attach a JSON metadata sidecar, and get an interactive classification report with PDF export. Includes a Methods page documenting the algorithm, limitations, and interpretation guidance. Demo samples are available on the landing page.

### Command line

```bash
# Classify a single sample
ec-molsubtype classify sample.maf --metadata sample_meta.json --output result.json

# Human-readable output
ec-molsubtype classify sample.maf -m meta.json --human

# Classify a batch
ec-molsubtype classify-batch --input-dir ./samples/ --metadata-tsv metadata.tsv --output results.tsv

# Inspect a POLE variant
ec-molsubtype check-pole "p.P286R" --tmb 250.0
```

### Input format

**MAF file** (minimum required columns):

| Column | Required | Description |
|--------|----------|-------------|
| `Hugo_Symbol` | Yes | Gene symbol (e.g., POLE, TP53) |
| `Variant_Classification` | Yes | Effect type (Missense_Mutation, Frame_Shift_Del, etc.) |
| `HGVSp_Short` | Yes | Protein change (e.g., p.P286R) |
| `Reference_Allele` | No | Reference allele (used for spectrum analysis) |
| `Tumor_Seq_Allele2` | No | Alternate allele (used for spectrum analysis) |
| `t_alt_count` | No | Tumor alt allele count (for VAF) |
| `t_ref_count` | No | Tumor ref allele count (for VAF) |

**Sample metadata** (JSON sidecar, optional but recommended):

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

If `msi_pct` is unavailable, the tool automatically computes an MSI proxy from coding microsatellite frameshifts in the MAF.

### Output

JSON per sample with classification, evidence, and clinical notes. Batch output is TSV. PDF reports available from the web interface.

### Demo samples

Seven synthetic samples in `demo/` cover all subtypes and edge cases:

```bash
ec-molsubtype classify demo/polemut_clear.maf -m demo/polemut_clear.json --human
ec-molsubtype classify-batch --input-dir demo/ --metadata-tsv demo/metadata.tsv
```

See `demo/README.md` for descriptions of each case.

## Project structure

```
src/ec_molsubtype/
  classify.py     # WHO 2020 hierarchical classifier
  pole.py         # POLE EDM assessment (tiers, VUS flagging)
  mmr.py          # MMR deficiency (MSI markers + gene mutations)
  tp53.py         # TP53 variant interpretation (hotspots + DBD missense)
  msi.py          # MSI proxy from coding microsatellite frameshifts
  spectrum.py     # Substitution spectrum analysis (SBS10 for POLE, indels for MMRd)
  evidence.py     # Secondary evidence aggregation + confidence scoring
  tmb.py          # Tumor mutational burden
  cna.py          # Copy number alteration burden
  signatures.py   # Pre-computed signature concordance (SigProfiler weights)
  models.py       # Pydantic data models
  io.py           # MAF parsing (file and in-memory)
  report.py       # JSON/TSV/human-readable output
  figo.py         # FIGO 2023 molecular staging annotations
  cli.py          # Typer CLI (classify, classify-batch, check-pole, serve)
  data/           # Curated POLE/TP53/MMR variant lists and signature profiles
  web/            # FastAPI web interface with Jinja2 templates and PDF export
demo/             # 7 synthetic samples (POLEmut, MMRd, p53abn, NSMP, edge cases)
scripts/          # GENIE extraction, TCGA/CPTAC validation, figure generation
figures/          # Publication-grade validation figures (PDF + PNG)
tests/            # 100 tests covering all classification modules
```

## Important caveats

1. **Research and clinical decision support only.** Final classification must be confirmed by a molecular pathologist in the context of histopathology and clinical information.
2. **MMR status from sequencing is less reliable than IHC.** IHC for MLH1/MSH2/MSH6/PMS2 is the clinical standard. Panel sequencing may miss promoter methylation, large deletions, and other mechanisms of MMR loss.
3. **POLE VUS require secondary evidence** (TMB >100, SBS10a/b signature) before classification as POLEmut. The tool flags these but does not auto-classify.
4. **TP53 sequencing does not replace p53 IHC** (~92% concordance). Tumors with TP53 deletions/LOH may appear wild-type by sequencing but show aberrant p53 IHC. This accounts for ~6% of p53abn cases that cannot be detected from mutations alone.
5. **Substitution spectrum from small panels** (<300 genes) should be interpreted with caution. The tool uses 6-class spectrum analysis as a robust alternative to full signature decomposition.
6. **Multiple classifiers** (~4% of cases): some tumors meet criteria for multiple subtypes. Published evidence shows prognosis follows the highest-priority class (e.g., POLEmut+p53abn behaves like POLEmut).

## References

- Swedish national guidelines: Nationellt vardprogram livmoderkroppscancer v3.2 (2026-03-10)
- ESGO-ESTRO-ESP guidelines update 2025 (Lancet Oncology)
- FIGO 2023 staging system for endometrial cancer
- Leon-Castillo A, et al. Interpretation of somatic POLE mutations in endometrial carcinoma. *J Pathol*. 2020;250(3):323-335.
- Leon-Castillo A, et al. Clinicopathological and molecular characterisation of multiple-classifier endometrial carcinomas. *J Pathol*. 2020;250(4):413-422.
- Cancer Genome Atlas Research Network. Integrated genomic characterization of endometrial carcinoma. *Nature*. 2013;497(7447):67-73.
- WHO Classification of Female Genital Tumours, 5th edition (2020)
- Middha S, et al. Reliable pan-cancer microsatellite instability assessment by using targeted next-generation sequencing data. *JCO Precis Oncol*. 2017;1:1-17.
- Cortes-Ciriano I, et al. A molecular portrait of microsatellite instability across multiple cancers. *Nat Commun*. 2017;8:15180.

## License

MIT
