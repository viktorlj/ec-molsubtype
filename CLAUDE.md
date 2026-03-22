# CLAUDE.md — ec-molsubtype

## Project overview

A Python CLI (and eventually web app) for **molecular subtyping of endometrial cancer** from panel sequencing data, following Swedish national guidelines (Nationellt vårdprogram livmoderkroppscancer, version 3.2, 2026-03-10) and the WHO 2020 / FIGO 2023 molecular classification framework.

The tool classifies tumors into four molecular subtypes using a hierarchical sequential algorithm, with an added secondary evidence layer for confidence scoring.

## Clinical background

### The four molecular subtypes (in priority order)

1. **POLEmut** — pathogenic POLE exonuclease domain mutation (exons 9–14). Best prognosis. Swedish guidelines: no adjuvant treatment for stage I–II regardless of other risk factors.
1. **MMRd** — mismatch repair deficient. Intermediate prognosis. Predictive for checkpoint inhibitor benefit. Also triggers Lynch syndrome screening.
1. **p53abn** — TP53 mutation (Swedish guidelines use TP53 mutation, NOT copy number alteration count). Worst prognosis. Swedish guidelines: postoperative chemotherapy for stage I–II regardless of other risk factors.
1. **NSMP** — no specific molecular profile. Diagnosis of exclusion. Intermediate prognosis. ER status is an emerging prognostic sub-stratifier within this group.

### Hierarchical classification algorithm (WHO 2020)

```
Step 1: POLE pathogenic EDM found? → YES → POLEmut (regardless of MMR/TP53 status)
Step 2: MMR deficiency? → YES → MMRd (regardless of TP53 status)
Step 3: TP53 pathogenic mutation? → YES → p53abn
Step 4: None of the above → NSMP
```

Multiple classifiers (~3–6% of cases): report all features, but assign primary class per hierarchy. Evidence shows POLEmut+p53abn behaves like POLEmut, and MMRd+p53abn behaves like MMRd.

### POLE pathogenic variant list

Confirmed pathogenic hotspot EDMs (exonuclease domain, exons 9–14):

- **Tier 1 (established hotspots):** P286R, V411L, S297F, A456P, S459F
- **Tier 2 (additional confirmed pathogenic):** F367S, L424I, L424V, M444K, A456P, D368N, P286H (from León-Castillo et al. scoring system and subsequent studies)
- **POLE VUS in exonuclease domain:** flag but do NOT auto-classify as POLEmut. Require secondary evidence confirmation (TMB >100 mut/Mb, SBS10a/b signature, C>A >20%).
- Mutations OUTSIDE the exonuclease domain are NOT classified as POLEmut.

### TP53 interpretation

Accept pathogenic/likely pathogenic TP53 variants per ClinVar/ACMG criteria. Consider:

- Missense in DNA-binding domain (hotspots: R175H, G245S, R248W/Q, R249S, R273H/C, R282W)
- Truncating variants (nonsense, frameshift, splice-site)
- Filter out: benign polymorphisms (especially P72R), low-VAF subclonal variants below threshold
- Note: in clinical practice, p53 IHC (overexpression or null pattern) is the standard surrogate. TP53 sequencing is used here because the input is panel data.

### MMR assessment from panel data

Two independent sources of MMR evidence are available, and the tool should accept either or both:

**1. MSI marker percentage (primary MSI evidence)**
The local panel includes ~100 MSI markers and reports the percentage of unstable markers. This is a strong direct readout of MSI status.

- **MSI-H:** ≥20% unstable markers (configurable threshold via `--msi-threshold`, default 20%)
- **MSI-L:** >0% but <20% unstable markers
- **MSS:** 0% unstable markers
- Input: a float `msi_pct` (0.0–100.0) in the per-sample metadata
- MSI-H from marker panel is sufficient to call MMRd with high confidence.

**2. MMR gene mutations (supporting / alternative evidence)**

- Pathogenic/likely pathogenic variants in MLH1, MSH2, MSH6, PMS2, EPCAM
- Biallelic inactivation (two hits: mutation + LOH, or two mutations) strongly supports MMRd
- MLH1 promoter hypermethylation (if available from panel)
- Single heterozygous MMR variants alone do NOT confirm somatic MMR deficiency — flag as "possible Lynch carrier, IHC confirmation recommended"

**Combined logic:**

- `msi_pct ≥ threshold` → MMRd (high confidence)
- `msi_pct < threshold` BUT biallelic MMR gene inactivation → MMRd (moderate confidence, flag discordance)
- `msi_pct` not provided AND pathogenic biallelic MMR hits → MMRd (moderate confidence)
- `msi_pct` not provided AND only single heterozygous MMR hit → NOT classified as MMRd; flag as "possible MMRd — IHC recommended"
- No MSI data and no MMR mutations → MMR proficient (for classification purposes)

**Note:** even with MSI marker data, IHC for MLH1/MSH2/MSH6/PMS2 is recommended for Lynch syndrome screening (identifies which protein is lost → guides germline testing).

## Secondary evidence layer

For each classification, compute concordance with secondary genomic features:

|Subtype|Expected TMB                       |Expected CNA               |Expected signatures                                              |Other                                        |
|-------|-----------------------------------|---------------------------|-----------------------------------------------------------------|---------------------------------------------|
|POLEmut|Very high (>100 mut/Mb, often >200)|Low                        |SBS10a, SBS10b dominant; C>A >20%, T>G >4%, C>G <0.6%, indels <5%|                                             |
|MMRd   |High (10–100 mut/Mb)               |Low–intermediate           |SBS6, SBS15, SBS21, SBS26; high indel fraction                   |                                             |
|p53abn |Low (<10 mut/Mb)                   |High (many arm-level SCNAs)|SBS1, SBS5 (clock signatures)                                    |Frequent co-mutations: PIK3CA, PPP2R1A, FBXW7|
|NSMP   |Low (<10 mut/Mb)                   |Low                        |SBS1, SBS5                                                       |Often PTEN, PIK3CA, ARID1A mutated; ER+/PR+  |

### Confidence scoring

- **High:** primary classifier + all concordant secondary evidence
- **Moderate:** primary classifier present but some discordance (e.g., POLE hotspot but TMB only moderately elevated)
- **Low:** borderline cases (POLE VUS, single heterozygous MMR hit, TP53 VUS)
- **Discordant:** primary classifier contradicted by secondary evidence (e.g., "POLE hotspot" but TMB <10 → possible passenger/artifact)

## Technical conventions

### Language & tooling

- **Python only** (≥3.11)
- **uv** for package management (NOT pip, NOT conda)
- **Marimo** notebooks exclusively (NEVER Jupyter). Follow Marimo conventions: functional reactive cells, `mo.md()` for markdown, no global mutable state.
- **polars** for dataframes (NOT pandas)
- **pydantic** for data models and validation
- **typer** for CLI
- **pytest** for testing

### Project structure

```
ec-molsubtype/
├── CLAUDE.md                    # this file
├── pyproject.toml               # uv project, all deps here
├── src/
│   └── ec_molsubtype/
│       ├── __init__.py
│       ├── models.py            # pydantic: SampleInput, ClassificationResult, Evidence
│       ├── classify.py          # main orchestrator: sequential algorithm
│       ├── pole.py              # POLE EDM curation, variant lookup, VUS flagging
│       ├── mmr.py               # MMR gene mutation logic, biallelic assessment
│       ├── tp53.py              # TP53 variant interpretation, hotspot/truncating
│       ├── signatures.py        # mutational signature concordance checks
│       ├── cna.py               # CNA burden scoring (fraction genome altered)
│       ├── tmb.py               # TMB calculation/thresholding
│       ├── evidence.py          # secondary evidence aggregator + confidence scoring
│       ├── io.py                # MAF/VCF parsing (annotated VCF or MAF input)
│       ├── report.py            # structured output: JSON, TSV, human-readable
│       ├── figo.py              # FIGO 2023 molecular staging annotation
│       ├── cli.py               # typer CLI entrypoints
│       └── data/
│           ├── pole_pathogenic.json      # curated POLE EDM list with tiers
│           ├── tp53_pathogenic.json      # TP53 hotspots and functional domains
│           ├── mmr_genes.json            # MMR gene list + known pathogenic variants
│           └── signature_profiles.json   # expected SBS profiles per subtype
├── tests/
│   ├── test_classify.py         # end-to-end classification tests
│   ├── test_pole.py
│   ├── test_mmr.py
│   ├── test_tp53.py
│   ├── test_evidence.py
│   ├── test_io.py
│   ├── conftest.py              # shared fixtures: example samples per subtype
│   └── data/                    # test fixture files (mini MAF/VCF)
├── notebooks/
│   ├── 01_explore_classification.py     # marimo: interactive walkthrough
│   ├── 02_validate_tcga_ec.py           # marimo: validate against TCGA EC cohort
│   └── 03_secondary_evidence.py         # marimo: explore concordance patterns
├── logbook.md                   # experiment log (append-only)
└── findings.md                  # key findings and decisions
```

### CLI interface

```bash
# Classify a single sample
ec-molsubtype classify sample.maf --output result.json

# Classify a batch
ec-molsubtype classify-batch --input-dir ./samples/ --output results.tsv

# Inspect POLE variant
ec-molsubtype check-pole "p.P286R" --tmb 250.0 --signatures sbs_weights.json

# Validate against known TCGA classifications
ec-molsubtype validate --tcga-manifest manifest.tsv --predictions predictions.tsv
```

### Input expectations

**MAF format** (minimum required columns):

- Hugo_Symbol, Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2
- Variant_Classification (Missense_Mutation, Nonsense_Mutation, Frame_Shift_Del, etc.)
- HGVSp_Short (e.g., p.P286R)
- t_alt_count, t_ref_count (for VAF calculation)
- Optionally: ClinVar classification, SIFT, PolyPhen

**Annotated VCF**: expect VEP or SnpEff annotations in INFO/CSQ field.

**Supplementary per-sample metadata** (JSON sidecar or TSV):

- `sample_id`
- `tmb` (mut/Mb) — if not provided, compute from MAF variant count + panel size
- `panel_size_mb` — for TMB calculation
- `msi_pct` (float, 0.0–100.0) — percentage of unstable MSI markers from the panel's ~100 MSI markers. This is the primary input for MSI/MMRd assessment.
- `msi_status` — categorical override (MSI-H / MSS / MSI-L). If both `msi_pct` and `msi_status` are provided, `msi_pct` takes precedence and `msi_status` is used for validation.
- `fraction_genome_altered` — if available from CNA caller
- `signature_weights` — dict of SBS signature name → weight (from SigProfiler/deconstructSigs)

At minimum, either `msi_pct` or MMR gene mutations in the MAF must be present for MMRd assessment. If neither is available, the MMRd step is skipped with a warning.

### Output format

JSON per sample:

```json
{
  "sample_id": "SAMPLE_001",
  "primary_subtype": "POLEmut",
  "confidence": "high",
  "classification_path": [
    {"step": 1, "test": "POLE_EDM", "result": "positive", "variant": "POLE p.P286R", "tier": 1}
  ],
  "multiple_classifier": {
    "is_multiple": true,
    "secondary_features": ["p53abn"],
    "tp53_variant": "TP53 p.R248W"
  },
  "secondary_evidence": {
    "tmb": {"value": 285.3, "expected_range": ">100", "concordant": true},
    "cna_burden": {"value": 0.05, "expected": "low", "concordant": true},
    "signatures": {"SBS10a": 0.45, "SBS10b": 0.22, "concordant": true},
    "substitution_profile": {"C>A_fraction": 0.32, "concordant": true}
  },
  "figo_molecular_annotation": "IAmPOLEmut",
  "clinical_notes": [
    "POLE P286R is a tier 1 established pathogenic hotspot.",
    "Concurrent TP53 R248W detected — multiple classifier. Per hierarchy, classified as POLEmut.",
    "All secondary evidence concordant with POLE ultramutated phenotype.",
    "Swedish guidelines: no adjuvant treatment recommended for stage I-II POLEmut."
  ],
  "flags": [],
  "version": "0.1.0"
}
```

## Development workflow

1. Start with `models.py` — define all pydantic models
1. Build `pole.py`, `mmr.py`, `tp53.py` — the three classifiers with curated data files
1. Implement `classify.py` — the sequential orchestrator
1. Add `io.py` — MAF/VCF parsing
1. Build `evidence.py` + `tmb.py` + `signatures.py` + `cna.py` — secondary evidence
1. Wire up `cli.py`
1. Add `report.py` for output formatting
1. Tests throughout — aim for >90% coverage on classification logic
1. Notebook exploration with TCGA EC data

## Phase 2: Web app (future)

- FastAPI backend reusing the core classification engine
- HTMX + minimal JS frontend (or Marimo as app mode)
- File upload for MAF/VCF
- Interactive result display with evidence visualization
- Batch processing with progress
- PDF report generation

## Key references

- Swedish national guidelines: Nationellt vårdprogram livmoderkroppscancer v3.2 (2026-03-10)
- ESGO-ESTRO-ESP guidelines update 2025 (Lancet Oncology)
- FIGO 2023 staging system for endometrial cancer
- León-Castillo et al. 2020 "Interpretation of somatic POLE mutations in endometrial carcinoma" (J Pathol)
- León-Castillo et al. 2020 "Clinicopathological and molecular characterisation of multiple-classifier endometrial carcinomas" (J Pathol)
- TCGA 2013 endometrial cancer landmark study
- WHO Classification of Female Genital Tumours, 5th edition (2020)

## Important caveats to document in the tool

1. This tool is for **research and clinical decision support only**. Final classification should be confirmed by a molecular pathologist.
1. MMR status from panel sequencing (mutation-only) is LESS reliable than IHC. Always recommend IHC confirmation.
1. POLE VUS require secondary evidence confirmation before classification as POLEmut.
1. The tool does NOT replace p53 IHC — TP53 sequencing and p53 IHC have ~92% concordance but are not identical.
1. Signature decomposition quality depends on the panel size and variant count. Small panels may yield unreliable signature fits.
