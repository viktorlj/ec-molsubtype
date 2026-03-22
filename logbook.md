# Lab Notebook — ec-molsubtype

## 2026-03-22: GENIE V18 endometrial cancer data extraction and repo setup

**Session context:** Claude Opus 4.6, ~15 min. Initial repo preparation.

### Goal
Prepare the ec-molsubtype repository for GitHub (README, LICENSE, .gitignore) and extract endometrial cancer cases from AACR GENIE V18 MSK-IMPACT panels to serve as a validation dataset for the molecular subtyping algorithm.

### Method
- Explored full codebase: 13 source modules (1,737 LOC), 7 test files (820 LOC), 4 curated JSON data files
- Verified all 85 tests pass
- Explored GENIE V18 data directory (`/Users/viklj600/Bioinformatik/Data/Genie_V18/`): 250K samples, 2.7M mutations, 64-column MAF
- Filtered `data_clinical_sample.txt` for 14 endometrial cancer OncoTree codes on 4 MSK-IMPACT panels (341/410/468/505)
- Extracted matching mutations from `data_mutations_extended.txt` (837 MB) using polars lazy scan
- Computed TMB from non-silent mutation counts / panel coverage size (Mb). Raw TMB values in the GENIE public release are anonymized (near-zero placeholders); only `tmb_bin` categories are informative
- Panel sizes computed from `genomic_information.txt` includeInPanel regions: MSK-IMPACT341=0.932 Mb, MSK-IMPACT410=1.056 Mb, MSK-IMPACT468=1.171 Mb, MSK-IMPACT505=1.261 Mb

### Results
- **4,041 EC samples** on MSK-IMPACT panels (3,955 with at least one mutation)
- **97,474 total mutations** (96,704 non-silent)
- OncoTree distribution: UEC 2,174 (54%), USC 630 (16%), UCS 438 (11%), UCEC 417 (10%), UMEC 146, UCCC 118, UDDC 34, + 7 rarer subtypes
- Panel distribution: MSK-IMPACT505 1,958 (48%), MSK-IMPACT468 1,791 (44%), MSK-IMPACT410 222, MSK-IMPACT341 70
- Key gene mutations: PTEN 3,300 | ARID1A 2,663 | PIK3CA 2,551 | TP53 1,902 | KRAS 790 | POLE 617 | MSH6 472 | MSH2 300 | MLH1 169 | PMS2 92
- Computed TMB: median 5.6, mean 19.9, max 653.4 mut/Mb
- TMB >= 10: 1,131 samples (28.0%); TMB >= 100: 141 samples (3.5%)
- GENIE TMB bins: Mid (2–16) 2,099 | Low (<2) 1,205 | High (>16) 737

### Interpretation
- The 3.5% with TMB >= 100 mut/Mb aligns well with expected POLEmut prevalence in endometrial cancer (~5–8% in TCGA; slightly lower here likely due to threshold differences and panel-based TMB estimation)
- 617 POLE mutations across 4,041 samples — not all will be pathogenic EDMs. Need to filter for exonuclease domain hotspots in the next analysis step
- 1,902 TP53 mutations — USC cases should be enriched for these (serous histology correlates with p53abn)
- MSI status not directly available in GENIE public data (no MSI marker % or MSI calls). MMRd classification will rely on MMR gene mutation patterns unless we can derive MSI from indel signatures
- The high PTEN/ARID1A/PIK3CA mutation rates are consistent with endometrioid histology dominance (54% UEC)

### Decisions & Next Steps
- Run the ec-molsubtype classifier on extracted GENIE data — will need to convert MAF + TMB to the tool's expected input format
- Investigate MMRd detection: without MSI marker data, need biallelic MMR gene inactivation or signature-based approaches
- Cross-validate subtype assignments against known TCGA EC molecular subtypes (for overlapping samples)
- Consider adding CNA data extraction (data_CNA.txt) for secondary evidence scoring
- The extraction script is at `scripts/extract_genie_ec.py`; output in `data/genie/` (gitignored)

### Output files
- `data/genie/ec_msk_samples.tsv` — 4,041 rows, clinical metadata
- `data/genie/ec_msk_mutations.maf` — 97,474 rows, 31 MB
- `data/genie/ec_msk_tmb.tsv` — 4,041 rows, computed TMB + GENIE TMB bins
- `data/genie/extraction_summary.txt` — summary statistics

## 2026-03-22: CNA extraction and TP53-CIN correlation analysis

**Session context:** Claude Opus 4.6, ~5 min. Continues from 2026-03-22: GENIE V18 data extraction.

### Goal
Extract CNA segmentation data for EC MSK-IMPACT samples and quantify the relationship between TP53 mutation status and chromosomal instability (CIN).

### Method
- Filtered `data_cna_hg19.seg` (424 MB, 7.6M segments) for EC sample IDs
- Computed per-sample CIN metrics: FGA (fraction genome altered, |log2R| > 0.2), wGII (weighted genome instability index), segment counts
- Classified TP53 variants: hotspot (curated list) > truncating > missense > other
- Merged with OncoTree codes and computed TMB for joint analysis
- hg19 chromosome sizes used for genome fraction calculations

### Results
- **4,037 / 4,041 samples** have CNA data (197,258 segments)
- **FGA overall:** median 0.044, mean 0.137, max 0.983
- **Segments per sample:** median 35, mean 49, max 323
- **1,709 samples (42%) have TP53 mutations:** 564 hotspot, 665 missense, 444 truncating, 36 other

**TP53 mutant vs wild-type CIN (key finding):**

| Status | n | FGA median | FGA mean | FGA IQR | wGII median |
|--------|---|-----------|----------|---------|-------------|
| TP53 mutated | 1,708 | 0.189 | 0.229 | 0.020-0.384 | 0.1221 |
| TP53 wild-type | 2,329 | 0.009 | 0.069 | 0.000-0.070 | 0.0369 |

**CIN by histological subtype:**
- USC (serous): n=629, FGA median=0.185, TP53mut 90% — highest CIN + TP53 co-occurrence
- UCS (carcinosarcoma): n=438, FGA median=0.328, TP53mut 84% — highest CIN overall
- UEC (endometrioid): n=2,171, FGA median=0.008, TP53mut 16% — low CIN, low TP53
- UCEC (NOS): n=417, FGA median=0.054, TP53mut 53% — intermediate

**CIN by TP53 variant class (among TP53 mutants):**
- Hotspot: n=563, FGA median=0.191
- Missense (non-hotspot): n=665, FGA median=0.205
- Truncating: n=444, FGA median=0.168
- All classes show similar CIN levels; no major difference between hotspot and truncating

### Interpretation
- Strong TP53-CIN correlation confirmed: TP53-mutated samples have ~21x higher median FGA (0.189 vs 0.009). This supports using FGA as secondary evidence for the p53abn subtype.
- The wide IQR for TP53 mutants (0.020-0.384) suggests heterogeneity — some TP53-mutated tumors have low CIN (possibly subclonal TP53 or early events before CIN accumulation).
- USC/UCS showing >84% TP53 mutation rate and high CIN is consistent with the serous-like/p53abn biology.
- UEC low TP53 rate (16%) aligns with endometrioid biology (expected enrichment in POLEmut, MMRd, NSMP subtypes).
- Non-hotspot missense slightly higher FGA than hotspot — possibly because some of these are gain-of-function variants not yet in the hotspot list, or because the hotspot group includes some lower-penetrance variants.

### Decisions & Next Steps
- Use FGA threshold for secondary evidence: p53abn expects FGA >= 0.3, NSMP/POLEmut expect FGA <= 0.2
- Consider the TP53 wild-type outliers with high FGA — could be p53-null by IHC (large deletions not detectable by panel sequencing)
- Integrate CNA metrics into the ec-molsubtype evidence module
- Next: run full classifier on GENIE data

### Output files
- `data/genie/ec_msk_cna_segments.tsv` — 197,258 rows, 11 MB
- `data/genie/ec_msk_cin_metrics.tsv` — 4,037 rows
- `data/genie/ec_msk_tp53_cin.tsv` — 4,037 rows (merged with TP53 status, OncoTree, TMB)
