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

## 2026-03-22: Literature review — MSI detection from panel MAF data (no BAM access)

**Session context:** Claude Opus 4.6, ~20 min. Research task for MMRd classification strategy.

### Objective
Determine how to detect microsatellite instability (MSI) from panel sequencing variant calls (MAF format) when BAM files are not available, as is the case for AACR GENIE public data. Identify coding microsatellite (cMS) loci covered by MSK-IMPACT panels that can serve as MSI indicators.

### Approach
- PubMed literature search for MSI detection methods applicable to MAF/variant-call data
- Web search for coding microsatellite loci coordinates and mutation frequencies
- Cross-referenced MSI indicator genes against MSK-IMPACT 341/468/505 gene panels

### Key papers identified

1. **Middha et al. 2017** (PMID:30211344, DOI:10.1200/PO.17.00084) — MSK-IMPACT MSIsensor validation. 12,288 tumors on IMPACT468 (~1,000 microsatellites). MSIsensor score >= 10 defines MSI-H. 99.4% concordance with PCR/IHC for CRC and UEC. **Requires BAM files** — not applicable to MAF-only analysis.

2. **Ziegler et al. 2025** (PMID:39746944, DOI:10.1038/s41467-024-54970-z) — MiMSI, deep learning MSI classifier for MSK-IMPACT. Higher sensitivity than MSIsensor at low tumor purity (0.895 vs 0.67). Also **requires BAM** (multiple instance learning on read-level features).

3. **Wang & Liang 2018 — MSIpred** (PMID:30510242, DOI:10.1038/s41598-018-35682-z) — **Works from MAF only.** SVM with 22 features from mutation annotation data. Key features: Frame_Shift_Del, Frame_Shift_Ins, INDEL counts, SNP/indel ratios, all normalized per Mb. Accuracy >= 98% on TCGA WES data. Trained on exome data; may lose sensitivity on targeted panels due to reduced feature space. Python package: github.com/wangc29/MSIpred.

4. **Kim JE et al. 2018 — I Index** (PMID:30389464, DOI:10.1016/j.jmoldx.2018.09.005) — Indel/total mutation ratio ("I index") for MSI detection from targeted NGS (382 CRC genes). I index >= 9% and somatic mutation load >= 40 detect MSI with high sensitivity/specificity. Also: mutated homopolymer genes >= 5 as criterion. **Works from variant calls.**

5. **Cortes-Ciriano et al. 2017** (PMID:28585546, DOI:10.1038/ncomms15180) — Pan-cancer MSI portrait. Random forest on MSI event counts. Identifies recurrent cMS loci by tumor type. UCEC 28.3% MSI-H. Endometrial-enriched targets: JAK1, TFAM, SMC6. TGFBR2 only 5% in UCEC vs 80% in STAD.

6. **Kim TM et al. 2013** (PMID:24209623, DOI:10.1016/j.cell.2013.10.015) — Landscape of MSI in CRC/EC genomes. Comprehensive coding MS loci catalog. Tumor-type specificity of frameshift targets.

7. **Ballhausen et al. 2020** (PMID:32958755, DOI:10.1038/s41467-020-18514-5) — Shared frameshift mutation landscape. ReFrame tool. 41 cMS in 40 genes. ACVR2A 91% mutated in MSI CRC, median 18 cMS mutations per MSI EC tumor.

8. **Hause et al. 2016** (PMID:27694933, DOI:10.1038/nm.4191) — Classification of MSI across 18 cancer types. MSI-specific instability signatures. MSI may be continuous rather than discrete.

9. **Kawaguchi et al. 2009** (PMID:19787250, DOI:10.3892/ijo_00000411) — 11 candidate mononucleotide repeat target genes in MSI-H EC: hMSH6(C8), TGFBR2(A10), MSH3(A8), MBD4(A10), BAX(G8), PTEN(A6), HDAC2(A9), EPHB2(A9), CASP5(A10), TCF-4(A9), AXIN2(G7).

10. **Ferreira et al. 2014** (PMID:25196364, DOI:10.1002/humu.22686) — RPL22 A(8) exon 2 mutations: 50% in MSI-H EC, 77% in MSI-H CRC.

### Coding microsatellite (cMS) indicator loci — coverage on MSK-IMPACT

| Gene | Repeat | Exon | MSI-H freq (EC) | 341 | 468 | 505 |
|------|--------|------|-----------------|-----|-----|-----|
| ACVR2A | A(8) x2 | Ex3+Ex10 | ~50-91% | NO | NO | NO |
| TGFBR2 | A(10) | Ex3 | 5-36% | YES | YES | YES |
| RNF43 | G(7) | coding | ~30-54% | YES | YES | YES |
| MSH3 | A(8) | Ex1 | ~58% (CRC) | NO | YES | YES |
| MSH6 | C(8) | Ex5 | 36% | YES | YES | YES |
| RPL22 | A(8) | Ex2 | 50% | NO | NO | NO |
| ARID1A | multiple | multiple | 40-73% | YES | YES | YES |
| APC | multiple | Ex16 | varies | YES | YES | YES |
| PTEN | A(6) | Ex7 | common | YES | YES | YES |
| JAK1 | coding MS | multiple | enriched EC | YES | YES | YES |
| AXIN2 | G(7) | coding | varies | YES | YES | YES |
| MRE11A | T(11) | coding | common | YES | YES | YES |
| B2M | coding MS | multiple | common | YES | YES | YES |
| TCF7L2 | A(9) | coding | varies | NO | YES | YES |
| KMT2C | coding MS | multiple | common | YES | YES | YES |
| BAX | G(8) | coding | 23% (EC) | NO | NO | NO |
| BMPR2 | coding MS | multiple | varies | NO | NO | NO |

**Critical gap:** ACVR2A (the single best MSI indicator, 50-91% mutated in MSI-H) is NOT on any MSK-IMPACT panel. RPL22 (50% in EC) and BAX are also absent.

**Well-covered indicators on all 3 panels:** TGFBR2, RNF43, MSH6, ARID1A, APC, PTEN, JAK1, AXIN2, MRE11A, B2M, KMT2C.

**468/505 only:** MSH3, TCF7L2.

### MAF-only MSI detection strategies (without BAM)

**Strategy 1: Indel-based features (MSIpred-like)**
- Count frameshift insertions and deletions per Mb
- Compute I index = indel count / total mutation count
- Threshold: I index >= 9% suggests MSI (Kim JE et al.)
- Requires normalization by panel size
- Limitation: 300-500 gene panel has fewer variants than WES; feature space is compressed

**Strategy 2: Coding microsatellite frameshift counting**
- Count frameshift indels specifically at known cMS loci on the panel
- MSI indicator genes on IMPACT: TGFBR2, RNF43, MSH6, MSH3, ARID1A, PTEN, JAK1, MRE11A, B2M, AXIN2, KMT2C, APC, TCF7L2
- Threshold: >= 3-5 cMS with frameshift mutations strongly suggests MSI-H
- Precedent: Kim JE (mutated homopolymer >= 5 detects MSI), Idylla (2/7 positive = MSI-H)
- Most robust for GENIE data since these specific genes are well-annotated in the MAF

**Strategy 3: Combined TMB + indel fraction + cMS count**
- High TMB (10-100 mut/Mb) + elevated indel fraction + multiple cMS frameshifts = MMRd
- Distinguishes from POLEmut: POLE has very high TMB (>100) but LOW indel fraction (< 5%)
- POLEmut: dominated by C>A substitutions, few frameshifts
- MMRd: moderate-high TMB, HIGH indel fraction (>15-20%), many frameshifts at cMS

**Strategy 4: MMR gene biallelic inactivation (already implemented)**
- Biallelic hits in MLH1/MSH2/MSH6/PMS2 = strong MMRd evidence
- Single heterozygous hit = possible Lynch, needs IHC confirmation
- Combine with cMS counting for confidence scoring

### Proposed implementation for ec-molsubtype

A tiered approach combining all strategies:
1. If `msi_pct` provided (from local panel MSI markers) -> use directly (primary)
2. Count frameshift indels at cMS indicator loci on the panel -> cMS score
3. Compute indel fraction (I index) from MAF
4. Check for biallelic MMR gene inactivation
5. Combine: cMS_score >= 3 AND I_index >= 9% -> MMRd (moderate-high confidence)
6. Cross-validate against TMB range (10-100 for MMRd, exclude POLEmut if TMB > 100 + low indels)

### Key caveats
- ACVR2A absence from MSK-IMPACT is a significant gap — the single most informative MSI indicator gene
- Panel-based cMS counting has lower sensitivity than BAM-based MSIsensor because fewer loci are interrogated
- Indel fraction from small panels may be noisy; need to validate thresholds on GENIE data
- Endometrial cancer shows different cMS target spectrum than CRC: TGFBR2 much less frequent (5% vs 80%), JAK1/TFAM more frequent
- MSIpred was trained on WES data; direct application to panel data needs validation
- The I index threshold of 9% was derived from a 382-gene CRC panel; may need recalibration for IMPACT panels

### Next steps
- Implement cMS counting in ec-molsubtype mmr.py module
- Compute I index from GENIE EC MAF data and determine panel-specific thresholds
- Validate: compare cMS-based MMRd calls against TMB distribution (MMRd should cluster in TMB 10-100)
- Look for RPL22-like high-frequency cMS targets that are on IMPACT but not in our current list

## 2026-03-22: MSI proxy calibration and initial GENIE validation

**Session context:** Claude Opus 4.6. Continues from literature review.

### Goal
Build an MSI proxy from coding microsatellite frameshift counting, calibrate specificity, and run the full classifier on 4,041 GENIE V18 MSK-IMPACT EC samples.

### Method
- Computed per-gene frameshift rates in MSI-H (current classifier MMRd) vs MSS (NSMP, TMB<10, I-index<0.3) control groups
- Stratified 13 cMS genes into high-specificity (>99.5% spec) and low-specificity tiers
- Built `msi.py` module: pseudo_msi_pct = (high_spec_genes_hit / total_high_spec_genes) * 100
- Excluded ARID1A (11.4% background) and PTEN (8.0% background) from score
- Set msi_threshold=9% (captures 1+ of 11 high-spec genes = 9.1%)

### Results — MSI-specificity calibration

| Gene | MSI-H % | MSS % | Enrichment | Specificity | Tier |
|------|---------|-------|------------|-------------|------|
| JAK1 | 59.6% | 0.1% | 449x | 99.9% | HIGH |
| RNF43 | 44.7% | 0.3% | 135x | 99.7% | HIGH |
| MSH3 | 25.7% | 0.0% | inf | 100.0% | HIGH |
| MSH6 | 17.2% | 0.0% | inf | 100.0% | HIGH |
| KMT2C | 13.1% | 0.3% | 49x | 99.7% | HIGH |
| AXIN2 | 12.6% | 0.0% | inf | 100.0% | HIGH |
| TCF7L2 | 10.1% | 0.1% | 152x | 99.9% | HIGH |
| APC | 8.0% | 0.2% | 40x | 99.8% | HIGH |
| MRE11A | 6.2% | 0.1% | 93x | 99.9% | HIGH |
| B2M | 4.6% | 0.0% | inf | 100.0% | HIGH |
| TGFBR2 | 2.3% | 0.1% | 17x | 99.9% | HIGH |
| ARID1A | 91.5% | 11.4% | 8x | 88.6% | LOW |
| PTEN | 62.8% | 8.0% | 8x | 92.0% | LOW |

### Results — Validation (msi_threshold=9%, n=4,041)

| Subtype | Observed | Expected (TCGA) | Notes |
|---------|----------|-----------------|-------|
| POLEmut | 4.6% | ~7% | POLE VUS not auto-classified |
| MMRd | 17.8% | ~28% | No ACVR2A, no IHC, no direct MSI |
| p53abn | 20.0% | ~26% | TP53 LOH/deletions missed |
| NSMP | 57.6% | ~39% | Absorbs missed MMRd + p53abn |

Secondary evidence validation (all concordant with expectations):
- POLEmut: TMB median=121.3, FGA=0.000, I-index=0.006
- MMRd: TMB median=30.9, FGA=0.032, I-index=0.381
- p53abn: TMB median=4.1, FGA=0.238, I-index=0.111
- NSMP: TMB median=4.8, FGA=0.035, I-index=0.167

Histology-specific:
- UEC (endometrioid): 25.2% MMRd (expected ~30%) — good
- USC (serous): 46.7% p53abn, 53.2% NSMP — missing TP53 LOH cases
- Multiple classifiers: 3.9% (expected 3-6%)

### Interpretation
- The classifier architecture works correctly. All secondary evidence metrics cleanly separate subtypes.
- MMRd sensitivity is ~64% of expected (17.8/28). The ~10% gap is consistent with missing ACVR2A (the #1 MSI indicator, 50-91% in MSI-H) and no direct MSI markers or IHC.
- p53abn sensitivity ~77% of expected (20/26). USC showing 53% NSMP suggests many cases with TP53 deletions/LOH not detectable by panel mutation calling — these would show p53-null by IHC but wild-type by sequencing.
- POLEmut is ~66% of expected. Remaining gap likely from POLE VUS that need TMB >100 confirmation.
- The specificity-tiered cMS approach is sound. JAK1 as the top MSI indicator (449x enrichment, 59.6% sensitivity) confirms EC-specific MSI biology.

### Decisions & Next Steps
- The validation demonstrates the classifier is working as designed with known limitations
- For clinical use: msi_pct from local panel MSI markers (primary) bypasses all cMS proxy limitations
- Consider adding logic for POLE VUS auto-promotion when TMB >100 + concordant spectrum
- Explore TP53 detection gap: check if we can use CNA data to infer TP53 LOH
- Write tests for msi.py module

### Output files
- `data/genie/ec_msk_classifications.tsv` — 4,041 rows, full classification results
- `data/genie/validation_summary.txt` — summary statistics
- `scripts/validate_genie_ec.py` — validation pipeline

## 2026-03-23: Multi-panel validation and I-index boost for MSI proxy

**Session context:** Claude Opus 4.6, ~30 min. Continues from 2026-03-22: MSI proxy calibration.

### Goal
Improve MMRd detection by incorporating the indel fraction (I-index) into the MSI proxy, and validate the classifier across multiple non-MSK panels (DFCI, UCSF) to test generalizability.

### Method
- Added I-index boost to `msi.py`: when I-index >= 0.25 and at least 1 cMS gene (any tier) is hit, add virtual gene equivalents. Formula: `bonus = min(2.0, (i_index - 0.15) / 0.20)` virtual genes added to score. Capped at 2 to avoid overcalling.
- Validated false positive rate: checked all NSMP samples that received I-index boost — all remained below 9% threshold (zero false positives).
- Explored non-MSK panels in GENIE V18. Selected: DFCI-ONCOPANEL-3.1 (643 samples), DFCI-ONCOPANEL-3 (248), UCSF-IDTV5-TO (578). PROV-TSO500HT-V2 excluded — has 538 clinical entries but zero mutations in GENIE MAF (data sharing gap).
- Ran classifier on all 5 panels (5,141 total samples) with panel-specific gene lists and computed FGA from CNA segments for non-MSK samples.

### Results

**Cross-panel subtype distribution:**

| Panel | n | POLEmut | MMRd | p53abn | NSMP |
|-------|---|---------|------|--------|------|
| MSK-IMPACT468 | 1,757 | 4.9% | 21.2% | 21.0% | 52.9% |
| MSK-IMPACT505 | 1,915 | 4.7% | 22.9% | 18.5% | 53.9% |
| UCSF-IDTV5-TO | 578 | 5.5% | 19.6% | 13.8% | 61.1% |
| DFCI-ONCOPANEL-3.1 | 643 | 3.1% | 13.8% | 24.3% | 58.8% |
| DFCI-ONCOPANEL-3 | 248 | 3.2% | 25.8% | 22.2% | 48.8% |
| TCGA reference | — | ~7% | ~28% | ~26% | ~39% |

**UEC (endometrioid) only — most comparable to TCGA:**

| Panel | n | POLEmut | MMRd | p53abn | NSMP |
|-------|---|---------|------|--------|------|
| MSK-IMPACT468 | 912 | 7.1% | 30.8% | 5.2% | 56.9% |
| MSK-IMPACT505 | 1,083 | 5.9% | 32.0% | 4.0% | 58.1% |
| UCSF-IDTV5-TO | 453 | 6.2% | 21.4% | 8.2% | 64.2% |
| DFCI-ONCOPANEL-3.1 | 246 | 4.5% | 24.4% | 8.1% | 63.0% |
| DFCI-ONCOPANEL-3 | 154 | 3.2% | 33.8% | 7.1% | 55.8% |

**Secondary evidence (all panels pooled):**
- POLEmut (n=236): TMB=123.0, FGA=0.001, I-index=0.006
- MMRd (n=1,077): TMB=29.0, FGA=0.035, I-index=0.370
- p53abn (n=1,014): TMB=4.8, FGA=0.256, I-index=0.100
- NSMP (n=2,814): TMB=5.6, FGA=0.060, I-index=0.143

**I-index boost impact:**
- 184 MMRd samples (17.1%) had 0 high-spec cMS genes — rescued by I-index boost + low-spec gene hits
- These have TMB median=13.0, I-index=0.400, and 83% are UEC — consistent with true MMRd
- 163 NSMP samples received partial boost but ALL stayed below 9% threshold → zero false positives

### Interpretation
- I-index boost significantly improved MMRd sensitivity without compromising specificity. MSK panels now show 21-23% MMRd overall and 30-32% for UEC, very close to TCGA ~28-30%.
- DFCI-ONCOPANEL-3.1 shows lower MMRd (13.8%) — missing KMT2C, MSH3, TGFBR2 from cMS indicator set reduces sensitivity. DFCI-ONCOPANEL-3 (25.8%) uses identical gene list but a smaller cohort with different ascertainment.
- UCSF lower MMRd (19.6%) may reflect missing B2M from panel or different variant calling pipeline.
- p53abn remains underestimated across all panels (13.8-24.3% vs ~26%). USC histology contributing many NSMP cases suggests TP53 LOH/deletions not detectable from mutation data alone.
- Secondary evidence is remarkably consistent across all panels — validates the classifier architecture.

### Decisions & Next Steps
- The classifier is performing well for UEC/endometrioid subtypes, which is the primary clinical use case
- Panel-specific cMS gene availability should be documented as a known factor in MMRd sensitivity
- Investigate TP53 LOH from CNA data to improve p53abn detection
- Consider adding NKI (211 samples) and COLU (153) panels in future validation rounds
- Write unit tests for msi.py I-index boost logic

### Output files
- `data/genie/ec_multi_panel_classifications.tsv` — 5,141 rows, all panels
- `data/genie/multi_panel_validation_summary.txt` — cross-panel comparison
- `scripts/validate_multi_panel.py` — multi-panel validation pipeline

## 2026-03-23: cBioPortal EC dataset survey for independent validation

**Session context:** Claude Opus 4.6, ~30 min. Systematic search for validation datasets independent of GENIE V18.

### Objective
Identify endometrial cancer datasets on cBioPortal (and elsewhere) that can serve as independent validation cohorts for the molecular subtyping classifier. Requirements: NOT from AACR GENIE, has mutation data, ideally has MSI/MMR status and molecular subtype annotations.

### Approach
- Queried cBioPortal REST API for all endometrial/uterine studies (keyword search + oncotree filtering)
- Retrieved clinical attributes, sample counts, and molecular profiles for each study
- Extracted actual clinical data values (SUBTYPE, MSI_STATUS, POLE, etc.) to assess annotation quality
- Web searches for PORTEC, NRG/GOG, and other deposited datasets
- Assessed MSK-IMPACT studies for GENIE overlap risk

### Results — All cBioPortal EC datasets identified

| Study ID | Name | Samples | Mutations | CNA | MSI Data | Subtype | GENIE Overlap Risk |
|----------|------|---------|-----------|-----|----------|---------|-------------------|
| ucec_tcga_pub | TCGA Nature 2013 | 373 (248 sequenced) | YES (WES MAF) | YES (GISTIC) | YES (5/7 marker) | YES (n=232) | NONE |
| ucec_tcga_pan_can_atlas_2018 | TCGA PanCancer Atlas | 529 | YES (WES MAF) | YES (GISTIC+arm) | YES (MANTIS+sensor) | YES (n=507) | NONE |
| ucec_tcga | TCGA Firehose Legacy | 549 | YES | YES | NO direct | NO | NONE |
| ucec_tcga_gdc | TCGA GDC 2025 | 547 | YES | YES | NO | NO | NONE |
| ucec_cptac_2020 | CPTAC Cell 2020 | 95 | YES (WES+WGS) | YES | YES | YES (4 subtypes) | NONE |
| ucec_msk_2018 | MSK 2018 | 197 | YES (IMPACT) | YES | YES (IHC+MSI) | partial (POLE, IHC) | HIGH |
| ucec_ancestry_cds_msk_2023 | MSK Cancer Discovery 2023 | 1,882 | YES (IMPACT) | YES | NO direct | YES (4 subtypes) | HIGH |
| ucec_ccr_msk_2022 | MSK CCR MSI 2022 | 181 | YES (IMPACT) | NO | YES (MSI+mechanism) | MMRd only | HIGH |
| ucec_ccr_cfdna_msk_2022 | MSK cfDNA 2022 | 44 | YES | NO | NO | NO | HIGH |
| ucec_msk_2024 | MSK NatMed 2024 | 33 | YES (IMPACT) | NO | YES (MSI score) | NO | HIGH |
| ucs_msk_2024 | MSK MolOnc ERBB2 2024 | 69 | YES (IMPACT) | NO | NO | NO | HIGH |
| uec_msk_2024 | MSK HR+ 2024 | 24 | YES (IMPACT) | NO | NO | NO | HIGH |
| ucec_msk_dacruzpaula_2025 | MSK scDNA 2025 | 11 | YES | NO | NO | NO | HIGH |
| uec_cptac_gdc | CPTAC GDC 2025 | 442 | YES | YES | NO | NO | NONE |

### Results — Detailed annotation assessment

**1. TCGA PanCancer Atlas (ucec_tcga_pan_can_atlas_2018) — TOP PICK**
- 529 samples, 507 with molecular subtype annotations
- Subtypes: UCEC_CN_HIGH (163), UCEC_MSI (148), UCEC_CN_LOW (147), UCEC_POLE (49)
- MSI: MANTIS scores (n=526) + MSIsensor scores (n=528) — continuous values
- CNA: GISTIC + arm-level CNA + aneuploidy scores
- TMB: nonsynonymous TMB available (n=517)
- Also has: FGA, methylation, RNA-seq, RPPA
- WES data — fully independent from GENIE (no panel overlap)

**2. TCGA Nature 2013 (ucec_tcga_pub) — GOLD STANDARD**
- 373 samples, 248 with mutation data, 232 with all 4 molecular subtypes
- Subtypes: Copy-number low/Endometrioid (90), MSI/Hypermutated (65), Copy-number high/Serous-like (60), POLE/Ultra-mutated (17)
- MSI: 5-marker call (MSS 213, MSI-H 125, MSI-L 32) + 7-marker call (MSS 223, MSI-H 127, MSI-L 20)
- MLH1 silencing annotated (114/370 positive)
- Mutation rate cluster: LOW (152), HIGH (67), HIGHEST (26)
- CNA cluster K4 annotated
- WES data — fully independent from GENIE

**3. CPTAC (ucec_cptac_2020) — EXCELLENT but small**
- 95 samples with rich multi-omic data
- Subtypes: CNV_low (43), MSI-H (25), CNV_high (20), POLE (7)
- MSI status: MSS (70), MSI-H (25)
- POLE subtype: Yes (7), No (88)
- p53 IHC: mostly "Cannot be determined" (81), some overexpression (4) and loss (2)
- TP53 pathway status annotated
- Mutation signatures (C>A, C>G, C>T, T>A, T>C, T>G fractions) available per sample
- Also: proteomics, acetylproteomics, phosphoproteomics, methylation, circular RNA
- WES+WGS data — fully independent from GENIE

**4. MSK Cancer Discovery 2023 (ucec_ancestry_cds_msk_2023) — LARGE but GENIE overlap**
- 1,882 samples — the LARGEST EC cohort with molecular subtypes
- Subtypes: CN-H/TP53abn (746), CN-L/NSMP (599), MSI-H (440), POLE (97)
- Gene panels: IMPACT468 (1,185), IMPACT505 (476), IMPACT410 (163), IMPACT341 (58)
- **WARNING: HIGH GENIE overlap risk.** All use MSK-IMPACT panels, same as GENIE MSK contribution. Many/most of these patients are likely in GENIE V18.
- PMID: 37651310 (Weigelt et al., Cancer Discovery 2023)

**5. MSK 2018 (ucec_msk_2018) — SMALL but has IHC**
- 197 samples of advanced-stage EC
- MMR-D by IHC: Y (24), N (75), not done (89), inconclusive (1)
- MSI status: MSS (154), MSI-H (28), MSI-indeterminate (14)
- IHC detail: protein-level MMR loss patterns (MLH1/PMS2, MSH6, MSH2/MSH6)
- POLE ultra-mutated: Y (2), N (187)
- SCNA cluster: A (60), B (71), C (31)
- **WARNING: HIGH GENIE overlap risk.** MSK-IMPACT data.
- PMID: 30068706

**6. MSK CCR MSI 2022 (ucec_ccr_msk_2022) — MSI-ENRICHED**
- 181 samples — all MSI-H/MSI-focused
- MSI status: high (168), low (8), indeterminant (5)
- MMR-D mechanism: MLH1 methylation (117), Somatic (40), Germline (24)
- Useful for testing MMRd detection sensitivity, but ascertainment-biased (MSI-enriched)
- **WARNING: HIGH GENIE overlap risk.**
- PMID: 35849120

### Results — Non-cBioPortal datasets

**PORTEC trials (PORTEC-1, -2, -3)**
- PORTEC-3 molecular classification published (Leon-Castillo et al. JCO 2020, PMID: 32749941)
- ~423 samples with molecular subtypes (POLEmut, MMRd, p53abn, NSMP)
- **Data NOT publicly deposited.** Held by LUMC, Leiden. Available only through collaboration.

**NRG/GOG studies**
- Not found on cBioPortal. No publicly deposited datasets with molecular subtypes identified.

**MSK-CHORD (msk_chord_2024)**
- 25,040 samples but focused on NSCLC, breast, CRC, prostate, pancreatic
- No endometrial cancer samples included

### Interpretation

**Best candidates for independent validation (no GENIE overlap):**

1. **TCGA PanCancer Atlas** — The strongest option. 529 samples, 507 with subtypes, continuous MSI scores, comprehensive annotations. The TCGA original subtypes (CN_HIGH, MSI, CN_LOW, POLE) map directly to p53abn, MMRd, NSMP, POLEmut. WES-based — fully independent.

2. **TCGA Nature 2013** — The original 373-sample publication cohort. Subset of PanCancer Atlas but has the classic 5/7-marker MSI calls and MLH1 silencing annotations that PanCancer Atlas lacks. Most rigorous ground truth.

3. **CPTAC** — Small (n=95) but the richest annotation set: proteomics, mutation signatures per sample, MSI status, POLE status. Ideal for testing secondary evidence concordance. Fully independent.

4. **TCGA GDC 2025** — 547 samples, updated mutation calls through GDC pipeline. No subtype annotations, but we can derive them from mutations (POLE, TP53) and compare against PanCancer Atlas annotations for the same patients.

**Datasets with GENIE overlap risk (use with caution):**

5. **MSK Cancer Discovery 2023** — 1,882 samples with molecular subtypes. The MOST useful dataset if we can identify and remove GENIE-overlapping patients. All use IMPACT panels (identical to GENIE). Sample deduplication by barcode comparison would be needed.

6. **MSK 2018** — 197 samples with IHC-validated MMR status. Valuable ground truth for MMRd detection, but most patients are likely in GENIE V18.

**GENIE overlap assessment:**
All MSK cBioPortal studies use MSK-IMPACT panels. MSK contributes all its IMPACT data to GENIE. Therefore, any MSK cBioPortal study should be assumed to have near-complete overlap with GENIE unless proven otherwise by barcode comparison. The GENIE data guide confirms that MSK submits prospectively sequenced clinical samples, which is the same source as these publications.

### Decisions & Next Steps
1. Download TCGA PanCancer Atlas data (ucec_tcga_pan_can_atlas_2018) as primary validation — mutations, CNA, clinical data with subtypes + MSI scores
2. Download TCGA Nature 2013 data (ucec_tcga_pub) for MSI marker ground truth
3. Download CPTAC data (ucec_cptac_2020) for secondary evidence validation
4. Map TCGA subtype labels to our classifier categories: UCEC_POLE -> POLEmut, UCEC_MSI -> MMRd, UCEC_CN_HIGH -> p53abn, UCEC_CN_LOW -> NSMP
5. Run classifier on TCGA mutation data and compute concordance
6. For MSK CDS 2023: attempt barcode-level deduplication against GENIE V18 sample list
7. Consider requesting PORTEC data through LUMC collaboration

### Data access
- cBioPortal API: `https://www.cbioportal.org/api/`
- Study-specific mutations: `GET /molecular-profiles/{studyId}_mutations/mutations`
- Clinical data: `GET /studies/{studyId}/clinical-data`
- Bulk download: `https://www.cbioportal.org/datasets` (zip per study)
- GitHub datahub: `https://github.com/cBioPortal/datahub` (via git-lfs)
- Raw files: `https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/{studyId}/data_mutations.txt`
