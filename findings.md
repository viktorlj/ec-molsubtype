# Findings — ec-molsubtype

## 2026-03-22: GENIE V18 EC cohort characteristics

**Evidence:** Extracted 4,041 MSK-IMPACT EC samples from GENIE V18.
**Confidence:** High (large N, well-annotated public dataset).

- Computed TMB (mut/Mb) distribution: median 5.6, mean 19.9, max 653.4
- 3.5% of samples have TMB >= 100 (candidate POLEmut), consistent with literature (~5-8% in TCGA)
- 617 POLE mutations total — majority likely passengers or non-EDM; filtering for pathogenic EDMs is essential
- GENIE public release anonymizes raw TMB; only binned categories (Low/Mid/High) are provided. TMB must be computed from mutation count / panel size
- MSI marker data is NOT available in GENIE public release. MMRd detection must rely on MMR gene mutations or be derived from mutational patterns
- TP53 mutations (n=1,902) should be enriched in USC (serous) cases

**Open questions:**
- What fraction of the 617 POLE mutations are true pathogenic EDMs (tier 1/2 hotspots)?
- Can we reliably detect MMRd without MSI markers using only biallelic MMR gene inactivation?
- How well does panel-computed TMB correlate with GENIE's official TMB bins?

---

## 2026-03-22: TP53 mutations strongly predict chromosomal instability in EC

**Evidence:** 4,037 EC samples with CNA segmentation data (GENIE V18 MSK-IMPACT).
**Confidence:** High (large N, effect size is very large).

- TP53-mutated EC has ~21x higher median FGA than wild-type (0.189 vs 0.009)
- USC (serous): 90% TP53 mutated, FGA median 0.185 — most CIN-high subtype
- UCS (carcinosarcoma): 84% TP53 mutated, FGA median 0.328 — highest absolute CIN
- UEC (endometrioid): 16% TP53 mutated, FGA median 0.008 — CIN-low
- All TP53 variant classes (hotspot, truncating, non-hotspot missense) show similar FGA levels (~0.17-0.21 median), suggesting functional equivalence for CIN association
- Wide IQR in TP53 mutants (0.020-0.384) indicates heterogeneity — some TP53-mutated tumors are CIN-low (possibly subclonal or early events)

**Implications for the classifier:**
- FGA is a strong secondary evidence marker for p53abn subtype
- Recommend FGA >= 0.3 threshold for "high CIN" concordant with p53abn, <= 0.2 for "low CIN" concordant with POLEmut/NSMP
- TP53 wild-type samples with high FGA warrant investigation — possible large TP53 deletions not captured by panel sequencing, or CIN through non-TP53 mechanisms

**Literature support:** Consistent with TCGA 2013 (Nature 497:67-73) showing serous-like/copy-number-high cluster enriched for TP53 mutations.

---

## 2026-03-22: MAF-only MSI detection — literature review and strategy

**Evidence:** Systematic literature review of MSI detection methods from panel variant calls (no BAM access).
**Confidence:** High for published methods; Medium for applicability to MSK-IMPACT panels.

**Key finding 1: BAM-based tools are the gold standard but inapplicable to GENIE MAF data.**
MSIsensor (Middha et al. PMID:30211344) and MiMSI (Ziegler et al. PMID:39746944) both require BAM files. MSIsensor achieves 99.4% concordance with IHC/PCR on MSK-IMPACT data using ~1,000 microsatellite loci. Not usable for our MAF-only analysis.

**Key finding 2: Three published MAF-only approaches exist.**
- MSIpred (Wang & Liang, PMID:30510242): SVM with 22 features from MAF. >= 98% accuracy on WES. Frame_Shift_Del and Frame_Shift_Ins are the most important features. Trained on exome data — panel applicability uncertain.
- I Index (Kim JE et al. PMID:30389464): indel/total mutation ratio >= 9% + mutation load >= 40 detects MSI. Also: >= 5 mutated homopolymer genes. Validated on 382-gene CRC panel.
- Coding microsatellite (cMS) counting: count frameshift indels at known MSI target genes. Idylla assay uses 7 markers (2/7 positive = MSI-H). Can be adapted to any gene panel.

**Key finding 3: ACVR2A (the #1 MSI indicator) is absent from all MSK-IMPACT panels.**
ACVR2A exon 10 A(8) is mutated in 50-91% of MSI-H tumors. Its absence is a significant gap. RPL22 (50% in EC) and BAX (23% in EC) are also absent. However, 13+ other informative cMS genes are well-covered across all panels.

**Key finding 4: Endometrial cancer has a different cMS target spectrum than colorectal cancer.**
TGFBR2 A(10) is mutated in only ~5% of MSI-H EC (vs 80% in gastric, 58% in CRC). EC-enriched targets include JAK1, TFAM, and SMC6 (Cortes-Ciriano et al. PMID:28585546). MSH6 C(8) at 36% is the top mononucleotide repeat target in EC (Kawaguchi et al. PMID:19787250).

**Key finding 5: Indel fraction distinguishes MMRd from POLEmut.**
MMRd: moderate TMB (10-100), high indel fraction (>15-20%), many cMS frameshifts. POLEmut: very high TMB (>100), very low indel fraction (<5%), dominated by C>A substitutions. This distinction is critical for the hierarchical classifier.

**Proposed strategy for ec-molsubtype:**
1. Count frameshift indels at cMS indicator loci covered by the panel (>= 3 = likely MMRd)
2. Compute I index (indel fraction) from MAF
3. Combine with biallelic MMR gene assessment
4. Cross-validate against TMB (10-100 range expected for MMRd)

**Literature:**
- Middha et al. 2017: [DOI](https://doi.org/10.1200/PO.17.00084)
- Ziegler et al. 2025: [DOI](https://doi.org/10.1038/s41467-024-54970-z)
- Wang & Liang 2018 (MSIpred): [DOI](https://doi.org/10.1038/s41598-018-35682-z)
- Kim JE et al. 2018 (I Index): [DOI](https://doi.org/10.1016/j.jmoldx.2018.09.005)
- Cortes-Ciriano et al. 2017: [DOI](https://doi.org/10.1038/ncomms15180)
- Kim TM et al. 2013: [DOI](https://doi.org/10.1016/j.cell.2013.10.015)
- Ballhausen et al. 2020: [DOI](https://doi.org/10.1038/s41467-020-18514-5)
- Kawaguchi et al. 2009: [DOI](https://doi.org/10.3892/ijo_00000411)
- Ferreira et al. 2014 (RPL22): [DOI](https://doi.org/10.1002/humu.22686)
- Hause et al. 2016: [DOI](https://doi.org/10.1038/nm.4191)

**Open questions:**
- What I index threshold is appropriate for MSK-IMPACT 468/505 panels specifically (vs the 382-gene CRC panel where 9% was derived)?
- How many cMS frameshift mutations do MSI-H EC samples in GENIE actually have at covered loci?
- Can we use the cMS counting approach to stratify confidence levels for MMRd calls?

---

## Emerging Hypotheses

### H1: Panel-based molecular subtyping can reproduce TCGA EC classifications with >85% concordance
- **Falsifiable test:** Run classifier on GENIE samples, cross-reference with any available TCGA overlap or published MSK cohort subtype calls
- **Priority:** High — this is the primary validation goal

### H2: MMRd detection from panel data without MSI markers will have lower sensitivity than MSI-based detection
- **Falsifiable test:** Compare MMRd calls using only MMR gene mutations vs. GENIE TMB-high + MMR-mutated as a proxy. If available, compare against published IHC data for MSK EC cohorts
- **Priority:** High — MMRd is the most clinically actionable subtype for immunotherapy

### H3: TP53 wild-type tumors with high CIN represent undetected p53abn cases (large deletions/LOH)
- **Falsifiable test:** Identify TP53-WT samples with FGA > 0.3; check if they cluster with serous histology (USC) and have low TMB. If so, these may be p53-null by IHC but undetectable by panel sequencing
- **Priority:** Medium — affects classifier sensitivity for p53abn subtype

### H4: Coding microsatellite frameshift counting at >= 3 panel-covered loci can identify MMRd with >= 80% sensitivity in GENIE EC data
- **Falsifiable test:** Count cMS frameshifts at TGFBR2/RNF43/MSH6/MSH3/ARID1A/PTEN/JAK1/MRE11A/B2M/AXIN2/KMT2C/APC/TCF7L2 in GENIE EC samples. Compare against TMB 10-100 + high indel fraction as an orthogonal MMRd proxy. Sensitivity should be >= 80% against the proxy.
- **Priority:** High — this is the core question for MAF-only MMRd detection

### H5: The I index (indel/total mutation ratio) threshold for MSI-H on MSK-IMPACT panels is lower than the 9% published for a 382-gene CRC panel
- **Falsifiable test:** Compute I index distributions for GENIE EC samples stratified by TMB bin and cMS count. If the bimodal separation between MSI-H and MSS occurs at a lower threshold (e.g., 5-7%), the panel composition affects the ratio.
- **Priority:** Medium — calibration needed before clinical application
