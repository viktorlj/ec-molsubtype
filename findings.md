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
