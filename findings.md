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

## Emerging Hypotheses

### H1: Panel-based molecular subtyping can reproduce TCGA EC classifications with >85% concordance
- **Falsifiable test:** Run classifier on GENIE samples, cross-reference with any available TCGA overlap or published MSK cohort subtype calls
- **Priority:** High — this is the primary validation goal

### H2: MMRd detection from panel data without MSI markers will have lower sensitivity than MSI-based detection
- **Falsifiable test:** Compare MMRd calls using only MMR gene mutations vs. GENIE TMB-high + MMR-mutated as a proxy. If available, compare against published IHC data for MSK EC cohorts
- **Priority:** High — MMRd is the most clinically actionable subtype for immunotherapy
