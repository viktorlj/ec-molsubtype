# ec-molsubtype: automated molecular subtyping of endometrial cancer from panel sequencing data

## Authors

Viktor Ljungstrom^1^

^1^ Department of Immunology, Genetics and Pathology, Uppsala University, Uppsala, Sweden

## Abstract

Molecular classification of endometrial cancer (EC) into four prognostic subtypes (POLEmut, MMRd, p53abn, NSMP) is now integrated into WHO 2020 and FIGO 2023 guidelines and directly informs treatment decisions. However, translating panel sequencing data into molecular subtypes requires interpreting POLE variants across pathogenicity tiers, assessing mismatch repair deficiency from heterogeneous evidence sources, and evaluating TP53 mutations against curated databases --- tasks that are difficult to standardize across institutions. Here we present ec-molsubtype, an open-source tool that implements the WHO 2020 hierarchical classification algorithm with a secondary evidence layer for confidence scoring. The tool accepts MAF files from any panel, computes a substitution spectrum analysis for POLE signature support, estimates microsatellite instability from coding microsatellite frameshifts when direct MSI markers are unavailable, and integrates immunohistochemistry results when provided. Validated against 507 TCGA PanCancer Atlas samples with ground truth subtypes, ec-molsubtype achieved 86.8% overall accuracy with F1 scores of 96.8% (POLEmut), 90.5% (MMRd), 86.7% (p53abn), and 81.0% (NSMP). Independent replication on 95 CPTAC samples yielded 90.5% accuracy with 100% MMRd recall. Cross-panel validation on 5,141 AACR GENIE V18 samples across five sequencing platforms demonstrated consistent subtype distributions, with endometrioid-specific MMRd rates of 30--32% on MSK-IMPACT panels matching the TCGA reference of ~28--30%. The tool is available as a Python CLI and web application with PDF report generation.

## Introduction

Endometrial cancer is the most common gynecological malignancy in high-income countries, with incidence rising globally (Sung et al. 2018). The landmark 2013 TCGA study identified four distinct molecular subtypes with markedly different prognoses: POLE ultramutated (POLEmut), microsatellite instability hypermutated (MMRd), copy-number high (p53abn), and copy-number low (NSMP) (Cancer Genome Atlas Research Network 2013). This classification was adopted by the WHO in 2020 and integrated into the FIGO 2023 staging system, where it directly modifies clinical management (Concin et al. 2021; Berek et al. 2023). In Sweden, the national clinical guidelines specify that POLEmut stage I--II tumors can omit adjuvant therapy, while p53abn tumors require postoperative chemotherapy regardless of other risk factors.

The clinical implementation of molecular classification relies on a hierarchical sequential algorithm: POLE exonuclease domain mutations take priority, followed by mismatch repair deficiency, then TP53 aberrations, with NSMP as a diagnosis of exclusion. While conceptually simple, the practical interpretation of panel sequencing data presents several challenges. POLE variants must be distinguished across pathogenicity tiers, with variants of uncertain significance (VUS) in the exonuclease domain requiring secondary evidence for classification (Leon-Castillo et al. 2020). MMR deficiency assessment may combine MSI marker data, MMR gene mutation analysis, and immunohistochemistry, each with different sensitivities and specificities. TP53 interpretation must account for hotspot mutations, truncating variants, DNA-binding domain missense variants, and the ~8% discordance between sequencing and p53 immunohistochemistry (Kobel et al. 2020).

Several tools address aspects of this classification, but none provide a complete, panel-agnostic solution with integrated confidence scoring. Commercial platforms embed molecular classification within proprietary pipelines that are not customizable by clinical laboratories. Research tools such as the ProMisE algorithm (Talhouk et al. 2015) were developed for immunohistochemistry-first workflows rather than sequencing-first approaches. There is a need for an open, transparent, and validated tool that accepts standard variant call formats from any sequencing panel and produces clinically actionable molecular subtype assignments with explicit confidence levels.

Here we present ec-molsubtype, a Python tool that implements the full WHO 2020 hierarchical classification with curated variant databases, substitution spectrum analysis, a microsatellite instability proxy for panels lacking direct MSI markers, and optional integration of immunohistochemistry results. We validate the tool against two independent whole-exome cohorts with ground truth subtypes and demonstrate its generalizability across five different targeted sequencing panels.

## Methods

### Classification algorithm

ec-molsubtype implements the WHO 2020 hierarchical sequential classifier in four steps (Figure 1a). At each step, if the criterion is met, the sample receives that subtype and subsequent steps are skipped; otherwise the algorithm proceeds to the next step. All classifiers are evaluated upfront so that multiple-classifier cases (~4% of samples) can be reported even though only one primary subtype is assigned.

**Step 1: POLE assessment.** Variants in the POLE gene are checked against a curated database of pathogenic exonuclease domain mutations (exons 9--14, codons 268--471). Tier 1 hotspots (P286R, V411L, S297F, A456P, S459F) and tier 2 confirmed pathogenic variants (F367S, L424I, L424V, M444K, D368N, P286H) are auto-classified as POLEmut, based on the Leon-Castillo et al. scoring system (Leon-Castillo et al. 2020). Exonuclease domain VUS are flagged but not auto-classified; secondary evidence (TMB >100 mut/Mb, SBS10a/SBS10b signatures) is required.

**Step 2: MMR deficiency assessment.** Three independent evidence sources are evaluated in priority order. First, if direct MSI marker data is available (percentage of unstable microsatellite sites), MSI-H is called at a configurable threshold (default 20%). Second, when direct MSI data is unavailable, a coding microsatellite (cMS) proxy estimates MSI status by counting frameshift indels at 11 high-specificity indicator genes (JAK1, RNF43, MSH3, MSH6, KMT2C, AXIN2, TCF7L2, APC, MRE11A, B2M, TGFBR2), with an indel fraction (I-index) boost to compensate for missing markers. Third, biallelic MMR gene inactivation (two pathogenic hits in MLH1, MSH2, MSH6, PMS2, or EPCAM) directly supports MMRd. When MMR immunohistochemistry is provided, protein loss takes precedence as the clinical gold standard.

**Step 3: TP53 assessment.** TP53 variants are classified as pathogenic if they are curated hotspots (R175H, G245S, R248W/Q, R249S, R273H/C, R282W), missense mutations in the DNA-binding domain (codons 102--292), truncating variants (nonsense, frameshift, splice-site), or classified as pathogenic/likely pathogenic in ClinVar. The P72R/R72P polymorphism is excluded. Low-VAF variants (<5%) are flagged as potentially subclonal. When p53 immunohistochemistry is available and shows an aberrant pattern, this directly classifies as p53abn even without a detectable point mutation, capturing TP53 deletions and LOH.

**Step 4: NSMP.** Samples not meeting any of the above criteria are classified as NSMP (no specific molecular profile).

### Secondary evidence and confidence scoring

Each classification is evaluated against a secondary evidence layer comprising tumor mutational burden (TMB), copy number alteration burden (fraction genome altered, FGA), pre-computed mutational signature weights (if available), and a 6-class substitution spectrum computed directly from the variant list. Expected ranges are: POLEmut (TMB >100, low FGA, C>A >20%, indels <5%), MMRd (TMB 10--100, low-intermediate FGA, indel fraction >15%), p53abn (TMB <10, high FGA >0.3), NSMP (TMB <10, low FGA). Confidence is scored as high (all concordant), moderate (some discordance), low (borderline/VUS), or discordant (primary contradicted).

### MSI proxy from coding microsatellite frameshifts

For panels lacking direct MSI markers, we developed a proxy based on frameshift indels at coding microsatellite loci. We calibrated indicator gene specificity using 4,041 MSK-IMPACT endometrial cancer samples from AACR GENIE V18, comparing frameshift rates in classifier-defined MMRd versus MSS controls (NSMP with TMB <10 and I-index <0.3). Eleven genes showed >99.5% specificity for MSI-H: JAK1 (59.6% sensitivity, 449-fold enrichment), RNF43 (44.7%, 135x), MSH3 (25.7%), MSH6 (17.2%), KMT2C (13.1%), AXIN2 (12.6%), TCF7L2 (10.1%), APC (8.0%), MRE11A (6.2%), B2M (4.6%), and TGFBR2 (2.3%). ARID1A and PTEN were excluded due to high background frameshift rates in MSS endometrioid tumors (11.4% and 8.0%, respectively). The pseudo-MSI percentage is computed as (high-specificity genes with frameshifts / total genes on panel) x 100, with a bonus for elevated indel fraction (I-index >= 0.25) when at least one indicator gene is affected.

### Substitution spectrum analysis

The tool computes the 6-class substitution spectrum (C>A, C>G, C>T, T>A, T>C, T>G) and indel fraction directly from the variant list, following pyrimidine reference convention. This provides panel-compatible secondary evidence without requiring external signature decomposition tools. POLE ultramutated tumors show a distinctive spectrum (C>A >20%, T>G >4%, C>G <0.6%, indels <5%) that is assessed for concordance. MMRd tumors are expected to show high indel fraction (>15%). These thresholds were derived from the Leon-Castillo et al. POLE scoring system and validated on GENIE panel data.

### Validation cohorts

**TCGA PanCancer Atlas** (n=507 with subtypes). Ground truth subtypes from the integrated TCGA analysis: UCEC_POLE (49), UCEC_MSI (148), UCEC_CN_HIGH (163), UCEC_CN_LOW (147). Clinical data, MSIsensor scores, TMB, FGA, and whole-exome mutation data were downloaded from cBioPortal (study ID: ucec_tcga_pan_can_atlas_2018). MSIsensor scores were used directly as the MSI input with a threshold of 10 (Middha et al. 2017).

**CPTAC** (n=95). Independent cohort with subtypes: POLE (7), MSI-H (25), CNV_high (20), CNV_low (43). Categorical MSI status (MSI-H/MSS) and whole-exome mutations were downloaded from cBioPortal (study ID: ucec_cptac_2020).

**AACR GENIE V18** (n=5,141). Panel sequencing data from five platforms: MSK-IMPACT468 (1,757), MSK-IMPACT505 (1,915), UCSF-IDTV5-TO (578), DFCI-ONCOPANEL-3.1 (643), and DFCI-ONCOPANEL-3 (248). No ground truth subtypes are available; subtype distributions were compared against TCGA reference rates. TMB was computed from non-silent mutation counts divided by panel size. FGA was computed from CBS segmentation data (data_cna_hg19.seg).

### Implementation

ec-molsubtype is implemented in Python (>=3.11) using polars for data processing, pydantic for data validation, and FastAPI for the web interface. Curated variant databases (POLE tiers, TP53 hotspots, MMR genes, signature profiles) are stored as JSON. The tool accepts MAF files as input, with optional JSON metadata sidecars for TMB, MSI percentage, FGA, signature weights, and IHC results. Output is available as JSON, TSV, human-readable text, or PDF reports. The web interface provides file upload, interactive result display, and PDF export.

## Results

### Validation against TCGA ground truth subtypes

We first validated ec-molsubtype against 507 TCGA PanCancer Atlas endometrial cancer samples with established molecular subtypes (Figure 1b). The classifier achieved 86.8% overall accuracy (440/507 correct). Per-subtype performance was highest for POLEmut (F1 96.8%, recall 93.9%, precision 100%) and MMRd (F1 90.5%, recall 83.8%, precision 98.4%), reflecting the well-defined molecular features of these subtypes. p53abn achieved F1 86.7% (recall 82.2%, precision 91.8%), and NSMP achieved F1 81.0% (recall 92.5%, precision 72.0%) (Figure 1d).

Error analysis revealed two dominant misclassification patterns. First, 29 ground-truth p53abn samples were classified as NSMP due to TP53 loss of heterozygosity or deletion without a detectable point mutation (median FGA 0.493, confirming copy-number high biology). These cases are inherently undetectable from mutation data alone and would require p53 immunohistochemistry for correct classification. Second, 22 ground-truth MMRd samples were classified as NSMP; all had MSIsensor scores below the threshold of 10 (median 5.3) despite MANTIS scores above 0.4, representing genuinely borderline MSI cases.

A key methodological finding was that classifying all TP53 missense mutations in the DNA-binding domain (codons 102--292) as pathogenic --- rather than only curated hotspots --- improved accuracy from 74.8% to 86.8% (Supplementary Figure 1). This 12-percentage-point gain was driven by 66 p53abn cases with non-hotspot DBD missense variants (at 36 different codons) that were functionally deleterious per the TCGA ground truth.

### Independent replication on CPTAC

The classifier was independently validated on 95 CPTAC endometrial cancer samples (Figure 1c), achieving 90.5% accuracy (86/95 correct). MMRd recall was 100% (25/25), reflecting the availability of categorical MSI status in this cohort. POLEmut recall was 85.7% (6/7), p53abn 70.0% (14/20), and NSMP 95.3% (41/43). The six p53abn misclassifications again involved TP53 deletion/LOH cases without point mutations.

### Cross-panel validation on GENIE V18

To assess generalizability across sequencing platforms, we classified 5,141 endometrial cancer samples from five panels in AACR GENIE V18 (Figure 2). Since no ground truth subtypes are available in GENIE, we compared subtype distributions against TCGA reference rates and evaluated secondary evidence concordance.

Subtype distributions were consistent across panels for endometrioid histology (UEC), the subtype most comparable to the TCGA cohort. MSK-IMPACT468 showed 7.1% POLEmut, 30.8% MMRd, 5.2% p53abn, and 56.9% NSMP among 912 UEC samples; MSK-IMPACT505 showed 5.9%, 32.0%, 4.0%, and 58.1% among 1,083 UEC samples (Figure 2b). The MMRd rates of 30--32% closely match the TCGA reference of ~28--30% for endometrioid tumors, despite relying entirely on the cMS proxy for MSI estimation (GENIE public data lacks direct MSI markers).

The elevated NSMP proportion compared to TCGA (54--61% vs ~39%) and lower p53abn (14--24% vs ~26%) reflect the known limitation that TP53 deletions and LOH are undetectable from mutation data alone, particularly for serous carcinomas. This gap would be addressed by p53 immunohistochemistry, which the tool supports as an optional input.

### Secondary evidence validates subtype assignments

Across all 5,141 GENIE samples, secondary genomic features showed clear separation between subtypes (Figure 3). POLEmut samples had median TMB 123.0 mut/Mb (IQR 70--259) and median indel fraction 0.6%, consistent with the ultramutated, SNV-dominant phenotype. MMRd samples showed median TMB 29.0 (IQR 15--43) with median indel fraction 37.0%, reflecting microsatellite instability-driven indel accumulation. p53abn samples had low TMB (median 4.8) with high FGA (median 0.256), consistent with chromosomal instability. NSMP showed low TMB (5.6) and low FGA (0.060).

The TMB versus indel fraction scatter plot (Figure 3a) demonstrated clean separation of POLEmut (high TMB, very low indel fraction) from MMRd (moderate TMB, high indel fraction), confirming that the substitution spectrum analysis provides robust discrimination even from panel data.

### TP53 mutations and chromosomal instability

Analysis of 4,037 MSK-IMPACT samples with CNA segmentation data confirmed a strong association between TP53 mutation status and chromosomal instability (Figure 3b). TP53-mutated tumors had 21-fold higher median FGA (0.189) compared to wild-type (0.009). This association was consistent across histological subtypes: uterine serous carcinoma (USC) showed 90% TP53 mutation rate with median FGA 0.185, while endometrioid carcinoma (UEC) showed 16% TP53 mutation rate with median FGA 0.008. These findings support FGA as strong secondary evidence for the p53abn subtype and validate its use in confidence scoring.

## Discussion

We present ec-molsubtype, an open-source tool for molecular subtyping of endometrial cancer from panel sequencing data. The tool achieves 86.8% accuracy against TCGA ground truth and 90.5% on the independent CPTAC cohort, with consistent subtype distributions across five different sequencing panels.

Several design decisions warrant discussion. First, classifying all TP53 DNA-binding domain missense variants as pathogenic, rather than restricting to curated hotspots, substantially improved accuracy (+12 percentage points). This approach is consistent with the observation that the vast majority of TP53 missense mutations in the DBD are loss-of-function, as confirmed by functional assays and the IARC TP53 database (Bouaoun et al. 2016). The trade-off is a slight increase in false positive p53abn calls (9 NSMP incorrectly classified as p53abn on TCGA), which is mitigated by recommending p53 immunohistochemistry confirmation.

Second, the coding microsatellite frameshift proxy provides effective MMRd detection without requiring direct MSI markers or BAM-level analysis. The key innovation is specificity-tiered indicator gene selection: excluding ARID1A and PTEN (which have high background frameshift rates in MSS endometrioid tumors) and incorporating the indel fraction as a boosting signal. This achieved endometrioid-specific MMRd rates of 30--32% on MSK-IMPACT panels, matching the TCGA reference. However, sensitivity is inherently lower than direct MSI testing (categorical MSI achieves 100% MMRd recall on CPTAC), supporting the recommendation that direct MSI markers should be used when available.

Third, the tool's modular design allows integration of immunohistochemistry results. MMR protein loss and aberrant p53 staining take precedence over sequencing-based classification, addressing the two main sources of error: TP53 LOH (~6% of p53abn cases) and borderline MSI (~4% of MMRd cases). When both sequencing and IHC are available, the tool reports concordance or discordance, supporting quality assurance in molecular pathology workflows.

The main limitation is the ~13% error rate, dominated by biologically inherent challenges rather than algorithmic failures. TP53 deletions and LOH, which account for approximately 6% of total errors, are fundamentally undetectable from mutation data and require p53 IHC or copy-number analysis. Borderline MSI cases (MSIsensor 5--10, comprising ~4% of errors) reflect the continuous nature of microsatellite instability rather than a discrete classification boundary (Hause et al. 2016). The remaining errors are distributed across rare edge cases.

In comparison to existing approaches, ec-molsubtype offers several advantages: it is open-source and transparent, accepts any panel's MAF output, provides explicit confidence scoring with secondary evidence, handles multiple-classifier cases (~4% of EC), and includes both CLI and web interfaces with PDF report generation. The tool follows Swedish national guidelines and is consistent with ESGO-ESTRO-ESP recommendations.

In conclusion, ec-molsubtype provides a validated, standardized approach to molecular subtyping of endometrial cancer from panel sequencing data. It is suitable for research use and clinical decision support, with the recommendation that final classification be confirmed by a molecular pathologist in the context of histopathology and immunohistochemistry.

## Data availability

The ec-molsubtype source code is available at [GitHub URL]. TCGA data was obtained from cBioPortal (ucec_tcga_pan_can_atlas_2018). CPTAC data was obtained from cBioPortal (ucec_cptac_2020). AACR GENIE V18 data was accessed through the AACR Project GENIE consortium.

## References

Berek JS, et al. FIGO staging of endometrial cancer: 2023. *Int J Gynaecol Obstet*. 2023;162(2):383-394.

Bouaoun L, et al. TP53 variations in human cancers: new lessons from the IARC TP53 database and genomics data. *Hum Mutat*. 2016;37(9):865-876.

Cancer Genome Atlas Research Network. Integrated genomic characterization of endometrial carcinoma. *Nature*. 2013;497(7447):67-73.

Concin N, et al. ESGO/ESTRO/ESP guidelines for the management of patients with endometrial carcinoma. *Int J Gynecol Cancer*. 2021;31(1):12-39.

Cortes-Ciriano I, et al. A molecular portrait of microsatellite instability across multiple cancers. *Nat Commun*. 2017;8:15180.

Hause RJ, et al. Classification and characterization of microsatellite instability across 18 cancer types. *Nat Med*. 2016;22(11):1342-1350.

Kobel M, et al. Interpretation of p53 immunohistochemistry in endometrial carcinomas: toward increased reproducibility. *Int J Gynecol Pathol*. 2019;38(Suppl 1):S123-S131.

Leon-Castillo A, et al. Interpretation of somatic POLE mutations in endometrial carcinoma. *J Pathol*. 2020;250(3):323-335.

Leon-Castillo A, et al. Clinicopathological and molecular characterisation of 'multiple-classifier' endometrial carcinomas. *J Pathol*. 2020;250(4):413-422.

Middha S, et al. Reliable pan-cancer microsatellite instability assessment by using targeted next-generation sequencing data. *JCO Precis Oncol*. 2017;1:1-17.

Sung H, et al. Global cancer statistics 2020: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries. *CA Cancer J Clin*. 2021;71(3):209-249.

Talhouk A, et al. A clinically applicable molecular-based classification for endometrial cancers. *Br J Cancer*. 2015;113(2):299-310.

## Figure legends

**Figure 1. ec-molsubtype validation on independent whole-exome cohorts.** (a) WHO 2020 hierarchical classification algorithm. Samples are tested sequentially for POLE pathogenic exonuclease domain mutations, MMR deficiency, and TP53 pathogenic mutations, with NSMP as the diagnosis of exclusion. (b) Confusion matrix for TCGA PanCancer Atlas validation (n=507). Cell values show count and row-normalized percentage. Overall accuracy 86.8%. (c) Confusion matrix for CPTAC validation (n=95). Overall accuracy 90.5%. (d) Per-subtype F1 scores for TCGA and CPTAC cohorts. (e) Tumor mutational burden distribution by classified subtype across 5,141 GENIE V18 panel-sequenced samples. Horizontal lines at 10 and 100 mut/Mb indicate MMRd and POLEmut expected ranges. (f) Fraction genome altered by classified subtype. Horizontal line at 0.2 indicates the p53abn concordance threshold.

**Figure 2. Cross-panel generalizability.** Subtype proportions across five sequencing panels and the TCGA reference, shown for (a) all endometrial cancer histologies and (b) endometrioid carcinoma (UEC) only. Numbers within bars indicate percentages. The consistent UEC MMRd rates of 30--32% on MSK-IMPACT panels match the TCGA reference despite using the cMS proxy for MSI estimation.

**Figure 3. Secondary evidence validates subtype separation.** (a) TMB versus indel fraction (I-index) colored by classified subtype. Dashed lines indicate the POLEmut zone (TMB >100, indels <5%) and MMRd zone (TMB 10--100, indels >15%). (b) Fraction genome altered by TP53 mutation status, showing 21-fold higher median FGA in TP53-mutated tumors. (c) MSI proxy score (pseudo-MSI percentage) by classified subtype, with the 9% threshold indicated.

**Supplementary Figure 1. Impact of TP53 DNA-binding domain missense reclassification.** (a) Per-subtype recall before (hotspots only) and after (all DBD missense) the reclassification. The +41 percentage point improvement in p53abn recall drove a 12-point increase in overall accuracy. (b) Remaining misclassification patterns after reclassification.
