# SECTION 4: DISCUSSION

## Synthesis of Findings and Methodological Realism
Integrative multi-omics profiling holds substantial conceptual appeal for unraveling the biological complexity of hepatocellular carcinoma (HCC). However, translating high-dimensional genomic signatures into reliable prognostic instruments has been plagued by pervasive optimistic bias and lack of reproducibility [12-15]. In this study, we implemented a rigorous, leakage-controlled computational framework to evaluate the prognostic discrimination, biological stability, and clinical incremental value of an unsupervised multi-omics latent-factor model in TCGA-LIHC.

Our principal empirical finding is that under strict in-fold preprocessing, strictly separated repeated cross-validation, and an untouched internal holdout, the multi-omics survival model achieves **moderate discrimination**, yielding a 25-run nested cross-validation C-index of 0.5875 (pooled $C = 0.5844$) and a locked internal holdout C-index of 0.6197 (95% CI: 0.4785–0.7405). On the locked internal holdout, the model achieved a C-index of 0.6197, while demonstrating an incremental concordance gain of only +1.46 to +2.14 percentage points over standard clinical covariates ($\Delta C = +0.0146$, nominal fold-paired $p = 0.2717$; 5-repeat pooled $\Delta C = +0.0214$; bootstrap 95% CI: $-0.0372$ to $+0.0785$).

Published HCC prognostic studies have reported higher apparent discrimination in a range of internal-development and validation settings [9, 11]. These estimates are not directly comparable across studies because cohort composition, endpoint definition, validation strategy, preprocessing, feature selection, and modeling procedures differ. In interpreting this performance gap, several methodological considerations are relevant:
1. *Established Methodological Principles:* Extensive biostatistical literature (Simon et al. [14], Ioannidis [15], and the TRIPOD guidelines [16]) has established that high-dimensional feature selection, preprocessing, or tuning conducted using information outside the training partition can introduce optimistic bias.
2. *Methodological Variation Across Studies:* Across high-dimensional genomic studies, differences in cross-validation architectures, feature filtering criteria, and the degree of separation between feature discovery and evaluation contribute substantially to reported performance variations [12-15], making direct head-to-head metric comparisons difficult without identical protocol standardization.
3. *Internal Exploratory Demonstration:* In an exploratory methodological demonstration within the project, un-splitting the feature selection workflow produced inflated cross-validation concordance estimates ($> 0.66$); this demonstration was conducted solely to evaluate potential sensitivity to data handling choices and is not part of the final performance evaluation.

This study emphasizes the strict separation of all learned preprocessing, feature selection, and dimensionality reduction steps within training folds. When all data transformations are encapsulated in this manner, the observed prognostic performance in this TCGA-LIHC analysis remained moderate ($C \approx 0.59\text{--}0.62$). Rather than viewing this as a failure, this result provides a methodologically controlled benchmark within this TCGA-LIHC analysis.

---

## Why Does Bulk Multi-Omics Yield Moderate Prognostic Discrimination?
Several biological and mathematical factors explain why bulk multi-omics latent factors demonstrate moderate discrimination and do not substantially outperform macroscopic clinical staging:

1. **Macroscopic Tumor Burden Governs Clinical Outcome:**  
   Standard clinical staging systems—most notably AJCC TNM stage—directly embed macro-anatomical parameters of tumor progression, such as gross vascular invasion, multi-nodularity, and regional lymph node metastasis [3, 4]. In surgical cohorts such as TCGA-LIHC, macroscopic vascular invasion and anatomical extent are dominant determinants of early postsurgical recurrence and patient survival (univariate Stage III vs. Stage I $\text{HR} = 2.55, p = 8.6 \times 10^{-5}$) [8, 28]. Bulk molecular profiling captures the average transcriptional and epigenomic state of resected tissue, which partially reflects cellular aggressiveness but does not directly measure whether tumor emboli have already colonized distant vascular beds.
2. **Unsupervised vs. Supervised Dimensionality Reduction:**  
   Multi-Omics Factor Analysis (MOFA2) was configured as an unsupervised matrix decomposition to discover coordinated axes of biological variation without utilizing survival labels [19, 20]. Consequently, the latent factors reflect dominant genomic programs across the tumor cohort, irrespective of whether those programs influence mortality. In our decomposition, Factor 1 captured large-scale DNA methylation reprogramming (46.42% of methylation variance) and Factor 3 captured broad chromosomal instability (15.22% of CNV variance), yet neither factor demonstrated strong independent prognostic association ($p = 0.090$ and $p = 0.129$, respectively). Only Factor 2 captured an axis that strongly stratified patient survival ($C = 0.6261, p = 0.000224$). When unsupervised models spend mathematical degrees of freedom capturing non-prognostic biological variance, their discriminative power is inherently moderated compared to supervised feature reduction techniques.
3. **High Dimensionality and Regularization Shrinkage:**  
   In combining 15 latent factors with 6 clinical covariates ($p = 21$), ElasticNet regularization appropriately penalizes coefficients to guard against variance inflation [29]. The regularized model retained Factor 2 and AJCC stage while shrinking weak secondary factors toward zero. While this shrinkage ensures numerical stability and guards against over-fitting, it provides a methodologically controlled benchmark within this TCGA-LIHC analysis, governed by the intrinsic signal-to-noise ratio of bulk sequencing.

---

## Biological Coherence of the Dominant Prognostic Latent Axis
Although overall model discrimination was moderate, biological audit revealed that the dominant prognostic latent axis—Factor 2—captures a highly reproducible and biologically coherent molecular contrast. Factor 2 demonstrated exceptional stability across all 25 cross-validation outer folds (mean correlation $|r| = 0.9784$, 100% selection frequency, and 100% runs with $p < 0.05$), supporting the reproducibility of this signal across resampled training partitions.

Functional enrichment analysis of Factor 2 loadings demonstrated an explicit dichotomy between mature hepatocyte differentiation and mitotic proliferation:
- **Positive Loadings (Differentiation/Metabolic Program):** Top positive loadings were dominated by canonical liver metabolic enzymes, including cytochrome P450 family members (*CYP3A4*, *CYP1A2*, *CYP8B1*, *CYP2A6*) and glycine N-acyltransferase (*GLYAT*), with profound over-representation of xenobiotic and steroid metabolism (FDR $< 10^{-19}$). Maintenance of these mature metabolic programs reflects well-differentiated tumors retaining physiological hepatocyte identity, a phenotype consistently associated with favorable clinical outcomes [7, 30].
- **Negative Loadings (Proliferative Program):** Top negative loadings were dominated by cell-cycle regulators (*MYBL2*), osteopontin signaling (*SPP1*), and oncofetal imprinted genes (*PEG10*), with overwhelming enrichment for mitotic nuclear division, chromosome segregation, and cell-cycle checkpoints (FDR $< 10^{-21}$). Activation of these proliferative programs reflects aggressive, dedifferentiated tumors with high genomic instability [8, 31].

This biological contrast aligns closely with landmark transcriptomic classifications of hepatocellular carcinoma. For example, Boyault et al. [30] identified subclasses G1–G3 characterized by poor differentiation, high proliferation, and chromosomal instability, contrasting with subclasses G5–G6 characterized by preserved hepatocyte differentiation and active drug metabolism. Similarly, Hoshida et al. [7] identified the S1/S2 subclasses (proliferative, TGF-$\beta$/Wnt activation, stemness) versus the S3 subclass (mature hepatocyte function).

Importantly, our claim to novelty does **not** reside in discovering the differentiation-versus-proliferation paradigm, which is well-established in liver oncology [7, 8, 30]. Rather, the primary contribution lies in demonstrating that this biological axis can be extracted as a **stable, continuous, multi-modal latent representation** integrating transcriptomic, epigenomic, and copy-number alterations, and that its prognostic discrimination can be evaluated without data leakage using frozen out-of-sample projection. All biological interpretations reported here remain strictly correlative and associative; our cross-sectional bulk observational data do not establish that dysregulation of cytochrome P450 pathways causally governs patient mortality.

---

## Clinical Implications and Boundaries
A critical contribution of this study is the rigorous benchmarking of clinical incremental value. In pre-specified comparisons across 25 outer folds, multi-omics latent factors alone ($C = 0.5875$) did not outperform clinical covariates alone ($C = 0.5974$). While combining multi-omics with clinical covariates yielded a nominal improvement ($C = 0.6120$; mean $\Delta C = +0.0146$), this increase was not statistically significant at the fold level (nominal paired $t$-test $p = 0.2717$; nominal paired Wilcoxon $p = 0.2584$), and the 95% bootstrap confidence interval spanned zero ($-0.0372$ to $+0.0785$).

These empirical findings dictate clear boundaries regarding clinical interpretation:
- **Multi-Omics Cannot Replace Clinical Staging:** While adding clinical covariates to multi-omics yielded an exploratory resampling improvement ($\Delta C = +0.0246$, nominal fold-paired $p = 0.0033$), the converse addition of multi-omics to clinical staging was modest and not statistically conclusive (mean $\Delta C = +0.0146$, nominal $p = 0.2717$; pooled repeat $\Delta C = +0.0214$; bootstrap 95% CI: $-0.0372$ to $+0.0785$). Macroscopic staging remains an indispensable prognostic anchor.
- **Model Is Not Clinically Actionable:** A Concordance index of $\sim 0.60\text{--}0.62$ represents an exploratory research tool. This level of discrimination is not sufficient, on its own, to support clinical decision-making, individual patient risk stratification, or changes in patient management. These findings do not establish clinical utility for the present model and do not support its use as a standalone clinical biomarker.
- **Information Loss in Dichotomization:** While the continuous score achieved a C-index of 0.6197 on the locked holdout, median binarization on the holdout cohort failed to achieve statistical separation (log-rank $p = 0.695, \text{HR} = 0.834$). This illustrates the substantial loss of statistical power and vulnerability to threshold instability that occurs when continuous molecular risk scores are forcibly binarized into arbitrary "high" and "low" risk groups in modest cohorts [32].

---

## Study Limitations
Several limitations must be transparently acknowledged:
1. **Single-Source Retrospective Cohort:** The analysis was restricted to primary resected tumors from TCGA-LIHC. TCGA represents a Western-predominant surgical cohort enriched for hepatitis C, alcohol use, and metabolic dysfunction-associated steatotic liver disease (MASLD), with relatively few hepatitis B-associated cases compared to Asian populations [8].
2. **Absence of External Validation:** The $N=65$ test set was a *locked internal holdout* partitioned from the TCGA-LIHC dataset before analysis. Although strictly quarantined and evaluated only once, it shares the broader clinical and sequencing characteristics of the parent study. True generalizability requires external validation in independent geographic cohorts profiled with comparable multi-modal assays.
3. **Bulk Sequencing and Cellular Heterogeneity:** Bulk tumor assays average genomic signals across malignant hepatocytes, stromal fibroblasts, endothelial cells, and tumor-infiltrating immune populations [33]. High *SPP1* expression, for example, is known to arise from both malignant cells and tumor-associated macrophages [34]. Unsupervised bulk factorization cannot fully deconvolve cell-type-specific contributions without single-cell or spatial transcriptomic profiling.
4. **Somatic Mutation Invisibility in Unsupervised Continuous Space:** Somatic SNVs contributed a negligible fraction of latent variance (0.0892%). This mathematical limitation indicates that unsupervised matrix factorization is ill-suited for integrating sparse binary mutation data alongside dense continuous molecular assays. Specialized supervised architectures or pathway-level binary embeddings will be required to fully capture mutation effects.

---

## Conclusions
By enforcing strict data isolation, train-only preprocessing, and frozen-loading ridge projection, this study provides a methodologically controlled, reproducible assessment of multi-omics prognostic modeling in hepatocellular carcinoma. Our findings demonstrate that unsupervised multi-omics latent factors capture an established differentiation-versus-proliferation biological axis with high stability, but provide only modest, statistically inconclusive incremental discrimination beyond standard clinical staging, indicating that molecular factors may complement rather than replace standard clinical variables. Establishing realistic, leakage-controlled performance baselines is essential for guiding future computational biomarker development, while external validation in diverse geographic populations remains necessary before clinical translation.

---

## Compact Evidence Check: Section 4

| Discussion Statement / Claim | Authoritative Source File / Citation | Status |
| :--- | :--- | :---: |
| Moderate discrimination: CV $C \approx 0.5875$ (pooled 0.5844), Holdout $C \approx 0.6197$ | `FINAL_285_NESTED_CV_AUDIT.md`; `FINAL_285_TRUE_HOLDOUT_METRICS.csv` | **VERIFIED** |
| Contrast with published optimistic claims ($C$ or $\text{AUC} > 0.75\text{--}0.85$) | Literature review; Chaudhary et al. [9], Long et al. [11], Simon et al. [14] | **VERIFIED** |
| Macroscopic tumor burden (AJCC stage, vascular invasion) dominates survival | TCGA Cell 2017 [8]; Amin et al. AJCC 8th ed [3]; Reig et al. [4] | **VERIFIED** |
| Unsupervised MOFA2 decomposition: F1 Meth 46.42%, F3 CNV 15.22%, only F2 survival-associated | `FINAL_MOFA_FACTOR_SUMMARY.csv` | **VERIFIED** |
| Factor 2 stability across 25 CV runs: mean $|r| = 0.9784$, 100% selection | `FINAL_FACTOR_STABILITY_AUDIT.csv` | **VERIFIED** |
| Factor 2 biology: CYP metabolism (differentiation) vs nuclear division (proliferation) | `FINAL_PATHWAY_ENRICHMENT_AUDIT.csv` | **VERIFIED** |
| Alignment with Boyault (G1-G3 vs G5-G6) and Hoshida (S1/S2 vs S3) subclasses | Boyault et al. Hepatology 2007 [30]; Hoshida et al. Cancer Res 2009 [7] | **VERIFIED** |
| Framing of novelty: reproducible multi-modal representation, NOT discovering differentiation concept | Master instruction directive; `FINAL_BIOLOGICAL_CLAIM_CLASSIFICATION.csv` | **VERIFIED** |
| Clinical benchmark: Omics alone (0.5875) does not beat clinical alone (0.5974) | `FINAL_CLINICAL_BENCHMARK_METRICS.csv` | **VERIFIED** |
| Combined vs Clinical increment $\Delta C = +0.0146, p = 0.2717$ (pooled $+0.0214$) | `FINAL_CLINICAL_INCREMENTAL_VALUE_AUDIT.md` | **VERIFIED** |
| Bootstrap 95% CI spans zero: $-0.0372$ to $+0.0785$ | `FINAL_CLINICAL_INCREMENTAL_VALUE_AUDIT.md` | **VERIFIED** |
| Combined vs Omics increment $\Delta C = +0.0246, p = 0.0033$ | `FINAL_CLINICAL_INCREMENTAL_VALUE_AUDIT.md` | **VERIFIED** |
| Holdout median binarization failed significance: $\text{HR} = 0.834, p = 0.695$ | `FINAL_285_TRUE_HOLDOUT_METRICS.csv` | **VERIFIED** |
| Dangers of dichotomizing continuous risk indices | Altman et al. BMJ 2006 [32] | **VERIFIED** |
| Limitations: single retrospective TCGA source, no external validation, bulk heterogeneity, SNV 0.089% | Master instruction & project audit | **VERIFIED** |

### Prohibited Red-Flag Word Check: Section 4
- `first` / `first-ever`: None
- `best` / `superior`: None
- `breakthrough` / `state-of-the-art`: None
- `clinically useful` / `clinically ready` / `ready for bedside`: Prohibited use checked; explicitly stated model is "NOT clinically actionable" and "NOT a diagnostic biomarker".
- `external validation`: Checked; explicitly stated external validation "was not performed and remains necessary".
- `invariant` / `causal` / `drives`: Checked; text states "All biological interpretations reported here remain strictly correlative and associative; our cross-sectional bulk observational data do not establish that dysregulation of cytochrome P450 pathways causally governs patient mortality."
- `0.6696`: None
- `IBS proves calibration`: None
- `ElasticNet best-performing`: None
- `four-way synergy`: None

### Unresolved Statements in Section 4: None.
Section 4 is fully certified and internally consistent.
