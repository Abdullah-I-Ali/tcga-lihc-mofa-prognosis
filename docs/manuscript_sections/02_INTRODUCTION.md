# SECTION 2: INTRODUCTION

Hepatocellular carcinoma (HCC) represents approximately 85% to 90% of primary liver malignancies and remains the third leading cause of cancer-related mortality globally [1, 2]. Clinical management and therapeutic allocation are guided predominantly by macro-anatomical staging systems, including the American Joint Committee on Cancer (AJCC) TNM classification and the Barcelona Clinic Liver Cancer (BCLC) algorithm [3, 4]. These staging frameworks rely primarily on macroscopic tumor burden, vascular invasion, extrahepatic metastasis, and underlying hepatic reserve (e.g., Child-Pugh classification or cirrhosis status) [4]. Although macroscopic staging remains indispensable for triage, patients diagnosed within the same clinical stage frequently exhibit markedly divergent survival trajectories, ranging from early postsurgical recurrence to prolonged disease-free survival [5, 6]. This clinical heterogeneity reflects distinct underlying molecular programs that are not fully captured by macroscopic anatomical staging [5, 7].

The advent of high-throughput genomic technologies has enabled multi-layered profiling of human tumors, providing simultaneous views of transcriptomic activation, DNA methylation alterations, somatic copy number variations (CNVs), and single-nucleotide variants (SNVs) [8]. The landmark characterization of HCC by The Cancer Genome Atlas (TCGA) Research Network established that hepatic tumorigenesis involves complex genomic crosstalk, identifying dominant transcriptomic subclasses, large-scale epigenetic silencing, recurrent chromosomal arm gains (such as 1q and 8q) and losses (such as 4q, 13q, and 17p), and characteristic mutations in genes such as *TP53*, *CTNNB1*, and *ALB* [8]. Numerous computational investigations have subsequently sought to harness these multi-omic assays to build prognostic algorithms capable of refining survival risk estimation in HCC [9, 11].

Despite considerable enthusiasm, translation of multi-omics prognostic models into reproducible clinical research tools faces recognized methodological challenges [12, 13]. In high-dimensional settings where candidate biological features far outnumber available patient samples ($p \gg n$), prognostic modeling workflows are particularly vulnerable to optimistic performance estimates when preprocessing, batch correction, feature selection, or other learned transformations incorporate information from outside the training partition [14-16]. When feature filtering or parameter estimation is conducted across an entire dataset prior to partitioning, validation estimates can be substantially inflated [14, 16]. Establishing reproducible performance baselines in hepatocellular carcinoma therefore requires fully encapsulated validation frameworks with strict separation between training and evaluation data [16].

A related challenge concerns the mathematical strategy chosen for multi-modal dimensionality reduction [17]. Many published workflows concatenate raw or arbitrarily filtered feature sets from multiple platforms into a single design matrix before applying penalized regression or tree-based machine learning [17, 18]. This brute-force concatenation frequently allows dense, continuous data types (such as transcriptomics or methylation arrays) to completely overshadow sparser modalities (such as somatic mutation matrices), obscuring the true relative contribution of each genomic layer [17, 19]. Unsupervised latent factor models, such as Multi-Omics Factor Analysis (MOFA2), offer a principled Bayesian framework to decompose multi-modal assays into a shared, low-dimensional latent space regularized under the prior to capture distinct axes of coordinated multi-omic variation without relying on survival labels during factor learning [19, 20]. By modeling view-specific noise distributions (e.g., Gaussian likelihoods for continuous expression and Bernoulli likelihoods for binary mutation indicators), MOFA2 captures coordinated axes of multi-omic variance [19].

However, applying unsupervised factor models to out-of-sample prognostic validation presents an underappreciated technical bottleneck. In standard implementations (such as Bioconductor `MOFA2`), the native model framework is designed primarily for retrospective decomposition; out-of-sample projection functions (e.g., `MOFA2::predict()`) are formulated to reconstruct observed features for samples included in the original training fit rather than estimating latent factor scores for entirely new, unseen test patients [19, 21]. Evaluating an unsupervised MOFA2 survival model on independent cross-validation folds or a held-out test cohort therefore requires establishing an explicit, mathematically sound projection operator that transforms out-of-sample feature vectors into the frozen training latent space while guarding against numerical instability and batch artifacts [22].

Finally, the clinical incremental value of multi-omics survival models remains insufficiently benchmarked [23, 24]. Many genomic studies report model discrimination in isolation, failing to evaluate whether multi-omics signatures provide statistically meaningful prognostic gain over standard, low-cost clinical covariates such as age, gender, AJCC pathologic stage, and histological tumor grade [23, 25]. A clinically relevant assessment therefore requires evaluating whether multi-omics information provides incremental prognostic value beyond established clinical variables within a pre-specified, leakage-controlled validation framework [24, 26].

To address these challenges, we conducted an original model-development study using primary tumor profiles from the TCGA-LIHC cohort. Our primary objective was to construct, validate, and biologically dissect an integrative multi-omics survival model under strict methodological safeguards designed to eliminate optimistic bias. Specifically, we implemented:
1. Strict boundary separation between a verified development cohort ($N=285$) and a quarantined, untouched locked internal holdout cohort ($N=65$) evaluated exactly once.
2. Complete in-fold encapsulation of all data processing—including empirical Bayes batch correction with variance clamping, feature selection across four modalities (11,402 total features), and MOFA2 factorization (15 latent factors)—across a 25-run repeated nested cross-validation protocol.
3. A joint ridge-regularized out-of-sample projection operator ($\lambda = 1.0$) using frozen training loadings, accompanied by rigorous numerical stability diagnostics.
4. Multi-algorithm survival modeling comparing regularized Cox proportional hazards (ElasticNet-Cox, LASSO-Cox), Random Survival Forests, and parametric accelerated failure time models (XGBoost-AFT).
5. Detailed biological characterization of prognostic latent factors using over-representation pathway analysis and cross-validation stability auditing.
6. A pre-specified clinical incremental-value benchmark comparing clinical staging, multi-omics latent factors, and combined models.

This study therefore seeks to establish a methodologically controlled and reproducible performance benchmark for multi-omics prognosis in HCC, providing transparent empirical evidence regarding the degree to which molecular integration adds prognostic information beyond established clinical variables.

---

## Compact Evidence Check: Section 2

| Statement / Literature Claim in Introduction | Primary Source / Citation | Status |
| :--- | :--- | :---: |
| HCC epidemiology: ~85-90% primary liver cancers; 3rd leading cause of cancer death | Sung et al. CA Cancer J Clin 2021 [1]; Llovet et al. Nat Rev Dis Primers 2021 [2] | **VERIFIED** |
| Clinical staging: AJCC TNM and BCLC systems, macro-tumor burden, liver function | Amin et al. AJCC Cancer Staging Manual 8th ed 2017 [3]; Reig et al. J Hepatol 2022 [4] | **VERIFIED** |
| Survival heterogeneity within identical stages | Forner et al. Lancet 2018 [5]; European Association for the Study of the Liver (EASL) J Hepatol 2018 [6] | **VERIFIED** |
| Molecular heterogeneity: distinct subclasses not captured by anatomy | Hoshida et al. Cancer Res 2009 [7]; Boyault et al. Hepatology 2007 | **VERIFIED** |
| TCGA-LIHC multi-omics characterization (RNA, Meth, CNV, SNV; TP53, CTNNB1, ALB) | The Cancer Genome Atlas Research Network. Cell 2017 (PMID: 28622513) [8] | **VERIFIED** |
| Published multi-omics prognostic modeling efforts in HCC | Chaudhary et al. Clin Cancer Res 2018 [9]; Long et al. J Cell Mol Med 2018 [11] | **VERIFIED** |
| Data leakage and optimistic bias in high-dimensional genomic modeling ($p \gg n$) | Simon et al. JNCI 2003 [14]; Ioannidis. Lancet 2005 [15]; Collins et al. Ann Intern Med 2015 (TRIPOD) [16] | **VERIFIED** |
| Methodological vulnerability to optimistic bias from un-split preprocessing ($p \gg n$) | Simon et al. [14]; Collins et al. [16] | **VERIFIED** |
| Multi-modal concatenation limitations vs principled factor models | Picard et al. Comput Struct Biotechnol J 2021 [17]; Cantini et al. Nat Commun 2021 [18] | **VERIFIED** |
| MOFA2 Bayesian framework: likelihoods, unsupervised factor decomposition | Argelaguet et al. Mol Syst Biol 2018 [19]; Argelaguet et al. Genome Biol 2020 [20] | **VERIFIED** |
| Out-of-sample projection limitations in native MOFA2 implementation | `GATE_2B5A_PROJECTION_METHODOLOGY.md`; Argelaguet et al. [19] | **VERIFIED** |
| Need for clinical incremental-value benchmarking (Age, Gender, Stage, Grade) | Steyerberg et al. Epidemiology 2010 [23]; Vickers et al. BMC Med Inform Decis Mak 2008 [24] | **VERIFIED** |
| Study design parameters: Dev $N=285$, Holdout $N=65$, 11,402 features, 15 factors, 25 runs | `TRUE_DEV_285_IDS.csv`; `FINAL_285_NESTED_CV_AUDIT.md` | **VERIFIED** |

### Prohibited Red-Flag Word Check: Section 2
- `first` / `first-ever`: None; framed as "To address these challenges, we conducted an original model-development study...".
- `best` / `superior`: None
- `breakthrough` / `state-of-the-art`: None
- `clinically useful` / `clinically ready` / `ready for bedside`: None; framed as evaluating "whether multi-omics information provides incremental prognostic value beyond established clinical variables".
- `external validation`: Prohibited use checked; text explicitly terms holdout an "untouched locked internal holdout cohort ($N=65$)".
- `invariant` / `causal` / `drives`: None
- `0.6696`: None
- `IBS proves calibration`: None
- `ElasticNet best-performing`: None

### Unresolved Statements in Section 2: None.
Section 2 is fully certified and internally consistent.
