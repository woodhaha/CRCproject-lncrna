# Manuscript Outline — CRC lncRNA Biomarker Discovery (Revised)

> Based on reproducible pipeline (`pipeline.R`) · Nature-style figures · 2026-06-20

## Target Journal: Frontiers in Oncology / Scientific Reports

---

## Title
**Identification and Validation of Long Non-Coding RNA Signatures for Diagnosis and Prognosis of Colorectal Cancer via Machine Learning**

---

## Abstract (~300 words)
- Background: CRC burden, need for biomarkers, lncRNA promise
- Methods: TOIL TCGA+GTEx (639 samples) + TCGA-COAD clinical (428 patients), limma DEA → Elastic Net Cox → Stepwise Cox → 5 ML classifiers
- Results: 17 lncRNAs (α=0.20), 11 features final Cox (AIC=615.8), top lncRNA MIR210HG (HR=1.40, p=0.0003)
- Classification: SVM 92.1%, NNET 92.1%, RF 91.1% test accuracy
- Prognosis: Risk score AUC 0.725 at 8 years
- Keywords: CRC, lncRNA, machine learning, biomarker, diagnosis, prognosis, elastic net

---

## Introduction (~800 words)
### Paragraph 1: CRC burden & need for biomarkers
- Global statistics, survival by stage, early detection gap
- Rising incidence in young adults

### Paragraph 2: lncRNA biology & biomarker potential
- lncRNA definition, regulatory roles, cancer hallmarks
- Tissue-specific expression → biomarker candidacy
- Known CRC lncRNAs (MIR210HG, LINC00261, PCAT6)

### Paragraph 3: Public data & computational discovery
- TCGA, GTEx, TOIL harmonization
- Limitation of original TCGA-COAD (41 normals)
- GENCODE v22/v23 annotations

### Paragraph 4: Machine learning in CRC
- ML applications in oncology
- Classifiers (SVM, RF, NNET) vs regression (Cox, Elastic Net)
- Gap: dual diagnostic + prognostic lncRNA pipelines are rare

### Paragraph 5: Study objectives & novelty
- Integrated pipeline design
- Reproducible, transparent, fully documented

---

## Results (~3000 words)

### 2.1 Clinical Baseline (Table 1)
| Content | Source |
|---------|--------|
| 428 patients, demographics | CLINICAL_BASELINE.md |
| Stage distribution (I:17.1%, II:39.3%, III:29.4%, IV:14.3%) | pipeline output |
| Dead 22.0%, median OS 676 days | pipeline output |
| Stage is dominant prognostic factor (p<0.001) | Table 1 |
| **Figure**: Fig_S1 (survival by stage) | |

### 2.2 Differential Expression & Feature Selection (Fig 1-2)
| Content | Source |
|---------|--------|
| 639 samples, 60,498 → 15,255 genes after filtering | Part 1 |
| 1,289 DE lncRNAs | |
| Elastic Net: α=0.200, λ=0.230 | Part 3 |
| 17 lncRNAs selected | |
| Feature dot chart | Fig_2 |
| Deviance surface + CV curves | Fig_3a-d |
| **Table**: Top DE lncRNAs with logFC, FDR | |

### 2.3 Survival Analysis & Risk Model (Fig 3-8)
| Content | Source |
|---------|--------|
| Univariate Cox: 799 lncRNAs, 72 significant (p<0.05) | Part 2 |
| Top findings: MIR210HG (HR=1.40), CD27-AS1 (HR=1.80) | |
| Stepwise Cox: 11 features, AIC=615.76, C=0.82 | Part 3 |
| Forest plot of final Cox model | Fig_3 |
| Risk score stratification (log-rank p<0.0001) | Fig_6 |
| Time-dependent ROC (1/3/5-year) | Fig_4 |
| AUC convergence to 0.725 at 8 years | Fig_5 |
| Risk score distribution | Fig_7,8 |

### 2.4 Random Forest Survival (Fig 9-13)
| Content | Source |
|---------|--------|
| Optimal RF: mtry=1, nodesize=2, ntree=500 | Part 4 |
| OOB error: train 25.0%, test 32.2% | |
| Variable importance | Fig_S2, S3 |
| Brier score comparison | Fig_12 |
| C-index comparison | Fig_13 |
| CoxPH outperforms RF for time-dependent prediction | |

### 2.5 Classification Models (Fig 14-16)
| Content | Source |
|---------|--------|
| 5 models: SVM, RF, NNET, Logistic, Elastic Net | Part 5 |
| 8 lncRNA features | |
| Best: SVM 92.1%, NNET 92.1%, RF 91.1% | model_performance.csv |
| Linear models underperform (64-65%) → non-linear boundaries | |
| ROC curves | Fig_14 |
| Lift curves | Fig_15 |
| Model comparison bwplot | Fig_16 |
| **Table 2**: Full performance metrics | |

---

## Discussion (~1500 words)

### Paragraph 1: Summary of findings
- Dual diagnostic + prognostic pipeline
- 17 → 11 features, 8 lncRNAs + stage
- SVM/NNET best classifiers

### Paragraph 2: Biological relevance of top lncRNAs
| lncRNA | Mechanism | Literature |
|--------|-----------|------------|
| MIR210HG | Ferroptosis/PCBP1, miR-1226-3p sponge | Wang 2025, Jiang 2024 |
| CD27-AS1 | Immune-related, STAT3, miR-224-5p/PBX3 | Tao 2021, Ma 2022 |
| EIF3J-AS1 | Autophagy → chemoresistance via ATG14 | Luo 2021 |
| LINC00261 | Tumor suppressor, Slug degradation | Yu 2017 |
| PCAT6 | Pan-cancer ceRNA, 5-FU resistance | Ghafouri-Fard 2021 |

### Paragraph 3: Novel findings
- CTB-25B13.12, AP006621.5, RP11-549B18.1 → no prior functional studies
- These represent genuine discovery candidates

### Paragraph 4: Methodological considerations
- TOIL harmonization advantage
- Elastic Net vs Lasso → correlated features retained
- Stepwise Cox refinement → parsimony
- RF overfitting (OOB gap 7.2%) → small feature set
- CoxPH > RF for this signature size

### Paragraph 5: Comparison with prior work
- Previous lncRNA signatures in CRC
- Machine learning studies in CRC diagnosis
- Reproducibility as a distinguishing feature

### Paragraph 6: Limitations
- Single cohort (TCGA), no external validation
- Missing clinical data (perineural invasion 60% missing)
- Retrospective design
- Need for prospective validation
- Functional validation needed for novel lncRNAs

### Paragraph 7: Future directions
- External validation (GEO cohorts)
- Functional studies on CTB-25B13.12, AP006621.5
- Multi-omics integration
- Clinical nomogram development
- Prospective cohort testing

---

## Methods (~1200 words)

### 4.1 Data Acquisition
- TCGA-COAD HTSeq-FPKM (479 tumor, 41 normal)
- TOIL TCGA+GTEx RSEM FPKM (log₂ normalized)
- GENCODE v22/v23 lncRNA annotations
- Clinical + follow-up from TCGA

### 4.2 Differential Expression Analysis
- limma pipeline, eBayes (trend=TRUE, robust=TRUE)
- BH multiple testing correction
- Filter: FPKM < 0.1 in >60% samples removed

### 4.3 Survival Analysis & Feature Selection
- Univariate Cox PH + Kaplan-Meier log-rank
- Elastic Net Cox: α grid search (0-1 by 0.1), 10-fold CV
- Stepwise Cox: AIC minimization (bidirectional)
- Train/test split: 70:30, stratified by vital status

### 4.4 Random Forest Survival
- randomForestSRC, hyperparameter tuning
- Brier score: riskRegression::Score, bootcv B=100
- C-index: pec::cindex, bootcv B=100

### 4.5 Classification Models
- 5 algorithms: SVM (radial), RF (ranger), NNET, Logistic, Elastic Net
- 10-fold × 5-repeat CV, Box-Cox normalization
- Metrics: Accuracy, AUC, Precision, Recall, F1, Kappa

### 4.6 Reproducibility
- R 4.6.0, seed=12345
- Pipeline available at: pipeline.R
- All figures regenerable with single command

---

## Figures

| Figure | Content |
|--------|---------|
| Fig 1 | Study flowchart (Fig_2_feature_lncRNAs adapted) |
| Fig 2 | Elastic Net tuning (Fig_3a-d combined) |
| Fig 3 | Forest plot of multivariate Cox (Fig_3_forestplot) |
| Fig 4 | Time-dependent ROC (Fig_4_roc panels) |
| Fig 5 | AUC over time (Fig_5_auc_over_time) |
| Fig 6 | KM curves by risk group (Fig_6_survival panels) |
| Fig 7 | ROC curves — 5 classifiers (Fig_14_roc panels) |
| S1 | Survival by stage (Fig_S1) |
| S2 | RF variable importance (Fig_S2, S3) |
| S3 | Feature density/box plots (Fig_S4, S5) |
| S4 | Model tuning profiles (Fig_S6) |

## Tables

| Table | Content |
|-------|---------|
| Table 1 | Clinical baseline (stratified alive/dead, with p-values) |
| Table 2 | Model performance comparison (5 models, train/test, all metrics) |
| Table S1 | 17 candidate lncRNAs: logFC, FDR, HR, CI, p-values |
| Table S2 | Full Cox regression results (799 lncRNAs) |

---

## Evidence Map

| Claim | Evidence | Status |
|-------|----------|--------|
| 17 lncRNAs by Elastic Net | pipeline.R Part 3 | ✅ Reproduced |
| MIR210HG top hit (HR=1.40) | coxph_os.csv | ✅ Verified |
| SVM 92.1% test accuracy | model_performance.csv | ✅ Verified |
| Clinical baseline (428 pts) | CLINICAL_BASELINE.md | ✅ Verified |
| Novel lncRNAs (CTB-25B13.12 etc.) | lncRNA_FUNCTION_REPORT.md | ✅ Literature-confirmed |
| Nature-style figures | outputs/ (36 PDFs) | ✅ Generated |
