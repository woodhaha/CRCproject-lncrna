# CLAUDE.md — CRC lncRNA Biomarker Discovery

> TCGA-COAD + GTEx · Machine Learning · Diagnostic & Prognostic lncRNA Signatures

## Project Summary
Discovery of diagnostic and prognostic long non-coding RNA (lncRNA) biomarkers in colorectal cancer using machine learning algorithms. Uses TCGA-COAD and GTEx (TOIL-recomputed) RNA-seq data to identify lncRNA signatures for tumor classification and survival prediction.

- **Corresponding Author**: Guanyu Wang (wangguanyu@zju.edu.cn)
- **Institution**: Zhejiang University, Sir Run Run Shaw Hospital
- **Target Journal**: Scientific Reports (originally), Frontiers (alternate)
- **Analysis Language**: R

## Quick Start
```r
# In R or RStudio, source the pipeline:
source("scripts/loadpackage.R")       # Load all required packages
source("scripts/CRC data mining.R")   # Main analysis pipeline
```

## Project Structure
```
CRCproject/
├── data/                  # TCGA-COAD expression, clinical, annotations, results
│   ├── TCGA-COAD_Merge.txt              # Merged expression + clinical data
│   ├── COAD_HTSeq_FPKM.txt              # Raw TCGA-COAD expression (FPKM)
│   ├── TcgaTargetGtex_rsem_gene_fpkm     # TOIL-recomputed TCGA+GTEx expression
│   ├── tcga_coad_clinical_followup_data.csv  # Clinical + follow-up
│   ├── gencode.v22.long_noncoding_RNAs.gtf   # GENCODE v22 lncRNA annotations
│   ├── gencode.v23.long_noncoding_RNAs.gtf   # GENCODE v23 lncRNA annotations
│   ├── limma_DEA.txt                    # Differential expression results
│   ├── coxph_os.csv / coxph_os.txt      # Cox regression results
│   ├── log_rank_os.txt                  # Log-rank test results
│   └── Machine_learning_data.csv        # Feature matrix for ML
│
├── scripts/               # R analysis pipeline (20 scripts)
│   ├── loadpackage.R               # Package loader (run first)
│   ├── CRC data mining.R           # Main pipeline
│   ├── cllinical data process.R    # Clinical data preprocessing
│   ├── Bier.R                      # Brier score calculation
│   ├── CALAUC.R / Cindex.R         # Calibration + C-index
│   ├── CaliPlot.R                  # Calibration plots
│   ├── GGsurvplot.R                # Kaplan-Meier plots
│   ├── Plot_model_ROC.R            # Model ROC curves
│   ├── Plot multiple ROC Curve.R   # Multi-model ROC comparison
│   ├── PlotsurvROC.R               # Time-dependent survival ROC
│   ├── PlotRFtune.R                # Random Forest hyperparameter tuning
│   ├── multiplot.R                 # Multi-panel figure assembly
│   ├── RiskScore summary table.R   # Risk score summary
│   ├── Model evaluation.R          # Comprehensive model evaluation
│   ├── Plot_Lift_Curve.R           # Lift curves
│   ├── rndr.R / rndr2.R / rndr3.R  # Rendering/output scripts
│   └── CaliPlot.R                  # Calibration plot helper
│
├── paper/                 # Manuscript and figures
│   ├── manuscript.tex               # Main manuscript (SciRep format)
│   ├── manuscript_frontiers.tex     # Frontiers-format version
│   ├── manuscript_rewritten.tex     # Rewritten version
│   ├── Figure-1.png through Figure-7.png  # Main figures
│   ├── supplementary.pdf            # Supplementary materials
│   ├── table/                       # Tables (xlsx)
│   ├── FIG/                         # Intermediate figure PDFs
│   ├── figs/                        # Combined figures
│   ├── Ref/                         # Reference style files
│   ├── 20191015/                    # First submission draft
│   └── BACKUP/                      # Backup manuscript version
│
└── outputs/               # Generated outputs
```

## Analysis Pipeline
1. **Data Preparation** — Merge TCGA-COAD expression with clinical/follow-up data; annotate lncRNAs via GENCODE v22/v23
2. **Differential Expression** — limma: tumor vs normal, FDR-adjusted
3. **Survival Analysis** — Univariate Cox PH + log-rank test
4. **Feature Selection** — Elastic Net regularized Cox regression (α=0.172, λ=0.267, 10-fold CV)
5. **ML Classification** — SVM, Random Forest, Neural Network for tumor vs normal
6. **Survival Prediction** — Cox PH risk score model + Random Forest survival (randomForestSRC)
7. **Model Evaluation** — C-index, Brier score, time-dependent ROC/AUC, calibration plots

## Key Results
- **17 feature lncRNAs** selected via Elastic Net Cox regression
- **9-lncRNA + stage** Cox PH risk score stratifies patients; AUC ~0.725 at 8 years
- **SVM / RF / Neural Network** classifiers highly sensitive/specific for tumor vs normal
- Key lncRNAs: RP11-440D17.3, EIF3J-AS1, TNRC6C-AS1, MIR210HG, AC113189.5, LINC00261, CTB-25B13.12, AC114730.3

## R Package Dependencies
```r
# Core: survival, survminer, glmnet, randomForestSRC, caret, limma
# Visualization: ggplot2, survminer, pROC, timeROC, risksetROC
# Data: TCGAutils, biomaRt, edgeR, DESeq2
```
Install missing packages before sourcing: `install.packages("<pkg>")` or `BiocManager::install("<pkg>")`

## Environment
- **OS**: Windows
- **R**: RStudio recommended (MCP server available)
- **Shell**: PowerShell (primary), Bash available

## Git Convention
- `feat:` — new analysis or figure
- `fix:` — bug fix in scripts/results
- `docs:` — manuscript or documentation updates
- `refactor:` — code reorganization
- `data:` — data file changes
- `submission:` — submission-related changes

## Notes
- The TOIL-recomputed TCGA+GTEx expression data is already log₂(FPKM+0.001) normalized
- Original TCGA-COAD has only 41 adjacent non-tumor samples — GTEx normal samples supplement this
- Two survival models: Cox PH risk score + Random Forest survival prediction
- Five classification models for diagnostic prediction
- Scripts assume data files are in `data/` — adjust paths if sourcing individually
