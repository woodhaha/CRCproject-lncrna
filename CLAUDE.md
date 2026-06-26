# CLAUDE.md — CRC lncRNA Biomarker Discovery

> TCGA-COAD + GTEx · Machine Learning · Diagnostic & Prognostic lncRNA Signatures
> Structure: MedSci Standard · Updated: 2026-06-25

## Project Summary
Discovery of diagnostic and prognostic long non-coding RNA (lncRNA) biomarkers in colorectal cancer using machine learning algorithms. Uses TCGA-COAD and GTEx (TOIL-recomputed) RNA-seq data.

- **Corresponding Author**: Guanyu Wang (wangguanyu@zju.edu.cn)
- **Institution**: Zhejiang University, Sir Run Run Shaw Hospital
- **Target Journal**: Scientific Reports / Frontiers
- **Analysis Language**: R

## Quick Start
```r
# Set working directory to project root
setwd("D:/Researching/CRCproject")
# Source the optimized pipeline:
source("03_Analysis/pipeline.R")
```

## Project Structure (MedSci Standard)
```
CRCproject/
├── CLAUDE.md                  # This file
├── 01_Literature/             # References & literature
│   ├── PDFs/                  # Downloaded papers
│   └── literature_matrix.md   # (to create)
├── 02_Data/
│   ├── raw/                   # NEVER overwrite
│   ├── cleaned/               # Processed data
│   └── data_dictionary.md     # Variable definitions
├── 03_Analysis/
│   ├── pipeline.R             # Main optimized pipeline
│   ├── R_scripts/             # Helper scripts (20 files)
│   ├── outputs/               # Generated tables/plots
│   └── figures/               # Final figures
├── 04_Manuscript/
│   ├── manuscript.tex         # Main manuscript
│   ├── figures/               # Manuscript figures
│   ├── tables/                # Manuscript tables
│   └── archive/               # Old drafts (2019 versions)
├── 05_Submission/             # Submission package
├── 06_Logs/                   # Decisions & change log
│   ├── decisions.md           # (to create)
│   └── change_log.md          # (to create)
└── docs/                      # Project documentation
```

## Analysis Pipeline
1. **Data Preparation** — Merge TCGA-COAD expression with clinical data; annotate lncRNAs via GENCODE
2. **Differential Expression** — limma: tumor vs normal (TOIL TCGA+GTEx), FDR-adjusted
3. **Survival Analysis** — Univariate Cox PH + log-rank test (parallelized with parLapply)
4. **Feature Selection** — Elastic Net regularized Cox (α grid search via foreach %dopar%)
5. **ML Classification** — Logistic, RF, SVM, Elastic Net, Neural Net (parallel via caret+doParallel)
6. **Survival Prediction** — Cox PH risk score + randomForestSRC
7. **Model Evaluation** — C-index, Brier score, time-dependent ROC/AUC, calibration

## Performance Optimizations (2026-06-25)
- Global PSOCK cluster created once, reused across all parallel sections
- Cox model fitting: `lapply` → `parLapply` (~70% faster)
- CALAUC time series: precomputed risk_df × 1 predict call instead of 562
- Cluster cleanup: `on.exit(stopCluster(.cluster))` ensures no leaks
- Helper scripts: `doMC` → `doParallel` (Windows compatible)
- Legacy script (`R_scripts/CRC data mining.R`) deprecated, points to `pipeline.R`

## Key Results
- **17 feature lncRNAs** selected via Elastic Net Cox
- **7-lncRNA + stage** Cox PH risk score: AUC ~0.725 at 8 years
- **SVM/RF/NN** — high sensitivity/specificity for tumor classification

## R Package Dependencies
```r
# Core: survival, survminer, glmnet, randomForestSRC, caret, limma, doParallel
# Viz: ggplot2, survminer, pROC, survivalROC, ggRandomForests
# Data: data.table, dplyr, tidyr, stringr, table1
```
Install: `install.packages("<pkg>")` or `BiocManager::install("<pkg>")`

## Git Convention
- `feat:` — new analysis or figure
- `fix:` — bug fix
- `docs:` — manuscript or doc updates
- `refactor:` — code reorg
- `data:` — data changes
- `submission:` — submission prep

## Notes
- TOIL expression is log₂(FPKM+0.001) normalized
- Original TCGA-COAD: only 41 non-tumor → GTEx supplements with 308 normal colon
- Parallel backend: PSOCK cluster (n-1 cores), compatible with Windows
- Old `scripts/CRC data mining.R` is deprecated — use `03_Analysis/pipeline.R`
