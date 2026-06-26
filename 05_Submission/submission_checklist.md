# Submission Checklist — CRC lncRNA Biomarker Discovery

**Target Journal:** Scientific Reports (wlscirep)  
**Date:** 2026-06-26  
**Corresponding Author:** Guanyu Wang (wangguanyu@zju.edu.cn)

---

## Manuscript

| # | File | Status |
|---|------|--------|
| 1 | `04_Manuscript/manuscript.tex` (551 lines, 32 refs) | ✅ |
| 2 | `05_Submission/cover_letter.tex` | ✅ |

## Main Figures (7)

| Fig | File | Status |
|-----|------|--------|
| Fig 1 | Pipeline flowchart (Figure-1-Guanyu) | ⬜ custom |
| Fig 2 | `Fig_2_feature_lncRNAs.pdf` + `Fig_3a/b_*.pdf` | ✅ |
| Fig 3 | `Fig_3_forestplot_cox.pdf` | ✅ |
| Fig 4 | `Fig_S1_survival_by_stage.pdf` + `Fig_6_survival_*.pdf` | ✅ |
| Fig 5 | `Fig_4_roc_1/3/5year.pdf` + `Fig_5_auc_over_time.pdf` | ✅ |
| Fig 6 | `Fig_12_brier_*.pdf` + `Fig_13_cindex_*.pdf` | ✅ |
| Fig 7 | `Fig_14_roc_*.pdf` + `Fig_15_lift_*.pdf` + `Fig_16_model_comparison.pdf` | ✅ |

## Main Tables (2)

| Tab | In manuscript | Status |
|-----|-------------|--------|
| Tab 1 | Top 17 lncRNAs (DE + survival) | ✅ |
| Tab 2 | 5-model classification performance | ✅ |

## Supplementary Materials

| # | Description | File | Status |
|---|-------------|------|--------|
| Tab S1 | Clinical baseline (train/test) | In manuscript | ✅ |
| Tab S2 | lncRNA features in patient groups | In manuscript | ✅ |
| Tab S3 | Cross-cohort clinical comparison | `tables/table_supp_clinical_comparison.tex` | ✅ |
| Fig S1 | Survival by stage | `Fig_S1_survival_by_stage.pdf` | ✅ |
| Fig S2-3 | RF variable importance | `Fig_S2/3_rf_vimp*.pdf` | ✅ |
| Fig S4-5 | Feature distribution | `Fig_S4/5_feature_*.pdf` | ✅ |
| Fig S6 | Model tuning | `Fig_S6_model_tuning.pdf` | ✅ |
| Fig EV1-2 | READ external validation | `Fig_EV1/2_READ_*.pdf` | ✅ |
| Fig EV3-4 | GEO external validation | `Fig_EV3/4_GEO_*.pdf` | ✅ |

## Code & Data

| Item | URL |
|------|-----|
| Analysis code | https://github.com/woodhaha/CRCproject-lncrna |
| TOIL expression | https://toil.xenahubs.net |
| TCGA clinical | https://portal.gdc.cancer.gov |
| GEO datasets | GSE39582, GSE17536 (NCBI GEO) |

## Pre-submission

- [x] All numerical claims cross-checked against pipeline outputs (14/14 match)
- [x] All 36 pipeline figures valid (no empty files)
- [x] References cleaned (duplicate RN19 removed, unused RN26 removed)
- [x] Discussion includes external validation + limitations
- [x] Conclusions quantify main model specifications
- [x] Author contributions, competing interests, funding declared
- [x] Code repository public with complete pipeline
- [ ] Compile manuscript.tex → PDF
- [ ] Convert figures to required format if journal specifies EPS/TIFF
- [ ] All authors confirm submission

---

## Submission Package

```
05_Submission/
├── cover_letter.tex
├── manuscript/
│   └── manuscript.tex
├── figures/          (main 7 figures)
├── supplementary/    (S1-S6 + EV1-EV4)
├── tables/           (Table S3)
└── submission_checklist.md
```
