# Clinical Baseline — TCGA-COAD Cohort

> **Cohort**: 428 colorectal adenocarcinoma patients with complete lncRNA expression + survival data
>
> Source: TCGA-COAD · Extracted: 2026-06-20

---

## Table 1. Baseline Characteristics

| Characteristic | Overall (n=428) | Alive (n=334, 78.0%) | Dead (n=94, 22.0%) | P-value |
|---------------|-----------------|----------------------|--------------------|---------|
| **Demographics** | | | | |
| Age (years), mean ± SD | 66.7 ± 13.0 | 65.8 ± 12.8 | 69.7 ± 13.4 | **0.015** |
| Age, median (range) | 68 (31–90) | 67 (31–90) | 71 (38–88) | — |
| Male, n (%) | 228 (53.3%) | 170 (50.9%) | 58 (61.7%) | 0.081 |
| Female, n (%) | 200 (46.7%) | 164 (49.1%) | 36 (38.3%) | |
| **Survival** | | | | |
| OS (days), mean ± SD | 891 ± 786 | 940 ± 802 | 717 ± 703 | **0.009** |
| OS, median (range) | 676 (6–4,502) | 735 (6–4,502) | 518 (14–3,147) | — |
| OS (years), mean | 2.4 | 2.6 | 2.0 | — |
| **Tumor Stage** | | | | **<0.001** |
| Stage I, n (%) | 73 (17.1%) | 67 (20.1%) | 6 (6.4%) | |
| Stage II, n (%) | 168 (39.3%) | 140 (41.9%) | 28 (29.8%) | |
| Stage III, n (%) | 126 (29.4%) | 96 (28.7%) | 30 (31.9%) | |
| Stage IV, n (%) | 61 (14.3%) | 31 (9.3%) | 30 (31.9%) | |
| **Tumor Marker** | | | | |
| Preoperative CEA (ng/ml), mean | 38.4 | 21.0 | 115.7 | 0.119 |
| Preoperative CEA, median | 3.0 | 2.7 | 6.8 | — |
| CEA available, n | 272 | 218 | 54 | — |
| **Pathological Features** | | | | |
| Venous invasion, n (%) | 89/372 (23.9%) | 69/296 (23.3%) | 20/76 (26.3%) | 0.647 |
| — Missing | 56 | 38 | 18 | |
| Lymphatic invasion, n (%) | 151/387 (39.0%) | 118/306 (38.6%) | 33/81 (40.7%) | 0.750 |
| — Missing | 41 | 28 | 13 | |
| Perineural invasion, n (%) | 43/172 (25.0%) | 33/140 (23.6%) | 10/32 (31.3%) | 0.456 |
| — Missing | 256 | 194 | 62 | |

---

## Key Observations

### Significant Differences (Alive vs Dead)

| Variable | Direction | P-value |
|----------|-----------|---------|
| **Tumor Stage** | Higher stage → higher mortality | **<0.001** |
| **Age** | Dead patients ~4 years older | **0.015** |
| **Overall Survival** | Alive patients followed longer | **0.009** |

### Non-Significant

| Variable | P-value |
|----------|---------|
| Gender | 0.081 |
| Preoperative CEA | 0.119 |
| Venous invasion | 0.647 |
| Lymphatic invasion | 0.508 |
| Perineural invasion | 0.456 |

### Notable

- **High missingness in perineural invasion** (59.8% missing) — limits statistical power
- **CEA is highly skewed** (mean 38.4 vs median 3.0 ng/ml) — a few extreme values drive the mean
- **Stage IV patients are 3× more likely to be in the Dead group** (31.9% vs 9.3%)
- **Stage distribution** is typical for CRC: Stage II > III > I > IV

---

## Statistical Methods

- **Continuous variables**: Independent t-test (reported as mean ± SD)
- **Categorical variables**: Chi-square test (reported as count + percentage)
- **Missing data**: Excluded pairwise; no imputation performed
- **CEA**: Preoperative carcinoembryonic antigen

*Generated from `tcga_coad_clinical_followup_data.csv` via R 4.6.0*
