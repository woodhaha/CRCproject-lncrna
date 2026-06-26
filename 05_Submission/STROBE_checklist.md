# STROBE Checklist Audit — Colorectal Cancer lncRNA Study
**Study design**: Cross-sectional / Cohort (TCGA-COAD, TCGA--GTEx colon data)
**Manuscript**: `manuscript_rewritten.tex` — v20 pipeline
**Review date**: 2026-06-25

---

## STROBE 22-Item Checklist

### Title and Abstract (Item 1)
| Item | Status | Notes |
|------|--------|-------|
| 1a — Title indicates study design | ⚠️ PARTIAL | Title states "Machine Learning-Based Identification" but does not explicitly say "cross-sectional study" or "cohort study." The design is inferential from context. |
| 1b — Abstract structured, balanced | ✅ PRESENT | Abstract includes background, methods, results (with key numbers: 37 lncRNAs, 443 patients, 22 features, AUC=0.725, accuracies), and conclusion. |

### Introduction (Items 2-3)
| Item | Status | Notes |
|------|--------|-------|
| 2 — Scientific background and rationale | ✅ PRESENT | Background covers CRC burden, lncRNA biomarker potential, TCGA/GTEx/TOIL resources, and ML in CRC research (paragraphs 1-4). Well-referenced with citations. |
| 3 — Specific objectives, hypotheses | ✅ PRESENT | Last paragraph of Introduction: "identify lncRNA signatures with dual diagnostic and prognostic value in CRC." Objectives clearly stated. |

### Methods (Items 4-12)
| Item | Status | Notes |
|------|--------|-------|
| 4 — Study design, setting, data sources | ✅ PRESENT | TCGA-COAD and TCGA--GTEx described in detail (UCSC Xena, GDC Data Portal). TOIL-recomputed data, GENCODE annotations, sample counts provided. |
| 5 — Participants: eligibility criteria | ⚠️ PARTIAL | TCGA-COAD and GTEx colon samples described, but no explicit eligibility/inclusion/exclusion criteria are stated. Gene filtering thresholds are given (FPKM < 0.1 in >60% of samples). |
| 6a — Variables: clearly defined outcomes | ✅ PRESENT | Overall survival (OS) is the endpoint. Diagnostic endpoint is tumor vs. non-tumor classification. |
| 6b — Variables: diagnostic vs. prognostic clearly distinguished | ✅ PRESENT | Two distinct analysis streams clearly separated (Results subsections "Prognostic model" and "Diagnostic classification"). |
| 7 — Sources of bias addressed | ⚠️ PARTIAL | TOIL batch-effect reduction mentioned. Random split validation is acknowledged. However, no systematic discussion of potential selection bias (TCGA sample selection), detection bias, or information bias. |
| 8 — Sample size justification | ❌ MISSING | No sample size calculation or power analysis. The sample (458 patients, 443 with complete data) is described as "among the larger studies of its kind" but no formal power justification is provided. |
| 9 — Quantitative variables: how handled | ✅ PRESENT | Pathological stage treated as ordinal with dummy coding, stage I reference. Expression values Box--Cox transformed. Log2(FPKM+1) transformation described. |
| 10a — Statistical methods: all described | ✅ PRESENT | limma for DE, univariate Cox, elastic-net Cox, stepwise Cox, RSF, and 5 ML classifiers all described. R 4.6.0, caret, glmnet, survivalROC, pec, randomForestSRC specified. |
| 10b — Subgroup and interaction analyses | ❌ MISSING | No subgroup analyses (e.g., by stage, sex, age) reported or discussed. |
| 10c — Missing data handling | ⚠️ PARTIAL | Missing clinical covariates mentioned: "a substantial proportion contained missing values." But no explicit method for handling missing data is described (e.g., imputation, complete-case analysis). |
| 11 — Sensitivity analyses | ❌ MISSING | No sensitivity analyses reported (e.g., different thresholds, different splits, leave-one-out). "Cox model with stage alone" serves as a reference but is not a formal sensitivity analysis. |
| 12a — Sample size for outcome events | ✅ PRESENT | 94 deaths among 458 patients (event rate ~20.5%) described. 443 patients with complete data. |
| 12b — Description of split-sample validation | ✅ PRESENT | 70/30 train/test split with balanced status described. |

### Results (Items 13-17)
| Item | Status | Notes |
|------|--------|-------|
| 13a — Participants: numbers at each stage | ⚠️ PARTIAL | Sample numbers are given at different pipeline stages (58,387 genes, 27,560 retained, 15,900 lncRNAs, 799 with survival data, 443 patients), but no formal flow diagram and no explicit count of how many were excluded at each filtering step. |
| 13b — Reasons for non-participation | ✅ N/A | Secondary database analysis; non-participation not applicable. |
| 13c — Flow diagram | ❌ MISSING | No participant flow diagram is included. A pipeline overview (Fig.1) shows the computational steps but not patient/sample attrition. |
| 14a — Descriptive data: demographics | ⚠️ PARTIAL | Table S1 is referenced for baseline characteristics of training/testing sets. But key demographics (age, sex distribution) are not discussed in the main text explicitly. |
| 14b — Number of outcome events | ✅ PRESENT | 94 deaths (20.5%) reported. |
| 14c — Follow-up time | ❌ MISSING | Median follow-up time is not reported. |
| 15 — Outcome data | ✅ PRESENT | Survival analysis results, C-index, AUC, accuracy metrics all reported with confidence intervals. |
| 16a — Main results: unadjusted estimates | ✅ PRESENT | Univariate Cox results in Table 1. |
| 16b — Adjusted estimates and confounders | ⚠️ PARTIAL | Multivariate Cox with stage adjustment. But many confounders excluded due to missing data — stated but not addressed. |
| 16c — Continuous variables categorization rationale | ⚠️ PARTIAL | RRS dichotomized at median — rationale is common practice but not explicitly justified. |
| 17 — Other analyses (subgroups, interactions) | ❌ MISSING | No subgroup analyses, interaction tests, or sensitivity analyses reported. |

### Discussion (Items 18-22)
| Item | Status | Notes |
|------|--------|-------|
| 18 — Key results summarized with reference to objectives | ✅ PRESENT | First paragraph of Discussion recaps all main findings with numbers. |
| 19 — Limitations addressed | ✅ PRESENT | Five limitations acknowledged: (1) no external validation, (2) missing clinical covariates, (3) modest sample size, (4) uncharacterized lncRNAs, (5) black-box model interpretability. |
| 20 — Cautious interpretation | ✅ PRESENT | Language is measured ("consistent with," "suggest," "potential"). Comparisons to literature are fair. |
| 21 — Generalizability | ⚠️ PARTIAL | Acknowledges single-institution (TCGA) limitation. No discussion of generalizability to different populations (non-European ancestry, different healthcare settings). |
| 22 — Funding source | ✅ PRESENT | Funding: NSFC (81272493, 81472213), Health Commission of Zhejiang Province (2019331258, 2019335600), Natural Sciences Foundation of Zhejiang (LY17H220001). |

---

## Summary

| Category | Present | Partial | Missing | Score |
|----------|---------|---------|---------|-------|
| Title/Abstract | 1 | 1 | 0 | 50% |
| Introduction | 2 | 0 | 0 | 100% |
| Methods | 4 | 4 | 4 | 33% |
| Results | 3 | 3 | 3 | 33% |
| Discussion | 3 | 1 | 0 | 75% |
| Other | 1 | 0 | 0 | 100% |
| **Overall** | **14** | **9** | **7** | **47%** |

### Critical Missing Items (Pre-Submission Priority)
1. **Sample size justification (Item 8)** — Add a statement explaining that the sample was determined by TCGA availability, with event count reference.
2. **Missing data handling (Item 10c)** — Explicitly state how missing covariates were handled (e.g., complete-case analysis, and state the assumption).
3. **Participant flow diagram (Item 13c)** — Add a diagram showing stepwise exclusion of samples/patients.
4. **Follow-up time (Item 14c)** — Report median follow-up with IQR.
5. **Sensitivity/subgroup analyses (Items 11, 17)** — At minimum, add a sensitivity analysis using different random splits.
