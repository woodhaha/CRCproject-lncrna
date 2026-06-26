# Self-Review Report — CRC lncRNA Manuscript
**Manuscript**: `manuscript_rewritten.tex` (v20 pipeline, Nature-style)
**Review date**: 2026-06-25

---

## 1. Internal Consistency — Abstract vs. Results

| Claim in Abstract | Source in Results | Verdict |
|---|---|---|
| 37 lncRNAs with non-zero coefficients (elastic-net) | Line 56: "identified 37 lncRNAs with non-zero coefficients" | ✅ MATCH |
| 443 patients with complete survival data | Line 54: "443 had complete survival, stage, and lncRNA expression information" | ✅ MATCH |
| 22-feature final model | Line 56: "final model of 22 features" | ✅ MATCH |
| AIC = 593.46 | Line 56: "AIC of 593.46" | ✅ MATCH |
| SVM testing accuracy 91.6% | Table 2, Line 192: "Testing & 0.916" | ✅ MATCH |
| RF testing accuracy 92.2% | Table 2, Line 195: "Testing & 0.922" | ✅ MATCH |
| NN testing accuracy 90.6% | Table 2, Line 198: "Testing & 0.906" | ✅ MATCH |
| Logistic regression: 72.8% | Table 2, Line 186: "Testing & 0.728" | ✅ MATCH |
| Elastic-net logistic: 73.3% | Table 2, Line 189: "Testing & 0.733" | ✅ MATCH |
| AUC converges to 0.725 at 8 years | Line 140: "converged to approximately 0.725 after 8 years" | ✅ MATCH |
| TCGA--GTEx: 349 non-tumor, 290 tumor | Line 52 & 233 | ✅ MATCH |
| "458 patients, 443 with complete data" in abstract | Line 54: "458 patients, of whom 443 had complete" | ✅ MATCH |
| C-index = 0.82 | Line 56: "C-index of 0.82" | ✅ MATCH |

**PASS**: All key numbers are internally consistent between abstract and results.

---

## 2. Figure References vs. Figure Environments

| \ref{FigX} in text | \begin{figure} present? | \includegraphics file | File on disk? |
|---|---|---|---|
| Fig.~\ref{Fig1} | ✅ Line 60 | Fig1_candidates.pdf | ✅ 24.9K |
| Fig.~\ref{Fig2} | ✅ Line 118 | Fig2a_deviance.pdf | ✅ 46.3K |
| Fig.~\ref{Fig3} | ✅ Line 131 | Fig3_forest.pdf | ✅ 46.8K |
| Fig.~\ref{Fig4} | ✅ Line 144 | Fig4a_km_train.pdf | ✅ 24.3K |
| Fig.~\ref{Fig5} | ✅ Line 151 | Fig5a_roc.pdf | ✅ 34.3K |
| Fig.~\ref{Fig6} | ✅ Line 158 | Fig6a_brier.pdf | ✅ 26.7K |
| Fig.~\ref{Fig7} | ✅ Line 204 | Fig7a_roc_train.pdf | ✅ 37.6K |

### Issues Found

**ISSUE 1 (Medium): Fig4 references c,d but only one image file included**
- Text (line 142-143): References "Fig.~\ref{Fig4}c, d" — RSF Kaplan--Meier panels
- But \includegraphics (line 146) only includes `Fig4a_km_train.pdf`
- **The figure caption also only mentions "training (left) and testing (right)"**, implying a single combined image. If Fig4a_km_train.pdf contains both panels, this is fine. But the text reference to "Fig4c,d" suggests panels that don't exist in the included file.
- **Recommendation**: Either combine all Fig4 panels into one image file, or add separate \includegraphics for Fig4b_km_test.pdf, Fig4c_rsf_train, Fig4d_rsf_test.

**ISSUE 2 (Medium): Fig5 only shows ROC (a) but text references time-dependent AUC**
- Text (line 140-141): References Fig5 for time-dependent AUC converging to 0.725
- \includegraphics shows `Fig5a_roc.pdf` (ROC curves at 3 and 5 years)
- `Fig5b_auc.pdf` exists on disk (23.0K) but is **NOT included** in the manuscript
- **Recommendation**: Add Fig5b to show the time-dependent AUC trajectory over 8 years.

**ISSUE 3 (High): Fig6 only references a,b in image but text mentions c,d**
- \includegraphics (line 160): `Fig6a_brier.pdf`
- Text (line 142): References "Fig.~\ref{Fig6}c, d" for C-index and stage-alone comparison
- `Fig6b_cindex.pdf` exists on disk (21.4K) but is **NOT included** in the \includegraphics
- Text claims panels c,d show C-index with stage-alone comparison, but those panels are NOT in the figure file.
- **Recommendation**: Create a combined Fig6 image with all four panels (a: Brier train, b: Brier test, c: C-index, d: C-index stage-alone), or adjust \includegraphics to include separate panels.

**ISSUE 4 (Medium): Fig7 image only includes (a) ROC train**
- \includegraphics shows `Fig7a_roc_train.pdf`
- Text (line 169-171) references "Fig.~\ref{Fig7}a, b" for ROC curves and "Fig.~\ref{Fig7}c, d" for lift curves
- `Fig7b_roc_test.pdf` exists on disk (37.0K) but is **NOT included**
- No lift-curve PDFs exist in the directory
- **Recommendation**: Include Fig7b_roc_test.pdf and create lift curve images.

**ISSUE 5 (Low): Fig2b unused file on disk**
- `Fig2b_coef_paths.pdf` (33.4K) exists on disk but is never referenced in the manuscript
- The Fig2 caption mentions coefficient shrinkage paths "(b)" so perhaps it was intended to be combined. If it is a combined image, this is fine.

---

## 3. Missing Figure Files

| \includegraphics reference | File exists? |
|---|---|
| `Fig1_candidates` | ✅ `Fig1_candidates.pdf` (24.9K) |
| `Fig2a_deviance` | ✅ `Fig2a_deviance.pdf` (46.3K) |
| `Fig3_forest` | ✅ `Fig3_forest.pdf` (46.8K) |
| `Fig4a_km_train` | ✅ `Fig4a_km_train.pdf` (24.3K) |
| `Fig5a_roc` | ✅ `Fig5a_roc.pdf` (34.3K) |
| `Fig6a_brier` | ✅ `Fig6a_brier.pdf` (26.7K) |
| `Fig7a_roc_train` | ✅ `Fig7a_roc_train.pdf` (37.6K) |

**All \includegraphics files exist on disk.** However, as noted above, supplementary panels (Fig4b, Fig5b, Fig6b, Fig7b, Fig7c-d lift curves) exist on disk but are not included, or don't exist as separate files.

### Supplementary Figure Files

| Reference | File exists? |
|---|---|
| Fig.~S1 (stage survival, risk dist.) | ✅ FigS1_stage_survival.pdf |
| Fig.~S2, S3 (RF variable importance) | ✅ FigS2_rf_vimp.pdf, FigS3_rf_vimp_bar.pdf |
| Fig.~S4, S5 (risk score) | ✅ FigS4_riskscore.pdf, FigS5_risk_boxplot.pdf |
| Fig.~S6 (RF survival) | ✅ FigS6_rf_survival.pdf |
| Table S1 (baseline chars.) | Referenced but not in main .tex file — presumably in supplementary.pdf |
| Table S2 (risk classification) | Referenced but not in main .tex file |

---

## 4. Statistical Reporting

### p-values
| Pattern | Examples | Verdict |
|---|---|---|
| Exact p-values | 0.003, 0.004, 0.007, 0.012, 0.022 | ✅ Consistent formatting |
| `$<0.001$` | Used throughout Table 1 where warranted | ✅ Proper threshold formatting |
| Scientific notation | `$8.97 \times 10^{-9}$` in Table 2 | ✅ Correct |
| Log-rank p-values | 0.008, 0.037, 0.002, etc. in Table 1 | ✅ Present |

### CI Ranges
| Issue | Location | Verdict |
|---|---|---|
| Hyphenation of CI ranges | All CIs in Table 1 use `--` (en-dash): (1.21--2.65) | ✅ Correct — en-dashes for ranges |
| CI in Table 2 accuracy | (0.63--0.72) etc. | ✅ Correct |
| HR CI in text | "0.717" and "0.802" in Discussion (RN17, RN18 citations) | ✅ But these are reported without their own CIs — minor issue |

### Issues Found

**ISSUE 6 (Medium): Several lncRNAs in Table 1 with non-significant p-values are still listed**
- RP11-440D17.3: FDR = 0.152, Cox p = 0.064, log-rank p = 0.200 — NONE reach significance
- CTB-25B13.12: FC = 0.014, FDR = 0.888 — essentially not differentially expressed
- MIR210HG: FC = 0.129, FDR = 0.352 — not significant for DE
- TNRC6C-AS1: FC = -0.140, FDR = 0.213 — not significant for DE
- **These are included because they passed the Cox-based selection, but their FC/FDR values are non-significant.** Table 1 caption should clarify the selection criteria more precisely (the Methods section states stricter criteria of |log2FC| >= 0.5 and FDR < 0.05; these lncRNAs would not have passed those).

**ISSUE 7 (Minor): Inconsistency in "8 lncRNAs" vs "22 features"**
- Methods state 8 lncRNAs for RSF and classification (lines 142, 217, 261, 269)
- But the final Cox model has 22 features including stage (line 56 and 129)
- The relationship between the 22 features (21 lncRNAs + 1 stage) and the 8 lncRNAs used in RSF and classifiers is **never explicitly explained**.
- **Recommendation**: Add a sentence explaining that the 8 lncRNAs are a subset of the 22 features retained for downstream models, or clarify the selection logic.

---

## 5. Discussion Claims vs. Results

| Discussion Claim | Results Evidence | Verdict |
|---|---|---|
| "37 lncRNAs with non-zero coefficients" | Line 56 confirms | ✅ Supported |
| "22-feature stepwise Cox model" | Line 56 confirms | ✅ Supported |
| "Parallelized computation accelerates by ~70%" | Line 213: "approximately 70% compared to serial execution" | ⚠️ CLAIMED — but no data/results section table or supplemental figure showing this benchmark. Should add benchmark details. |
| "RSF underperformed" | Line 142: OOB 25.011% vs. 32.169% | ✅ Supported |
| "Time-dependent AUC stabilized at 0.725" | Line 140 confirms | ✅ Supported |
| Diagnostic classifiers > 90% | Table 2 confirms (SVM 91.6%, RF 92.2%, NN 90.6%) | ✅ Supported |
| "Strong agreement between training and testing performance" | SVM: 97.8% vs. 91.6%, NN: 98.2% vs. 90.6% | ⚠️ PARTIAL — 6-7% gap still exists. This is reasonable but "strong agreement" may be slightly overstated for NN (7.6% gap). |
| "Several lncRNAs independently implicated" | Literature citations provided (RN20-RN24) | ✅ Supported |
| "Pipeline (v20) publicly available" | Line 293: GitHub URL provided | ✅ Supported |

### ISSUE 8 (Minor): 70% speedup claim is unverifiable
- The Discussion claims "approximately 70% compared to serial execution"
- No benchmark data is presented in Results, Methods, or Supplementary
- **Recommendation**: Add a small table or sentence showing the benchmark.

---

## 6. Limitations Assessment

| Original Limitation | Retained in manuscript? | Location |
|---|---|---|
| 1. No external validation (TCGA-only) | ✅ Yes | Line 223 |
| 2. Missing clinical covariates | ✅ Yes | Line 223 |
| 3. Modest sample size, potential optimism | ✅ Yes | Line 223 |
| 4. Uncharacterized lncRNAs need functional validation | ✅ Yes | Line 223 |
| 5. Black-box model interpretability | ✅ Yes | Line 223 |

**All 5 limitations from the original analysis are present and well-articulated in the Discussion.**

---

## 7. Author Information

| Field | Status | Notes |
|---|---|---|
| Author names | ✅ Present | Zhuha Zhou, Yongyu Bai, Qigang Xue, Zhuxian Zhou, Shaolian Han |
| Affiliations | ✅ Present | 3 affiliations listed with full addresses |
| Equal contribution | ✅ Present | "Zhuha Zhou and Yongyu Bai contributed equally" |
| Contributions | ✅ Present | Z.Z. and Y.B.: conceptualization, methodology, software, formal analysis, writing---original draft. S.H.: supervision, writing---review and editing. |
| Competing interests | ✅ Present | "The authors declare no competing interests." |
| Corresponding author | ✅ Present | Correspondence to Shaolian Han |

### ISSUE 9 (Medium-High): Correspondence email mismatch
- Line 24: `wangguanyu@zju.edu.cn` — This email address does NOT match Shaolian Han's name or affiliation.
- Shaolian Han is at Wenzhou Medical University (not Zhejiang University).
- The email `wangguanyu@zju.edu.cn` would belong to a different person (possibly Qigang Xue or Zhuxian Zhou?).
- **This is a critical correspondence error that must be fixed before submission.**

---

## Summary of Issues by Severity

### Critical (Must Fix Before Submission)
1. **ISSUE 9** — Correspondence email `wangguanyu@zju.edu.cn` does not match Shaolian Han or Wenzhou Medical University. It belongs to Zhejiang University, suggesting a copy-paste error.

### High Priority
2. **ISSUE 3** — Fig6 references panels c,d that are not included in the image. Text discusses C-index and stage-alone comparison, but only `Fig6a_brier.pdf` is included.
3. **ISSUE 1** — Fig4 references panels c,d (RSF results) but only `Fig4a_km_train.pdf` is included.

### Medium Priority
4. **ISSUE 2** — Fig5b_auc.pdf exists on disk but is not included in the manuscript. The AUC trajectory over 8 years is a key result.
5. **ISSUE 4** — Fig7 references 4 panels (a-d) but only panel (a) ROC train is included. Missing ROC test and lift curves.
6. **ISSUE 6** — Table 1 includes lncRNAs with non-significant DE values (RP11-440D17.3, CTB-25B13.12, TNRC6C-AS1, MIR210HG). Caption should clarify selection criteria.
7. **ISSUE 7** — "8 lncRNAs" vs "22 features" gap not explained. Should clarify that 8 lncRNAs are a downstream subset.

### Minor
8. **ISSUE 8** — 70% speedup claim has no supporting benchmark data.
9. **ISSUE 5** — Fig2b_coef_paths.pdf on disk but unreferenced (possibly combined with Fig2a).

### STROBE Gaps (from companion report)
10. No formal sample size justification
11. Missing data handling method not explicitly stated
12. No participant flow diagram
13. No median follow-up time reported
14. No sensitivity/subgroup analyses
