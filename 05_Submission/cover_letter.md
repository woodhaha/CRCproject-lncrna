# Cover Letter

**Date:** June 25, 2026

**To:** The Editors
*Frontiers in Oncology*

**Re:** Submission of manuscript: "Machine Learning-Based Identification of Long Non-Coding RNA Signatures for Dual Diagnostic and Prognostic Assessment in Colorectal Cancer"

Dear Editors,

We are pleased to submit our manuscript entitled "Machine Learning-Based Identification of Long Non-Coding RNA Signatures for Dual Diagnostic and Prognostic Assessment in Colorectal Cancer" for consideration for publication in *Frontiers in Oncology*.

## Summary

Colorectal cancer (CRC) remains a leading cause of cancer-related mortality worldwide, yet clinically actionable molecular biomarkers for early detection and prognosis are still limited. In this study, we developed an integrated computational pipeline (v20) to identify lncRNA signatures with both diagnostic and prognostic value in CRC. Using harmonized TCGA--GTEx expression data (639 samples) and TCGA-COAD clinical follow-up data (443 patients), we applied elastic-net regularized Cox regression for feature selection followed by stepwise multivariate Cox modeling and five machine learning classifiers for diagnostic prediction.

## Key Findings

1. **Prognostic model**: An 8-lncRNA panel plus pathological stage stratified patients into distinct risk groups (log-rank p < 0.0001), with time-dependent AUC converging to 0.725 at 8 years and a C-index of 0.82.

2. **Diagnostic classifiers**: Support vector machine (SVM), random forest, and neural network achieved testing-set accuracies of 91.6%, 92.2%, and 90.6%, respectively, substantially outperforming logistic regression variants (72.8--73.3%).

3. **Independent validation**: Several identified lncRNAs (EIF3J-AS1, MIR210HG, LINC00261) have independently validated roles in cancer biology, while others represent potentially novel CRC biomarkers.

4. **Reproducibility**: The parallelized pipeline (v20) with fixed random seeds and a single executable R script (1,250+ lines) is publicly available at https://github.com/woodhaha/CRC_data_mining.

## Significance

This work addresses a critical gap by simultaneously assessing both diagnostic and prognostic utility of lncRNA signatures using a standardized, reproducible machine learning framework. Our findings demonstrate that machine learning-driven analysis of lncRNA expression profiles can yield clinically relevant biomarkers for CRC, and the fully reproducible pipeline serves as a resource for future biomarker discovery studies.

## Declarations

- This manuscript has not been previously published and is not under consideration elsewhere.
- All analyses use publicly available datasets (TCGA, GTEx) under appropriate guidelines.
- The authors declare no competing interests.
- All authors have approved the manuscript and agree with its submission.

We believe this work will interest the readership of *Frontiers in Oncology* as a systematic, reproducible framework for transcriptomic biomarker discovery in colorectal cancer. We appreciate your consideration and look forward to hearing from you.

Sincerely,

Shaolian Han, M.D., Ph.D.
Department of Gastroenterology Surgery
The First Affiliated Hospital of Wenzhou Medical University
Nanbaixiang Street, Ouhai District, 325000, Wenzhou, Zhejiang, China
Email: wangguanyu@zju.edu.cn
