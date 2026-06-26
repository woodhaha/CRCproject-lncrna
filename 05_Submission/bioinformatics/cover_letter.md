# Cover Letter

**Date:** June 25, 2026

**To:** The Editor-in-Chief
*Bioinformatics*
Oxford University Press

**Re:** Submission of Original Paper: "Machine Learning-Based Identification of Long Non-Coding RNA Signatures for Dual Diagnostic and Prognostic Assessment in Colorectal Cancer"

Dear Editor-in-Chief,

We are pleased to submit our manuscript entitled "Machine Learning-Based Identification of Long Non-Coding RNA Signatures for Dual Diagnostic and Prognostic Assessment in Colorectal Cancer" for consideration as an Original Paper in *Bioinformatics*.

## Why Bioinformatics

The core contribution of this work is a fully reproducible, parallelized computational pipeline (v20) for transcriptomic biomarker discovery that integrates differential expression analysis, elastic-net regularized Cox regression, stepwise multivariate Cox modeling, random survival forests, and five machine learning classifiers (SVM, random forest, neural network, logistic regression, elastic-net logistic regression) into a single executable R script (>1,250 lines) with fixed random seeds. This pipeline addresses a methodological gap in the field by simultaneously assessing both diagnostic and prognostic utility of lncRNA signatures within a standardized framework. The audience of *Bioinformatics* -- computational biologists, bioinformaticians, and method developers -- is the ideal readership for a reproducible pipeline architecture that can be readily adapted to other cancer types and transcriptomic biomarker discovery tasks.

## Key Findings

1. **Prognostic model**: A 22-feature model (14 lncRNAs + pathological stage) identified via elastic-net Cox regression followed by stepwise multivariate Cox modeling stratified patients into distinct risk groups (log-rank p < 0.0001), with time-dependent AUC converging to 0.725 at 8 years and a C-index of 0.82.

2. **Diagnostic classifiers**: Five classifiers were systematically compared. SVM, random forest, and neural network achieved testing-set accuracies of 91.6%, 92.2%, and 90.6%, respectively, substantially outperforming logistic regression variants (72.8--73.3%).

3. **Independent validation**: Several identified lncRNAs (EIF3J-AS1, MIR210HG, LINC00261) have independently validated roles in cancer biology, while others represent potentially novel CRC biomarkers.

4. **Reproducibility**: The complete pipeline with fixed random seeds is publicly available at https://github.com/woodhaha/CRC_data_mining, and a stable archival version is deposited on Zenodo (DOI to be assigned upon acceptance).

## Novelty and Significance

This work offers three methodological contributions: (i) a dual diagnostic-prognostic pipeline architecture that avoids the common pitfall of analyzing these two objectives in isolation; (ii) a parallelized, fully reproducible R implementation with explicit seed control; and (iii) a systematic five-classifier comparison on the same harmonized data to benchmark performance. The 37 lncRNA elastic-net panel and the refined 22-feature Cox model provide a concrete resource for the CRC biomarker community.

## Declarations

**AI/LLM Disclosure**: Claude Code (Anthropic, model: Claude Opus 4) was used for R code optimization (parallel loop refactoring, seed management, pipeline orchestration) and manuscript editing (language polishing, formatting). All scientific decisions, study design, data analysis, result interpretation, and final conclusions were made by the authors. The authors bear full responsibility for the integrity of the work.

**Related Manuscripts**: No related manuscripts by the authors are currently under consideration elsewhere.

**Editorial Board**: None of the authors serves on the editorial board of *Bioinformatics*.

**Competing Interests**: The authors declare no competing interests.

We believe the reproducible pipeline architecture and systematic benchmarking presented here will be of strong interest to the *Bioinformatics* readership. We appreciate your consideration and look forward to hearing from you.

Sincerely,

Shaolian Han, M.D., Ph.D.
Department of Gastroenterology Surgery
The First Affiliated Hospital of Wenzhou Medical University
Nanbaixiang Street, Ouhai District, 325000, Wenzhou, Zhejiang, China
Email: wangguanyu@zju.edu.cn

**Corresponding author**: Shaolian Han
