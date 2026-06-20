# CRC lncRNA Analysis Pipeline Flowchart

```mermaid
flowchart TB
    subgraph DATA["📦 Data Sources"]
        A1["TCGA-COAD<br/>HTSeq-FPKM<br/>479 tumor + 41 normal"]
        A2["TCGA + GTEx<br/>TOIL-recomputed<br/>RSEM FPKM (log₂)"]
        A3["GENCODE v22/v23<br/>lncRNA annotations<br/>15,931 lncRNAs"]
        A4["TCGA-COAD Clinical<br/>follow-up + survival<br/>441 patients"]
    end

    subgraph PREP["🔧 Data Preparation"]
        B1["Filter low-expression<br/>genes (>60% samples<br/>with FPKM < 0.1)"]
        B2["Extract colon samples<br/>TCGA=331 · GTEX=308<br/>Total=639"]
        B3["Annotate lncRNAs<br/>v22: 15,900 · v23: 15,931<br/>ENSG → Symbol mapping"]
        B4["Clean clinical data<br/>factor encoding<br/>OS/DFS calculation"]
    end

    subgraph DEA["🔬 Part 1: Differential Expression"]
        C1["limma pipeline<br/>Tumor vs NonTumor<br/>eBayes (trend+robust)"]
        C2["BH-adjusted p-values<br/>FDR < 0.05"]
        C3["1,289 lncRNAs<br/>differentially expressed"]
    end

    subgraph SURV["📊 Part 2: Survival Analysis"]
        D1["Merge clinical + lncRNA<br/>428 patients · 799 lncRNAs"]
        D2["Univariate Cox PH<br/>HR · 95% CI · Wald test"]
        D3["Log-rank test<br/>median-split high/low"]
        D4["Merge DEA + Cox + Log-rank<br/>459 lncRNAs with complete data"]
    end

    subgraph FEAT["🎯 Part 3: Feature Selection"]
        E1["Elastic Net Cox<br/>α-grid search 0→1<br/>10-fold CV"]
        E2["Optimal: α=0.200<br/>λ=0.230 · dev=2.588"]
        E3["17 lncRNAs selected<br/>non-zero coefficients"]
        E4["Train/Test split 70:30<br/>stratified by status"]
        E5["Stepwise Cox PH<br/>AIC minimization"]
        E6["11 features final<br/>8 lncRNAs + stage"]
        E7["Risk score model<br/>Time-dependent ROC<br/>1/3/5-year AUC"]
    end

    subgraph RF["🌲 Part 4: Random Forest Survival"]
        F1["Hyperparameter tuning<br/>mtry · nodesize<br/>OOB error minimization"]
        F2["Optimal: mtry=1<br/>nodesize=2 · ntree=500"]
        F3["RF survival prediction<br/>train/test evaluation"]
        F4["Variable importance<br/>permutation-based VIMP"]
        F5["Brier Score<br/>riskRegression::Score<br/>bootcv B=100"]
        F6["C-index<br/>pec::cindex<br/>bootcv B=100"]
    end

    subgraph ML["🤖 Part 5: Classification"]
        G1["8 lncRNA features<br/>TCGA+GTEx expression<br/>639 samples"]
        G2["Box-Cox normalization<br/>caret::preProcess"]
        G3["Train/Test 70:30<br/>10-fold × 5-repeat CV"]
    end

    subgraph MODELS["🧠 ML Models"]
        H1["SVM<br/>Radial Basis Kernel<br/>σ/C grid search"]
        H2["Random Forest<br/>ranger · gini split<br/>mtry tuning"]
        H3["Neural Network<br/>nnet · size/decay<br/>grid search"]
        H4["Elastic Net<br/>glmnet · α/λ<br/>grid search"]
        H5["Logistic Reg<br/>glmStepAIC<br/>stepwise"]
    end

    subgraph EVAL["📈 Part 6: Model Evaluation"]
        I1["ROC curves<br/>AUC comparison"]
        I2["Lift curves<br/>cumulative gain"]
        I3["Confusion matrix<br/>Acc · Prec · Rec · F1"]
        I4["Model comparison<br/>resamples · bwplot"]
    end

    subgraph OUT["📁 Outputs"]
        J1["36 PDF figures<br/>Nature color scheme"]
        J2["2 interactive HTML<br/>3D surfaces (plotly)"]
        J3["9 data exports<br/>CSV · TXT"]
        J4["RESULTS.md<br/>comprehensive report"]
    end

    A1 --> B1
    A2 --> B2
    A3 --> B3
    A4 --> B4
    B1 --> C1
    B2 --> C1
    B3 --> C1
    C1 --> C2 --> C3

    B4 --> D1
    B1 --> D1
    B3 --> D1
    D1 --> D2 & D3
    D2 & D3 --> D4

    C3 --> D4
    D4 --> E1
    E1 --> E2 --> E3
    E3 --> E4
    E4 --> E5 --> E6
    E6 --> E7

    E6 --> F1
    F1 --> F2 --> F3
    F3 --> F4 & F5 & F6

    E3 --> G1
    G1 --> G2 --> G3
    G3 --> H1 & H2 & H3 & H4 & H5
    H1 & H2 & H3 & H4 & H5 --> I1 & I2 & I3 & I4

    E7 & F3 & F5 & F6 & I1 & I2 & I3 & I4 --> J1 & J2 & J3 & J4

    style DATA fill:#1F4E79,color:#fff,stroke:#1F4E79
    style PREP fill:#2C3E50,color:#fff,stroke:#2C3E50
    style DEA fill:#2471A3,color:#fff,stroke:#2471A3
    style SURV fill:#1A5276,color:#fff,stroke:#1A5276
    style FEAT fill:#C0392B,color:#fff,stroke:#C0392B
    style RF fill:#27AE60,color:#fff,stroke:#27AE60
    style ML fill:#8E44AD,color:#fff,stroke:#8E44AD
    style MODELS fill:#7D3C98,color:#fff,stroke:#7D3C98
    style EVAL fill:#E67E22,color:#fff,stroke:#E67E22
    style OUT fill:#1F4E79,color:#fff,stroke:#1F4E79
```

## Pipeline Summary

| # | Part | Input | Method | Output |
|---|------|-------|--------|--------|
| 1 | **DEA** | 60,498 genes × 639 samples | limma + eBayes | 1,289 DE lncRNAs |
| 2 | **Survival** | 428 patients × 799 lncRNAs | Cox PH + Log-rank | 459 lncRNAs with HR + p |
| 3 | **Feature Selection** | 459 lncRNAs | Elastic Net + Stepwise Cox | 11 features (8 lncRNAs + stage) |
| 4 | **RF Survival** | Train=300, Test=128 | randomForestSRC | OOB error, VIMP, Brier, C-index |
| 5 | **Classification** | 639 samples × 8 lncRNAs | SVM · RF · NNET · EN · LR | 5 trained models |
| 6 | **Evaluation** | 5 models | ROC · Lift · Confusion | Performance comparison |
| 7 | **Export** | All above | Nature-style visualization | 36 PDFs + 2 HTML + 9 data files |

## Key Results

```
                    Train     Test      AUC
SVM                 97.1%  →  92.1%    0.969  ★
Neural Network      96.9%  →  92.1%    0.973  ★
Random Forest      100.0%  →  91.1%    0.974
─────────────────────────────────────────────
Cox Risk Score      C=0.82  →  AUC~0.725 (8yr)
Elastic Net         α=0.20  →  17 lncRNAs selected
```
