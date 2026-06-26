# Data Dictionary — CRC lncRNA Biomarker Project

> Auto-generated: 2026-06-25 · Source: TCGA-COAD + GTEx Colon

## Raw Data

| File | Source | Rows | Description |
|------|--------|------|-------------|
| `COAD_HTSeq_FPKM.txt` | GDC Portal | 60,483 × 520 | TCGA-COAD HTSeq FPKM gene expression |
| `TcgaTargetGTEX_phenotype.txt` | UCSC Xena / TOIL | 19,130 × 7 | TCGA + GTEx sample phenotype metadata |
| `TCGA-COAD.followup-fullData.csv` | GDC Portal | clinical | Clinical follow-up (OS, DFS, stage, etc.) |
| `TCGA-COAD_Merge.txt` | GDC Portal | merged | Merged expression + clinical |
| `TCGA_Exp_convertensg2symbol.txt` | Converted | large | Ensemble ID → gene symbol mapping |
| `gencode.v22.long_noncoding_RNAs.gtf` | GENCODE v22 | — | lncRNA annotation (GRCh38) |
| `gencode.v23.long_noncoding_RNAs.gtf` | GENCODE v23 | — | lncRNA annotation (GRCh38) |
| `gencode_v22_RNAs_annotation.txt` | GENCODE v22 | 3-col | extracted: ENSG_id, Type, Symbol |
| `gencode_v23_RNAs_annotation.txt` | GENCODE v23 | 3-col | extracted: ENSG_id, Type, Symbol |
| `RNA Annotation Extraction.txt` | Manual | — | Annotation extraction notes |

## Cleaned Data

| File | Rows × Cols | Description |
|------|-------------|-------------|
| `Clinic_lncRNA_Exprs.txt` | 428 × ~805 | Merged clinical + lncRNA expression (log2 FPKM+0.001) |
| `DEA_OS_analysis.csv` | ~800 × 10 | Differential expression + overall survival merged |
| `Machine_learning_data.csv` | 639 × 9 | 7 diagnostic lncRNAs + group (Tumor/NonTumor) |
| `coxph_os.csv` | ~800 × 5 | Univariate Cox PH results: beta, HR, CI, wald.test, pvalue |

## Key Variables

| Variable | Type | Encoding | Source |
|----------|------|----------|--------|
| `id` | chr | TCGA barcode (12-char) | clinical |
| `OS` | num | Overall survival (days) | clinical |
| `DFS` | num | Disease-free survival (days) | clinical |
| `status` | int | 0=censored, 1=dead | clinical |
| `stage` | factor | i, ii, iii, iv | clinical |
| `age` | num | Years | clinical |
| `gender` | factor | male, female | clinical |
| `group` | factor | Tumor / NonTumor / Training / Testing | derived |

## Diagnostic lncRNAs (7 selected by stepwise Cox)

| ENSG ID | Gene Symbol |
|---------|------------|
| ENSG00000259065.1 | RP5-1021I20.1 |
| ENSG00000272913.1 | — |
| ENSG00000247095.2 | — |
| ENSG00000233223.2 | — |
| ENSG00000259974.2 | — |
| ENSG00000267317.2 | — |
| ENSG00000224272.2 | — |
