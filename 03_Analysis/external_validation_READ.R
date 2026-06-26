###############################################################################
# external_validation_READ.R v3 — Clean, robust approach
###############################################################################

library(SummarizedExperiment)
library(dplyr); library(survival); library(survminer)
library(ggplot2); library(ggpubr); library(pROC); library(caret)
library(stringr); library(limma)

setwd("D:/Researching/CRCproject/03_Analysis")
OUT_DIR <- "outputs/external_validation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

nb <- "#1F4E79"; nr <- "#C0392B"; ng <- "#27AE60"; ngr <- "#95A5A6"; no <- "#E67E22"

cat("═══════ TCGA-READ External Validation v3 ═══════\n\n")

# ── 1. READ data ──
se <- readRDS("external_validation/READ/READ_SE.rds")
cat(sprintf("READ SE: %d genes × %d samples\n", nrow(se), ncol(se)))

expr <- assay(se, "fpkm_unstrand")
expr <- log2(expr + 0.001)
row_info <- as.data.frame(rowData(se))
col_info <- as.data.frame(colData(se))

# Genes: use gene_name for symbol matching
rownames(expr) <- row_info$gene_name

# Samples: tumor vs normal
is_tumor <- grepl("Tumor", col_info$sample_type, ignore.case = TRUE)
is_normal <- grepl("Normal", col_info$sample_type, ignore.case = TRUE)
expr <- expr[, is_tumor | is_normal]
group <- factor(ifelse(is_tumor[is_tumor | is_normal], "Tumor", "Normal"),
                levels = c("Normal", "Tumor"))
cat(sprintf("Samples: %d Tumor, %d Normal\n", sum(is_tumor), sum(is_normal)))

# ── 2. COAD model ──
ml_data <- read.csv("outputs/Machine_learning_data.csv")
cls_genes <- setdiff(colnames(ml_data), "group")
cat(sprintf("COAD lncRNAs: %d\n", length(cls_genes)))

# Match by gene symbol
common <- intersect(cls_genes, rownames(expr))
cat(sprintf("Mapped to READ: %d / %d\n", length(common), length(cls_genes)))

if (length(common) < 3) {
  cat("⚠ Too few genes. List of COAD genes:\n")
  print(cls_genes)
  cat("READ gene examples:\n")
  print(head(rownames(expr), 20))
  quit(status = 1)
}

cat("COAD genes found:", paste(common, collapse=", "), "\n")

# ── 3. Diagnostic validation ──
cat("\n─── Diagnostic Validation ───\n")

# Build datasets
read_df <- as.data.frame(t(expr[common, , drop = FALSE]))
read_df$group <- group

coad_df <- ml_data[, common, drop = FALSE]
coad_df$group <- factor(ifelse(ml_data$group == "NonTumor", "Normal", "Tumor"),
                        levels = c("Normal", "Tumor"))

cat(sprintf("READ: %d samples, COAD: %d samples, Features: %d\n",
            nrow(read_df), nrow(coad_df), length(common)))

# Scale independently
pp <- preProcess(rbind(coad_df[, common], read_df[, common]), method = c("center", "scale"))
coad_s <- predict(pp, coad_df[, common])
read_s <- predict(pp, read_df[, common])

# Train logistic on COAD
set.seed(12345)
fit <- glm(coad_df$group ~ ., data = coad_s, family = "binomial")

# Predict READ
pred <- predict(fit, newdata = read_s, type = "response")
roc_obj <- roc(read_df$group, pred, quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

cat(sprintf("*** READ Diagnostic AUC: %.3f ***\n", auc_val))

pred_class <- factor(ifelse(pred > 0.5, "Tumor", "Normal"), levels = c("Normal", "Tumor"))
cm <- confusionMatrix(pred_class, read_df$group, positive = "Tumor")
cat(sprintf("Acc=%.3f Sens=%.3f Spec=%.3f\n",
            cm$overall["Accuracy"], cm$byClass["Sensitivity"], cm$byClass["Specificity"]))

cairo_pdf(file.path(OUT_DIR, "Fig_EV1_READ_ROC.pdf"), width = 6, height = 6)
plot.roc(roc_obj, col = nr, lwd = 2.5,
         main = sprintf("External Validation: TCGA-READ\n%d lncRNAs, AUC=%.3f", length(common), auc_val))
abline(0, 1, lty = 2, col = ngr)
legend("bottomright",
       legend = c(sprintf("READ (n=%d, AUC=%.3f)", nrow(read_df), auc_val),
                  sprintf("COAD→READ Logistic Transfer")),
       col = c(nr, NA), lwd = c(2.5, NA), bty = "n")
dev.off()
cat("✓ Fig EV1: READ ROC\n")

# ── 4. Prognostic validation ──
cat("\n─── Prognostic Validation ───\n")

surv_info <- col_info[is_tumor | is_normal, ]
surv_info$OS.time <- ifelse(!is.na(surv_info$days_to_death),
                             surv_info$days_to_death,
                             surv_info$days_to_last_follow_up)
surv_info$OS <- ifelse(surv_info$vital_status == "Dead", 1, 0)
cat(sprintf("Deaths: %d / %d\n", sum(surv_info$OS == 1), nrow(surv_info)))

# Load COAD Cox
cox <- read.csv("outputs/coxph_os.csv", stringsAsFactors = FALSE)
# Get lncRNAs that are in READ AND have significant Cox
cox_sig <- cox %>% filter(HR.pvalue < 0.05) %>% arrange(HR.pvalue)

# Match Cox lncRNAs to READ by symbol
cox_read <- cox_sig[cox_sig$Symbol %in% rownames(expr), ]
# Also match by ENSG (strip version)
read_ensg <- row_info$gene_id[match(rownames(expr), row_info$gene_name)]
names(read_ensg) <- rownames(expr)

ensg_match <- match(sub("\\.[0-9]+$", "", cox_sig$ENSG_id),
                    sub("\\.[0-9]+$", "", read_ensg))
ensg_available <- !is.na(ensg_match)
if (sum(ensg_available) > nrow(cox_read)) {
  extra <- cox_sig[ensg_available & !cox_sig$Symbol %in% cox_read$Symbol, ]
  cox_read <- rbind(cox_read, extra)
}

cox_read <- cox_read %>% distinct(Symbol, .keep_all = TRUE)
cat(sprintf("Cox lncRNAs in READ: %d\n", nrow(cox_read)))

# Filter to symbols actually in expr
cox_read <- cox_read[cox_read$Symbol %in% rownames(expr), ]
cat(sprintf("Cox lncRNAs after symbol filter: %d\n", nrow(cox_read)))

if (nrow(cox_read) >= 3 && sum(surv_info$OS == 1) >= 5) {
  cox_expr <- as.data.frame(t(expr[cox_read$Symbol, , drop = FALSE]))

  # Risk score
  risk <- rep(0, nrow(cox_expr))
  for (i in seq_len(nrow(cox_read))) {
    risk <- risk + cox_read$beta[i] * cox_expr[[i]]
  }

  surv_df <- data.frame(
    OS.time = surv_info$OS.time,
    OS = surv_info$OS,
    risk = risk,
    stage = surv_info$ajcc_pathologic_stage,
    stringsAsFactors = FALSE)
  surv_df <- surv_df[surv_df$OS.time > 0 & !is.na(surv_df$OS.time), ]
  surv_df$risk_group <- factor(ifelse(surv_df$risk > median(surv_df$risk, na.rm = TRUE),
                                      "High", "Low"), levels = c("Low", "High"))
  cat(sprintf("Survival set: %d patients, %d deaths\n", nrow(surv_df), sum(surv_df$OS)))

  # Log-rank
  fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = surv_df)
  lr_p <- surv_pvalue(fit)$pval
  cat(sprintf("Log-rank p = %.4f\n", lr_p))

  # Cox
  cx <- coxph(Surv(OS.time, OS) ~ risk, data = surv_df)
  s <- summary(cx)
  hr <- s$conf.int[1]; hr_lo <- s$conf.int[3]; hr_hi <- s$conf.int[4]
  cox_p <- s$coefficients[5]
  cat(sprintf("Cox HR=%.2f (%.2f-%.2f) p=%.4f\n", hr, hr_lo, hr_hi, cox_p))

  cairo_pdf(file.path(OUT_DIR, "Fig_EV2_READ_KM.pdf"), width = 6, height = 5)
  print(ggsurvplot(fit, data = surv_df, pval = TRUE,
                   palette = c(nb, nr), ggtheme = theme_survminer(),
                   title = "External Validation: TCGA-READ\nCOAD Risk Score — Survival",
                   xlab = "Time (days)", legend = "top",
                   legend.title = "Risk Group", surv.median.line = "hv"))
  dev.off()
  cat("✓ Fig EV2: READ KM\n")
} else {
  cat("⚠ Insufficient genes/events for prognostic validation\n")
}

# Summary
cat("\n═══════ Done ═══════\n")

sink(file.path(OUT_DIR, "validation_summary.txt"))
cat("TCGA-READ External Validation\n")
cat("==============================\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("READ samples:", ncol(expr), "(", sum(is_tumor), "Tumor +", sum(is_normal), "Normal)\n")
cat("Genes mapped:", length(common), "/", length(cls_genes), "\n")
cat("Diagnostic AUC:", round(auc_val, 4), "\n")
if (exists("lr_p")) {
  cat("Prognostic log-rank p:", format.pval(lr_p), "\n")
  cat("Prognostic Cox HR:", round(hr, 2),
      sprintf("(%.2f-%.2f)", hr_lo, hr_hi),
      "p =", format.pval(cox_p), "\n")
}
sink()
cat("✓", normalizePath(OUT_DIR), "\n")
