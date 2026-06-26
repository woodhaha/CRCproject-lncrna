###############################################################################
# geo_external_validation.R — GSE39582 + GSE17536 external validation
# Uses GEOquery to download, map probes → symbols, validate COAD models
###############################################################################

library(GEOquery)
library(dplyr); library(survival); library(survminer)
library(ggplot2); library(ggpubr); library(pROC); library(caret)
library(stringr); library(limma)

setwd("D:/Researching/CRCproject/03_Analysis")
OUT_DIR <- "outputs/external_validation/GEO"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

nb <- "#1F4E79"; nr <- "#C0392B"; ng <- "#27AE60"; ngr <- "#95A5A6"; no <- "#E67E22"

cat("═══════ GEO External Validation ═══════\n")
cat("Datasets: GSE39582 (n=566), GSE17536 (n=177)\n\n")

# ── 1. Download & load GEO datasets ──
cat("─── Downloading from GEO ───\n")

# GSE39582: 566 colon cancer, Affy Human Genome U133 Plus 2.0
cat("  GSE39582...")
gse39582 <- getGEO("GSE39582", GSEMatrix = TRUE, getGPL = TRUE)
gse39582 <- gse39582[[1]]
cat(sprintf(" %d samples\n", ncol(gse39582)))

# GSE17536: 177 CRC, same platform
cat("  GSE17536...")
gse17536 <- getGEO("GSE17536", GSEMatrix = TRUE, getGPL = TRUE)
gse17536 <- gse17536[[1]]
cat(sprintf(" %d samples\n", ncol(gse17536)))

# ── 2. Process expression ──
cat("\n─── Processing expression ───\n")

process_geo <- function(gse) {
  # Extract expression
  expr <- exprs(gse)
  if (is.null(expr)) stop("No expression data")

  # Get platform annotation for probe → gene mapping
  fd <- fData(gse)
  pd <- pData(gse)

  # Try common probe annotation columns
  symbol_col <- grep("Symbol|symbol|GENE_SYMBOL|Gene Symbol", colnames(fd), value = TRUE)[1]
  if (is.na(symbol_col)) {
    # Fallback: use "Gene" or "gene" column
    symbol_col <- grep("^Gene$|^gene$|gene_assignment", colnames(fd), value = TRUE)[1]
  }
  cat(sprintf("  Symbol column: %s\n", symbol_col))

  # Get probe → symbol mapping
  if (!is.na(symbol_col)) {
    symbols <- fd[[symbol_col]]
    # Handle multi-gene probes: take first
    symbols <- str_split(symbols, " /// |;|,", simplify = TRUE)[, 1]
    symbols <- trimws(symbols)
  } else {
    symbols <- fd$ID  # fallback to probe ID
  }

  # Aggregate to gene-level: median of probes per gene
  probe_to_gene <- data.frame(probe = fd$ID, gene = symbols, stringsAsFactors = FALSE)
  probe_to_gene <- probe_to_gene[probe_to_gene$gene != "" & !is.na(probe_to_gene$gene), ]

  # Keep only probes mapping to a single gene
  gene_counts <- table(probe_to_gene$gene)

  # For each gene, take the probe with highest mean expression
  gene_list <- unique(probe_to_gene$gene)
  gene_expr <- matrix(NA, nrow = length(gene_list), ncol = ncol(expr))
  rownames(gene_expr) <- gene_list
  colnames(gene_expr) <- colnames(expr)

  for (i in seq_along(gene_list)) {
    g <- gene_list[i]
    probes <- probe_to_gene$probe[probe_to_gene$gene == g]
    probes <- intersect(probes, rownames(expr))
    if (length(probes) == 1) {
      gene_expr[i, ] <- expr[probes, ]
    } else if (length(probes) > 1) {
      gene_expr[i, ] <- apply(expr[probes, , drop = FALSE], 2, median, na.rm = TRUE)
    }
  }
  gene_expr <- gene_expr[rowSums(is.na(gene_expr)) < ncol(gene_expr), ]

  list(expr = gene_expr, pd = pd)
}

res_39582 <- process_geo(gse39582)
res_17536 <- process_geo(gse17536)

cat(sprintf("  GSE39582: %d genes\n", nrow(res_39582$expr)))
cat(sprintf("  GSE17536: %d genes\n", nrow(res_17536$expr)))

# ── 3. Extract survival data ──
cat("\n─── Extracting survival ───\n")

extract_survival <- function(pd, name) {
  cat(sprintf("  %s pheno columns (first 30):\n", name))
  print(head(colnames(pd), 30))

  # Look more specifically at the data
  cat(sprintf("  Checking key columns...\n"))

  # GSE39582: "os.delay (months):ch1" and "os.event:ch1"
  # GSE17536: "overall survival follow-up time:ch1" and "overall_event (death from any cause):ch1"

  # Try to find OS time
  os_time_col <- grep("^os\\.delay|^overall survival follow.up|os_time|^overall.*follow.*time",
                       colnames(pd), ignore.case = TRUE, value = TRUE)[1]
  # OS event
  os_event_col <- grep("^os\\.event|^overall_event|os_event|overall.*event.*death",
                        colnames(pd), ignore.case = TRUE, value = TRUE)[1]

  if (is.na(os_time_col)) {
    # Try to guess: look for "months" in column names
    os_time_col <- grep("months|years|days", colnames(pd), ignore.case = TRUE, value = TRUE)[1]
  }
  if (is.na(os_event_col)) {
    os_event_col <- grep("event|status|recurrence|relapse|death", colnames(pd), ignore.case = TRUE, value = TRUE)[1]
  }

  cat(sprintf("  %s: time='%s' event='%s'\n", name, os_time_col, os_event_col))

  # Print some values
  if (!is.na(os_time_col)) {
    cat(sprintf("  Time values sample: %s\n", paste(head(pd[[os_time_col]], 10), collapse=", ")))
  }
  if (!is.na(os_event_col)) {
    cat(sprintf("  Event values sample: %s\n", paste(head(pd[[os_event_col]], 10), collapse=", ")))
  }

  if (!is.na(os_time_col) && !is.na(os_event_col)) {
    OS.time <- as.numeric(as.character(pd[[os_time_col]]))
    # Handle text-based event coding
    raw_event <- tolower(as.character(pd[[os_event_col]]))
    # "death" / "dead" / "recurrence" → 1; "no death" / "alive" / "no recurrence" → 0
    OS <- ifelse(grepl("^death$|^dead$|^yes$|^1$|^event$|^recurrence$|^relapse$|^metastasis$",
                       raw_event, ignore.case = TRUE), 1, 0)
    # Also handle numeric coding
    if (all(OS == 0) && any(raw_event %in% c("0", "1"))) {
      OS <- as.numeric(raw_event)
    }

    cat(sprintf("    OS.time range: %.0f-%.0f, events: %d/%d\n",
                min(OS.time, na.rm=TRUE), max(OS.time, na.rm=TRUE),
                sum(OS==1, na.rm=TRUE), length(OS)))
    data.frame(OS.time = OS.time, OS = OS, stringsAsFactors = FALSE)
  } else {
    cat(sprintf("  ⚠ No clear survival columns for %s\n", name))
    cat(sprintf("  Available columns: %s\n", paste(head(colnames(pd), 20), collapse=", ")))
    data.frame(OS.time = NA_real_, OS = NA_integer_)
  }
}

surv_39582 <- extract_survival(res_39582$pd, "GSE39582")
surv_17536 <- extract_survival(res_17536$pd, "GSE17536")

# ── 4. Sample classification (tumor vs normal) ──
cat("\n─── Sample classification ───\n")

classify_samples <- function(pd, name) {
  # Try tissue/source columns
  tissue_col <- grep("tissue|source|disease|tumor|histolog|group|characteristics",
                     colnames(pd), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(tissue_col)) {
    vals <- as.character(pd[[tissue_col]])
    is_normal <- grepl("normal|adjacent|healthy|non.tumor", vals, ignore.case = TRUE)
    is_tumor <- grepl("tumor|cancer|carcinoma|adenocarcinoma|malignant", vals, ignore.case = TRUE)
    data.frame(normal = is_normal, tumor = is_tumor, stringsAsFactors = FALSE)
  } else {
    cat(sprintf("  ⚠ No tissue column for %s\n", name))
    data.frame(normal = FALSE, tumor = TRUE)  # assume all tumor
  }
}

type_39582 <- classify_samples(res_39582$pd, "GSE39582")
type_17536 <- classify_samples(res_17536$pd, "GSE17536")

cat(sprintf("  GSE39582: %d tumor, %d normal\n",
            sum(type_39582$tumor), sum(type_39582$normal)))
cat(sprintf("  GSE17536: %d tumor, %d normal\n",
            sum(type_17536$tumor), sum(type_17536$normal)))

# ── 5. Match COAD lncRNAs ──
cat("\n─── Matching COAD lncRNAs ───\n")

ml_data <- read.csv("outputs/Machine_learning_data.csv")
cls_genes <- setdiff(colnames(ml_data), "group")

common_39582 <- intersect(cls_genes, rownames(res_39582$expr))
common_17536 <- intersect(cls_genes, rownames(res_17536$expr))

cat(sprintf("  GSE39582: %d / %d lncRNAs\n", length(common_39582), length(cls_genes)))
cat(sprintf("  GSE17536: %d / %d lncRNAs\n", length(common_17536), length(cls_genes)))

# Print what we have
if (length(common_39582) > 0) cat("  GSE39582 matches:", paste(common_39582, collapse=", "), "\n")
if (length(common_17536) > 0) cat("  GSE17536 matches:", paste(common_17536, collapse=", "), "\n")

# ── 6. Diagnostic validation ──
cat("\n─── Diagnostic Validation ───\n")

run_diagnostic <- function(geo_expr, geo_type, common_genes, coad_ml, dataset_name) {
  if (length(common_genes) < 3) {
    cat(sprintf("  %s: ⚠ Too few genes (%d), skipping\n", dataset_name, length(common_genes)))
    return(NULL)
  }

  # Filter to normals + tumors
  use <- geo_type$normal | geo_type$tumor
  geo_expr <- geo_expr[common_genes, use, drop = FALSE]
  geo_group <- factor(ifelse(geo_type$normal[use], "Normal", "Tumor"),
                      levels = c("Normal", "Tumor"))

  cat(sprintf("  %s: %d samples (%d N, %d T), %d genes\n",
              dataset_name, ncol(geo_expr),
              sum(geo_group == "Normal"), sum(geo_group == "Tumor"),
              length(common_genes)))

  if (sum(geo_group == "Normal") < 3) {
    cat(sprintf("  %s: ⚠ Too few normals, skipping\n", dataset_name))
    return(NULL)
  }

  geo_df <- as.data.frame(t(geo_expr))
  coad_df <- coad_ml[, common_genes, drop = FALSE]
  coad_df$group <- factor(ifelse(coad_ml$group == "NonTumor", "Normal", "Tumor"),
                          levels = c("Normal", "Tumor"))

  # Scale
  pp <- preProcess(rbind(coad_df[, common_genes], geo_df), method = c("center", "scale"))
  coad_s <- predict(pp, coad_df[, common_genes])
  geo_s <- predict(pp, geo_df)

  # Train on COAD, test on GEO
  set.seed(12345)
  fit <- glm(coad_df$group ~ ., data = coad_s, family = "binomial")
  pred <- predict(fit, newdata = geo_s, type = "response")
  roc_obj <- roc(geo_group, pred, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))

  cat(sprintf("  *** %s Diagnostic AUC: %.3f ***\n", dataset_name, auc_val))

  list(auc = auc_val, roc = roc_obj, pred = pred, group = geo_group,
       genes = common_genes, n = ncol(geo_expr))
}

diag_39582 <- run_diagnostic(res_39582$expr, type_39582, common_39582, ml_data, "GSE39582")
diag_17536 <- run_diagnostic(res_17536$expr, type_17536, common_17536, ml_data, "GSE17536")

# ── 7. Prognostic validation ──
cat("\n─── Prognostic Validation ───\n")

cox <- read.csv("outputs/coxph_os.csv", stringsAsFactors = FALSE)
cox_sig <- cox %>% filter(HR.pvalue < 0.05) %>% arrange(HR.pvalue)

run_prognostic <- function(geo_expr, surv_df, cox_data, dataset_name) {
  cox_genes <- intersect(cox_sig$Symbol, rownames(geo_expr))
  cat(sprintf("  %s: %d Cox lncRNAs\n", dataset_name, length(cox_genes)))

  if (length(cox_genes) < 3) {
    cat(sprintf("  %s: ⚠ Too few Cox genes\n", dataset_name))
    return(NULL)
  }

  # Risk score
  risk <- rep(0, ncol(geo_expr))
  for (g in cox_genes) {
    beta <- cox_sig$beta[cox_sig$Symbol == g][1]
    if (!is.na(beta) && beta != 0) {
      risk <- risk + beta * geo_expr[g, ]
    }
  }

  df <- data.frame(
    OS.time = surv_df$OS.time,
    OS = surv_df$OS,
    risk = risk,
    stringsAsFactors = FALSE)
  df <- df[!is.na(df$OS.time) & df$OS.time > 0 & !is.na(df$OS), ]
  df$risk_group <- factor(ifelse(df$risk > median(df$risk, na.rm = TRUE),
                                    "High", "Low"), levels = c("Low", "High"))
  cat(sprintf("  %s: %d patients, %d events\n", dataset_name, nrow(df), sum(df$OS)))

  if (sum(df$OS) < 5) {
    cat(sprintf("  %s: ⚠ Too few events\n", dataset_name))
    return(NULL)
  }

  fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = df)
  # Use survdiff for log-rank p-value (avoids survminer scoping issues)
  lr_p <- 1 - pchisq(survdiff(Surv(OS.time, OS) ~ risk_group, data = df)$chisq, 1)

  cx <- coxph(Surv(OS.time, OS) ~ risk, data = df)
  s <- summary(cx)
  hr <- s$conf.int[1]; hr_lo <- s$conf.int[3]; hr_hi <- s$conf.int[4]
  cox_p <- s$coefficients[5]

  cat(sprintf("  Log-rank p=%.4f  Cox HR=%.2f (%.2f-%.2f) p=%.4f\n",
              lr_p, hr, hr_lo, hr_hi, cox_p))

  list(logrank_p = lr_p, cox_hr = hr, cox_hr_lo = hr_lo,
       cox_hr_hi = hr_hi, cox_p = cox_p, fit = fit, surv = df,
       genes = cox_genes, n = nrow(df), events = sum(df$OS))
}

prog_39582 <- run_prognostic(res_39582$expr, surv_39582, cox_sig, "GSE39582")
prog_17536 <- run_prognostic(res_17536$expr, surv_17536, cox_sig, "GSE17536")

# ── 8. Plots ──
cat("\n─── Generating plots ───\n")

# Combined ROC
cairo_pdf(file.path(OUT_DIR, "Fig_EV3_GEO_ROC.pdf"), width = 8, height = 6)
any_roc <- FALSE
if (!is.null(diag_39582)) {
  plot.roc(diag_39582$roc, col = nb, lwd = 2.5, main = "External Validation — GEO Datasets\nDiagnostic ROC Curves")
  any_roc <- TRUE
} else {
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "FPR", ylab = "TPR", main = "External Validation — GEO Datasets\nDiagnostic ROC Curves")
  abline(0, 1, lty = 2, col = ngr)
}
if (!is.null(diag_17536)) {
  lines.roc(diag_17536$roc, col = nr, lwd = 2.5)
}
if (!is.null(diag_39582) || !is.null(diag_17536)) {
  leg_text <- c()
  leg_col <- c()
  if (!is.null(diag_39582)) {
    leg_text <- c(leg_text, sprintf("GSE39582 (n=%d, AUC=%.3f)", diag_39582$n, diag_39582$auc))
    leg_col <- c(leg_col, nb)
  }
  if (!is.null(diag_17536)) {
    leg_text <- c(leg_text, sprintf("GSE17536 (n=%d, AUC=%.3f)", diag_17536$n, diag_17536$auc))
    leg_col <- c(leg_col, nr)
  }
  legend("bottomright", legend = leg_text, col = leg_col, lwd = 2.5, bty = "n")
}
dev.off()
cat("  ✓ Fig EV3: GEO ROC\n")

# Combined KM
cairo_pdf(file.path(OUT_DIR, "Fig_EV4_GEO_KM.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
for (ds in c("GSE39582", "GSE17536")) {
  obj <- if (ds == "GSE39582") prog_39582 else prog_17536
  if (!is.null(obj)) {
    p <- ggsurvplot(obj$fit, data = obj$surv, pval = TRUE,
                    palette = c(nb, nr), ggtheme = theme_survminer(),
                    title = sprintf("%s (n=%d)", ds, obj$n),
                    xlab = "Time", legend = "top",
                    legend.title = "Risk", surv.median.line = "hv")
    print(p)
  } else {
    plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
         main = paste(ds, "- insufficient data"))
  }
}
dev.off()
cat("  ✓ Fig EV4: GEO KM\n")

# ── 9. Summary ──
cat("\n═══════ Summary ═══════\n")

sink(file.path(OUT_DIR, "geo_validation_summary.txt"))
cat("GEO External Validation Summary\n")
cat("================================\n")
cat("Date:", as.character(Sys.time()), "\n\n")
for (ds in c("GSE39582", "GSE17536")) {
  cat("---", ds, "---\n")
  common <- if (ds == "GSE39582") common_39582 else common_17536
  diag <- if (ds == "GSE39582") diag_39582 else diag_17536
  prog <- if (ds == "GSE39582") prog_39582 else prog_17536
  cat("Genes mapped:", length(common), "/", length(cls_genes), "\n")
  if (!is.null(diag)) cat("Diagnostic AUC:", round(diag$auc, 4), "\n")
  if (!is.null(prog)) {
    cat("Prognostic log-rank p:", format.pval(prog$logrank_p), "\n")
    cat("Prognostic Cox HR:", round(prog$cox_hr, 2),
        sprintf("(%.2f-%.2f)", prog$cox_hr_lo, prog$cox_hr_hi),
        "p =", format.pval(prog$cox_p), "\n")
  }
  cat("\n")
}
cat("COAD lncRNAs:", length(cls_genes), "\n")
sink()

cat("✓ Results:", normalizePath(OUT_DIR), "\n")
cat("═══ GEO Validation Complete ═══\n")
