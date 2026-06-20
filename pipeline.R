###############################################################################
# CRC lncRNA Biomarker Discovery — Consolidated Pipeline
# Nature-style color scheme · Output to outputs/
# Refactored 2026-06-20
###############################################################################

# ── Setup ───────────────────────────────────────────────────────────────────
DATA_DIR   <- "data"
OUT_DIR    <- "outputs"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Nature-inspired color palette ───────────────────────────────────────────
# Based on Nature Publishing Group / Nature Reviews style
nature_blue   <- "#1F4E79"   # deep navy — primary
nature_red    <- "#C0392B"   # crimson — contrast
nature_green  <- "#27AE60"   # emerald
nature_orange <- "#E67E22"   # burnt orange
nature_purple <- "#8E44AD"   # deep purple
nature_teal   <- "#16A085"   # dark teal
nature_grey   <- "#95A5A6"   # slate grey

# Paired palettes
pal_train_test  <- c(nature_blue, nature_red)        # Training / Testing
pal_risk        <- c(nature_blue, nature_red)        # Low / High risk
pal_tumor       <- c(nature_blue, nature_red)        # NonTumor / Tumor
pal_status      <- c(nature_blue, nature_red)        # Censored / Dead
pal_models      <- c(nature_blue, nature_red, nature_green, nature_orange, nature_purple)
pal_gradient    <- c(nature_red, nature_orange, nature_blue)  # High → Low for heat

# ── Package loader ──────────────────────────────────────────────────────────
LoadRpak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

pkgs <- c("ggplot2", "dplyr", "reshape2", "RColorBrewer", "scales", "grid",
          "caret", "glmnet", "pROC", "ggpubr", "ranger", "ggfortify", "tidyr",
          "survminer", "survival", "stringr", "table1", "xtable",
          "randomForestSRC", "ggRandomForests", "lattice", "c060", "limma",
          "broom", "survivalROC", "pec", "prodlim", "riskRegression", "doParallel",
          "plotly", "htmlwidgets", "data.table", "ggsci")

LoadRpak(pkgs)

# ── Fallback functions (packages that may fail on some systems) ──────────────
# Simple lookup: maps terms to values via a 2-column match table
lookup <- function(terms, key.match, missing = NA) {
  if (is.data.frame(key.match) && ncol(key.match) >= 2) {
    idx <- match(terms, key.match[[1]])
    out <- key.match[[2]][idx]
    out[is.na(idx)] <- missing
    return(out)
  }
  return(rep(missing, length(terms)))
}

# Check if ggRandomForests loaded; if not, provide base-graphics fallbacks
has.ggRF <- requireNamespace("ggRandomForests", quietly = TRUE)

# Global graphical parameters
theme_set(theme_pubr(base_size = 10))
par(mfrow = c(1, 1), las = 1)
options(digits = 4, stringsAsFactors = FALSE, scipen = 0)

message("═══ Setup complete · Nature color scheme active ═══")

###############################################################################
# PART 1: Differential Gene Expression Analysis
# TOIL-recomputed TCGA-GTEx data · limma pipeline
###############################################################################

message("\n═══ PART 1: Differential Expression Analysis ═══")

# ── 1.1 Load phenotype and annotation data ──────────────────────────────────
Phenodata <- read.table(file.path(DATA_DIR, "TcgaTargetGTEX_phenotype.txt"),
                        sep = "\t", header = TRUE, stringsAsFactors = FALSE,
                        col.names = c("sample", "category", "tissue", "site",
                                      "type", "gender", "study"))

# GENCODE v23 annotations
Annotation23 <- data.table::fread(file.path(DATA_DIR, "gencode_v23_RNAs_annotation.txt"),
                                   sep = ";", header = FALSE, data.table = FALSE)[, 1:3]
colnames(Annotation23) <- c("ENSG_id", "Type", "Symbol")

# lncRNA v22 annotations
lncRNA_v22 <- read.table(file.path(DATA_DIR, "lnc_RNAs_annotation_v22.txt"),
                          strip.white = TRUE, blank.lines.skip = TRUE,
                          header = FALSE, stringsAsFactors = FALSE, sep = ";")[, c(1:3)]
lncRNA_v22 <- cbind(do.call(rbind, stringr::str_split(lncRNA_v22$V1, " ")), lncRNA_v22)
lncRNA_v22 <- lncRNA_v22[, -6]
colnames(lncRNA_v22) <- c("Chr", "start", "end", "strand", "ENSG_id", "Type", "Symbol")
lncRNA_v22$Position <- paste(lncRNA_v22$Chr, ":", lncRNA_v22$start, "-",
                              lncRNA_v22$end, " ", lncRNA_v22$strand)
lncRNA_v22 <- lncRNA_v22 %>% dplyr::select(ENSG_id, Symbol, Position, Type)

# lncRNA v23 annotations
lncRNA_v23 <- read.table(file.path(DATA_DIR, "lnc_RNAs_annotation_v23.txt"),
                          strip.white = TRUE, blank.lines.skip = TRUE,
                          header = FALSE, stringsAsFactors = FALSE, sep = ";")[, c(1:3)]
lncRNA_v23 <- cbind(do.call(rbind, stringr::str_split(lncRNA_v23$V1, " ")), lncRNA_v23)
lncRNA_v23 <- lncRNA_v23[, -6]
colnames(lncRNA_v23) <- c("Chr", "start", "end", "strand", "ENSG_id", "Type", "Symbol")
lncRNA_v23$Position <- paste(lncRNA_v23$Chr, ":", lncRNA_v23$start, "-",
                              lncRNA_v23$end, " ", lncRNA_v23$strand)
lncRNA_v23 <- lncRNA_v23 %>% dplyr::select(ENSG_id, Symbol, Position, Type)

message(sprintf("  lncRNA v22: %d  ·  lncRNA v23: %d",
                nrow(lncRNA_v22), nrow(lncRNA_v23)))

# ── 1.2 Extract colon samples from TCGA-GTEx ────────────────────────────────
colon.sample <- dplyr::filter(Phenodata, site == "Colon")
message(sprintf("  Colon samples: TCGA=%d  GTEX=%d  Total=%d",
                sum(colon.sample$study == "TCGA"),
                sum(colon.sample$study == "GTEX"),
                nrow(colon.sample)))

# Load expression matrix (TOIL-recomputed, log2(FPKM+0.001))
expr <- data.table::fread(file.path(DATA_DIR, "TcgaTargetGtex_rsem_gene_fpkm"),
                           sep = "\t", header = TRUE, data.table = FALSE)
rownames(expr) <- expr$sample
expr <- expr[, colnames(expr) %in% colon.sample$sample]
dim(expr)
message(sprintf("  Expression matrix: %d genes × %d samples", nrow(expr), ncol(expr)))

# ── 1.3 Filter low-expression genes ─────────────────────────────────────────
filter_fun <- function(x) { length(x[x < 0.1]) / length(x) > 0.6 }
expr <- expr[apply(expr, 1, filter_fun) == FALSE, ]
rownames(expr) <- gsub("-", ".", rownames(expr))
colnames(expr) <- gsub("-", ".", colnames(expr))
message(sprintf("  After filtering: %d genes", nrow(expr)))

# ── 1.4 Define groups ───────────────────────────────────────────────────────
group <- ifelse(grepl(".11$|GTEX", colnames(expr)), "NonTumor", "Tumor")
message(sprintf("  NonTumor=%d  Tumor=%d", sum(group == "NonTumor"), sum(group == "Tumor")))

# ── 1.5 limma differential expression ──────────────────────────────────────
design <- model.matrix(~ 0 + factor(group))
colnames(design) <- levels(factor(group))
contrast_matrix <- limma::makeContrasts(Tumor - NonTumor, levels = design)

fit <- lmFit(expr, design)
fit <- contrasts.fit(fit, contrast_matrix)
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
DEA <- topTable(fit, adjust.method = "BH", number = Inf, sort.by = "p")

# Map ENSG → Symbol
DEA$Symbol <- lookup(rownames(DEA), Annotation23[, -2])
write.table(DEA, file.path(OUT_DIR, "limma_DEA.txt"), quote = FALSE)

# lncRNA subset
lncRNA_DEA <- DEA[rownames(DEA) %in% lncRNA_v23$ENSG_id, ]
lncRNA_DEA$ENSG_id <- rownames(lncRNA_DEA)
write.table(lncRNA_DEA, file.path(OUT_DIR, "TCGA_GTEX_lncRNA_DEA.txt"), quote = FALSE)
message(sprintf("  DEA complete: %d lncRNAs differentially expressed", nrow(lncRNA_DEA)))

###############################################################################
# PART 2: Clinical Data & Survival Analysis
# TCGA-COAD raw data · Cox PH · Elastic Net feature selection
###############################################################################

message("\n═══ PART 2: Clinical & Survival Analysis ═══")

# ── 2.1 Load clinical data ──────────────────────────────────────────────────
clinical <- read.csv(file.path(DATA_DIR, "tcga_coad_clinical_followup_data.csv"),
                     header = TRUE, stringsAsFactors = FALSE)

clinical <- clinical %>% dplyr::select(-lymphnodes.positive.by.HE) %>% mutate(
  stage = factor(stage, levels = c("i", "ii", "iii", "iv"),
                 labels = c("i", "ii", "iii", "iv")),
  vital.status = factor(vital.status, levels = c("alive", "dead"),
                        labels = c("alive", "dead")),
  venous.invasion = factor(venous.invasion, levels = c("no", "yes"),
                           labels = c("no", "yes")),
  perineural.invasion = factor(perineural.invasion, levels = c("no", "yes"),
                               labels = c("no", "yes")),
  person.neoplasm.status = factor(person.neoplasm.status,
                                   levels = c("tumor free", "with tumor"),
                                   labels = c("tumor free", "with tumor")),
  lymphatic.invasion = factor(lymphatic.invasion, levels = c("no", "yes"),
                              labels = c("no", "yes")),
  gender = factor(gender, levels = c("male", "female"),
                  labels = c("male", "female")))

message(sprintf("  Clinical data: %d patients", nrow(clinical)))

# ── 2.2 Survival by stage ───────────────────────────────────────────────────
clinicdata <- clinical %>% dplyr::select(id, stage, OS, status, gender, age)
fit_stage <- survfit(Surv(OS / 365, status) ~ stage, data = clinicdata)

cairo_pdf(file.path(OUT_DIR, "Fig_S1_survival_by_stage.pdf"), width = 6, height = 5)
print(
  ggsurvplot(fit_stage, data = clinicdata, risk.table = FALSE, pval = TRUE,
             ggtheme = theme_survminer(),
             title = "Overall survival of stage group in TCGA-COAD",
             xlab = "Time (years)", legend = "top",
             surv.median.line = "hv",
             palette = ggsci::pal_npg()(4))
)
dev.off()
message("  ✓ Fig S1: Survival by stage")

# ── 2.3 Extract lncRNA expression from TCGA-COAD HTSeq-FPKM ─────────────────
Exprs <- data.table::fread(file.path(DATA_DIR, "COAD_HTSeq_FPKM.txt"),
                            sep = "\t", header = TRUE, data.table = FALSE,
                            stringsAsFactors = FALSE)
colnames(Exprs) <- gsub("-", ".", colnames(Exprs))

Annotation22 <- read.table(file.path(DATA_DIR, "gencode_v22_RNAs_annotation.txt"),
                            strip.white = TRUE, blank.lines.skip = TRUE,
                            header = FALSE, stringsAsFactors = FALSE, sep = ";")[, c(1:3)]
colnames(Annotation22) <- c("ENSG_id", "Type", "Symbol")

rownames(Exprs) <- Exprs$ID
Exprs <- Exprs[, -1]
Exprs <- as.data.frame(log2(as.matrix(Exprs) + 0.001))
Exprs <- Exprs[, -grep("\\.11", colnames(Exprs))]  # Remove non-tumor

# Collapse technical replicates
colnames(Exprs) <- substr(colnames(Exprs), 1, 12)
Exprs <- as.data.frame(limma::avearrays(Exprs))
message(sprintf("  TCGA-COAD tumor samples: %d", ncol(Exprs)))

# Filter low expression
Exprs <- Exprs[apply(Exprs, 1, filter_fun) == FALSE, ]

# lncRNA subset
lncRNA_Exprs <- Exprs[rownames(Exprs) %in% lncRNA_v22$ENSG_id, ]
lncRNA_Exprs <- as.data.frame(t(lncRNA_Exprs))
lncRNA_Exprs <- data.frame(id = rownames(lncRNA_Exprs), lncRNA_Exprs)
message(sprintf("  lncRNAs for survival: %d", ncol(lncRNA_Exprs) - 1))

# ── 2.4 Merge clinical + expression ─────────────────────────────────────────
clinicdata <- clinical %>% dplyr::select(id, stage, OS, status)
clinicdata <- na.omit(clinicdata)

Clinic_lncRNA_Exprs <- merge(clinicdata, lncRNA_Exprs, by = "id")
write.table(Clinic_lncRNA_Exprs, file.path(OUT_DIR, "Clinic_lncRNA_Exprs.txt"), quote = FALSE)
message(sprintf("  Merged data: %d patients × %d variables", nrow(Clinic_lncRNA_Exprs),
                ncol(Clinic_lncRNA_Exprs)))

# ── 2.5 Univariate Cox PH + Log-rank ────────────────────────────────────────
kaplan_dat <- Clinic_lncRNA_Exprs[, !colnames(Clinic_lncRNA_Exprs) %in%
                                    c("id", "stage", "status", "OS")]
OS_vec     <- Clinic_lncRNA_Exprs$OS
Status_vec <- as.numeric(as.vector(Clinic_lncRNA_Exprs$status))

# Batch log-rank
batch_log_rank_os <- apply(kaplan_dat, 2, function(values) {
  grp <- ifelse(values > median(values), "high", "low")
  dat <- data.frame(group = grp, time = OS_vec, status = Status_vec,
                    stringsAsFactors = FALSE)
  fit <- surv_fit(Surv(time, status) ~ group, data = dat)
  surv_pvalue(fit)
})

log_rank_os <- do.call(rbind, batch_log_rank_os)
log_rank_os$ENSG_id <- rownames(log_rank_os)
log_rank_os$Symbol  <- lookup(log_rank_os$ENSG_id, lncRNA_v22[, -c(3, 4)])
log_rank_os <- log_rank_os %>%
  dplyr::select(ENSG_id, pval, Symbol) %>%
  mutate(logrank.pvalue = pval) %>%
  dplyr::select(logrank.pvalue, ENSG_id, Symbol)
log_rank_os <- na.omit(log_rank_os)
write.table(log_rank_os, file.path(OUT_DIR, "log_rank_os.txt"), quote = FALSE)

# Univariate Cox PH (batch)
covariates <- colnames(Clinic_lncRNA_Exprs)[!colnames(Clinic_lncRNA_Exprs) %in%
                                              c("id", "status", "OS")]
coxdata <- Clinic_lncRNA_Exprs %>% dplyr::select(-id)

cox_os_formulas <- sapply(covariates,
  function(x) as.formula(paste("Surv(OS, status)~", x)))
cox_os_models <- lapply(cox_os_formulas,
  function(x) { coxph(x, data = coxdata) })

# Extract stage model separately (categorical)
cox_stage <- cox_os_models[[1]]
cox_os_models[[1]] <- NULL

coxph_os <- lapply(cox_os_models, function(x) {
  s <- summary(x)
  res <- c(beta = signif(s$coef[1], 2),
           HR   = signif(s$coef[2], 2),
           CI   = paste0("(", signif(s$conf.int[, "lower .95"], 3),
                         "-", signif(s$conf.int[, "upper .95"], 3), ")"),
           wald.test = signif(s$wald["test"], 3),
           HR.pvalue = signif(s$wald["pvalue"], 3))
  return(res)
})

coxph_os_res <- as.data.frame(do.call(rbind, coxph_os), stringsAsFactors = FALSE)
coxph_os_res <- data.frame(ENSG_id = rownames(coxph_os_res), coxph_os_res)
# Fix types: rbind coerced all to character
coxph_os_res$beta      <- as.numeric(coxph_os_res$beta)
coxph_os_res$HR        <- as.numeric(coxph_os_res$HR)
coxph_os_res$wald.test <- as.numeric(coxph_os_res$wald.test)
coxph_os_res$HR.pvalue <- as.numeric(coxph_os_res$HR.pvalue)
coxph_os_res$Symbol <- lookup(coxph_os_res$ENSG_id, lncRNA_v22)
coxph_os_res <- na.omit(coxph_os_res)
write.csv(coxph_os_res, file.path(OUT_DIR, "coxph_os.csv"), row.names = FALSE)

# Merge DEA + log-rank + Cox results
lncRNA_DEA_sel <- lncRNA_DEA %>% dplyr::select(logFC, adj.P.Val, Symbol, ENSG_id)
DEA_OS_analysis <- merge(
  merge(lncRNA_DEA_sel, log_rank_os, by = "ENSG_id"),
  coxph_os_res, by = "ENSG_id")
DEA_OS_analysis <- DEA_OS_analysis %>% dplyr::select(-Symbol.x, -Symbol.y)
write.csv(DEA_OS_analysis, file.path(OUT_DIR, "DEA_OS_analysis.csv"), row.names = FALSE)
message(sprintf("  DEA-OS merged: %d lncRNAs", nrow(DEA_OS_analysis)))

###############################################################################
# PART 3: Elastic Net Feature Selection & Risk Prediction
###############################################################################

message("\n═══ PART 3: Elastic Net Feature Selection & Risk Models ═══")

# ── 3.1 Prepare survival data ───────────────────────────────────────────────
survdata <- Clinic_lncRNA_Exprs[, colnames(Clinic_lncRNA_Exprs) %in%
  c("id", "OS", "stage", "status", DEA_OS_analysis$ENSG_id)]
survdata <- survdata %>% mutate(status = as.numeric(as.vector(status)))
message(sprintf("  Survival matrix: %d patients × %d lncRNAs",
                nrow(survdata), ncol(survdata) - 4))

# ── 3.2 Elastic Net tuning (grid search over α) ────────────────────────────
# Replace c060::EPSGO (broken in R 4.6) with manual alpha grid search
set.seed(12345)
alpha_grid <- seq(0, 1, by = 0.1)
cv_results <- list()
best_deviance <- Inf
best_alpha <- 0.5
best_lambda <- 0.1

cl <- makePSOCKcluster(min(3, parallel::detectCores() - 1))
registerDoParallel(cl)

for (a in alpha_grid) {
  set.seed(12345)
  fit_try <- tryCatch(
    cv.glmnet(y = Surv(survdata$OS, survdata$status),
              x = as.matrix(survdata[, 5:ncol(survdata)]),
              family = "cox", standardize = TRUE,
              alpha = a, nlambda = 50, nfolds = 10),
    error = function(e) NULL)
  if (!is.null(fit_try)) {
    cv_results[[as.character(a)]] <- fit_try
    min_dev <- min(fit_try$cvm)
    if (min_dev < best_deviance) {
      best_deviance <- min_dev
      best_alpha <- a
      best_lambda <- fit_try$lambda.min
    }
  }
}
stopCluster(cl)

opt_alpha  <- best_alpha
opt_lambda <- best_lambda
opt_error  <- best_deviance

message(sprintf("  Optimal: α=%.3f  λ=%.3f  deviance=%.3f",
                opt_alpha, opt_lambda, opt_error))

# ── 3.3 Cross-validation deviance surface ───────────────────────────────────
# Plot deviance curves for each alpha
cairo_pdf(file.path(OUT_DIR, "Fig_3a_crossvalidation_deviance.pdf"), width = 8, height = 6)
plot(cv_results[["0"]]$lambda, cv_results[["0"]]$cvm, type = "n",
     log = "x", ylim = range(sapply(cv_results, function(f) range(f$cvm))),
     xlab = expression(log(lambda)), ylab = "Partial Likelihood Deviance",
     main = "Elastic Net CV Deviance by Alpha")
alphas_used <- names(cv_results)
for (i in seq_along(cv_results)) {
  lines(cv_results[[i]]$lambda, cv_results[[i]]$cvm,
        col = colorRampPalette(c(nature_blue, nature_red))(length(cv_results))[i],
        lwd = 1.5)
}
legend("topright", legend = paste0("α=", alphas_used),
       col = colorRampPalette(c(nature_blue, nature_red))(length(cv_results)),
       lwd = 1.5, cex = 0.7, bty = "n", title = "Alpha")
abline(v = opt_lambda, lty = 2, col = nature_red)
dev.off()

# 3D deviance surface (using first few alphas for surface plot)
if (length(cv_results) >= 3) {
  surf_data <- do.call(rbind, lapply(names(cv_results), function(a_name) {
    f <- cv_results[[a_name]]
    data.frame(alpha = as.numeric(a_name), lambda = f$lambda, deviance = f$cvm)
  }))
  deviance_plot <- plotly::plot_ly(
    x = ~surf_data$alpha, y = ~surf_data$lambda, z = ~surf_data$deviance) %>%
    add_markers(color = ~surf_data$deviance,
                colors = colorRamp(c(nature_blue, nature_orange, nature_red)))
  htmlwidgets::saveWidget(deviance_plot,
    file.path(OUT_DIR, "deviance_plot.html"))
}
message("  ✓ 3D deviance surface")

# ── 3.4 Fit final Elastic Net Cox model ─────────────────────────────────────
set.seed(12345)
cv.glmfit <- cv.glmnet(
  y = Surv(survdata$OS, survdata$status),
  x = as.matrix(survdata[, 5:ncol(survdata)]),
  family = "cox", standardize = TRUE,
  alpha = opt_alpha, nlambda = 50)

# Plot: Elastic Net deviance curve
cairo_pdf(file.path(OUT_DIR, "Fig_3b_elastic_net_deviance.pdf"), width = 6, height = 5)
plot(cv.glmfit, main = "Elastic Net Regularized CoxPH Regression",
     xlab = expression(log(lambda)))
abline(h = opt_error, lty = 3, col = nature_red)
dev.off()
message("  ✓ Elastic Net deviance curve")

# ── 3.5 Extract non-zero coefficients ───────────────────────────────────────
	res.cv <- cv.glmfit$glmnet.fit
	# Extract coefficients as named vector for plotCoef compatibility
	cof_raw <- coef(res.cv, s = cv.glmfit$lambda.min)
	cof_vec <- as.numeric(cof_raw)
	names(cof_vec) <- rownames(cof_raw)
	cof_vec <- cof_vec[cof_vec != 0]
	# Also build dataframe for downstream use
	cof <- as.data.frame(as.matrix(cof_raw))
	cof$ENSG_id <- row.names(cof)
	cof$Symbol  <- lookup(cof$ENSG_id, lncRNA_v22)
	colnames(cof) <- c("beta", "ENSG_id", "Symbol")
	
	nonozero.coef <- dplyr::filter(cof, beta != 0) %>% arrange(desc(beta))
	message(sprintf("  %d lncRNAs with non-zero coefficients", nrow(nonozero.coef)))
	
	# Coefficient shrinkage paths -- use named vector from cof_vec
	bet <- res.cv$beta[match(names(cof_vec), rownames(res.cv$beta)), ]

cairo_pdf(file.path(OUT_DIR, "Fig_3c_coefficient_paths.pdf"), width = 6, height = 5)
glmnet:::plotCoef(bet, lambda = res.cv$lambda, df = res.cv$df,
  dev = res.cv$dev.ratio, xvar = "lambda", add = FALSE,
  col = nature_blue, label = TRUE,
  main = "Coefficient shrinkage paths", xlab = expression(log(lambda)))
abline(v = log(cv.glmfit$lambda.min), lty = 3, col = nature_red)
abline(v = log(cv.glmfit$lambda.1se), lty = 3, col = nature_grey)
dev.off()

cairo_pdf(file.path(OUT_DIR, "Fig_3d_fraction_deviance.pdf"), width = 6, height = 5)
glmnet:::plotCoef(bet, lambda = res.cv$lambda, df = res.cv$df,
  dev = res.cv$dev.ratio, xvar = "dev", add = FALSE,
  col = nature_blue, label = TRUE,
  main = "Coefficient shrinkage by Fraction Deviance Explained",
  xlab = "Fraction Deviance Explained")
abline(v = 0.059, lty = 3, col = nature_red)
dev.off()
message("  ✓ Coefficient path plots")

# ── 3.6 Feature summary plot (Nature style) ─────────────────────────────────
ggtextdata <- DEA_OS_analysis[DEA_OS_analysis$ENSG_id %in% nonozero.coef$ENSG_id, ]
ggtextdata$Type     <- lookup(ggtextdata$ENSG_id, lncRNA_v22[, -c(2, 3)])
ggtextdata$Position <- lookup(ggtextdata$ENSG_id, lncRNA_v22[, c(1, 3)])

cairo_pdf(file.path(OUT_DIR, "Fig_2_feature_lncRNAs.pdf"), width = 8, height = 6)
print(
  ggdotchart(ggtextdata, x = "Symbol", y = "logFC", dot.size = "HR",
             color = "HR.pvalue", add = "segments", sorting = "descending",
             ggtheme = theme_pubr()) +
    geom_hline(yintercept = 0, linetype = 2, color = nature_grey) +
    scale_color_gradientn(colors = c(nature_red, nature_orange, nature_blue)) +
    theme(legend.position = "right") +
    labs(title = "Features of candidate lncRNAs") +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) +
    geom_rug(alpha = 0.8, size = 0.5, color = nature_grey)
)
dev.off()
message("  ✓ Feature lncRNA dot chart")

# ── 3.7 Train/Test split ────────────────────────────────────────────────────
survdata$group <- ifelse(survdata$status == "0", "alive", "dead")
survdata <- survdata %>% mutate(group = factor(group))
set.seed(12345)
trainIndex <- createDataPartition(survdata$group, p = 0.7, list = FALSE, times = 1)

Train <- survdata[trainIndex, ] %>% dplyr::select(-group)
Test  <- survdata[-trainIndex, ] %>% dplyr::select(-group)

Train <- Train[, colnames(Train) %in% c("id", "OS", "stage", "status", nonozero.coef$ENSG_id)]
Test  <- Test[, colnames(Test) %in% c("id", "OS", "stage", "status", nonozero.coef$ENSG_id)]
message(sprintf("  Train=%d  Test=%d  Features=%d", nrow(Train), nrow(Test), ncol(Train) - 4))

# ── 3.8 Multivariate Cox PH (stepwise) ──────────────────────────────────────
cox_data_train <- Train %>% dplyr::select(-id) %>% mutate(status = as.numeric(as.vector(status)))
cox_data_test  <- Test  %>% dplyr::select(-id) %>% mutate(status = as.numeric(as.vector(status)))

cox_full <- coxph(Surv(OS, status) ~ ., data = cox_data_train)
stepcox  <- step(cox_full, trace = 0, direction = "both", steps = 1000, k = 2)
cox_vip_list <- names(stepcox$means)
message(sprintf("  Stepwise Cox: %d features, AIC=%.2f",
                length(cox_vip_list),
                AIC(stepcox)))

# Forest plot (Nature style)
cairo_pdf(file.path(OUT_DIR, "Fig_3_forestplot_cox.pdf"), width = 7, height = 5)
print(
  ggforest(stepcox, main = "Feature importance in Overall Survival (TCGA-COAD)",
           data = cox_data_train)
)
dev.off()
message("  ✓ Forest plot")

# ── 3.9 Time-dependent survival ROC ──────────────────────────────────────────
PlotsurvROC_nature <- function(predict.time, Train, Test, model = stepcox) {
  newdata_t <- Train[, colnames(Train) %in% c(cox_vip_list, "status", "stage", "OS", "id")]
  riskScore_t <- predict(model, type = "risk", newdata = newdata_t)
  risk_df_t <- cbind.data.frame(newdata_t[, c("id", "OS", "status")],
                                 riskScore = riskScore_t)

  roc_t <- survivalROC(Stime = risk_df_t$OS / 365, status = risk_df_t$status,
                        marker = risk_df_t$riskScore, predict.time = predict.time,
                        method = "KM")

  plot(roc_t$FP, roc_t$TP, type = "l", xlim = c(0, 1), ylim = c(0, 1),
       col = nature_blue, lwd = 2, xlab = "False positive rate",
       ylab = "True positive rate",
       main = paste("ROC of CoxPH Risk Prediction (Year =", predict.time, ")"))

  abline(0, 1, col = nature_grey, lty = 2)

  newdata_v <- Test[, colnames(Test) %in% c(cox_vip_list, "status", "stage", "OS", "id")]
  riskScore_v <- predict(model, type = "risk", newdata = newdata_v)
  risk_df_v <- cbind.data.frame(newdata_v[, c("id", "OS", "status")],
                                 riskScore = riskScore_v)

  roc_v <- survivalROC(Stime = risk_df_v$OS / 365, status = risk_df_v$status,
                        marker = risk_df_v$riskScore, predict.time = predict.time,
                        method = "KM")

  lines(roc_v$FP, roc_v$TP, type = "l", col = nature_red, lwd = 2)
  legend("bottomright", col = c(nature_blue, nature_red), bty = "n",
         cex = 0.8, lwd = 2, title = "ROC AUC",
         legend = c(paste("Training:", round(roc_t$AUC, 3)),
                    paste("Testing:",  round(roc_v$AUC, 3))))
}

for (yr in c(1, 3, 5)) {
  cairo_pdf(file.path(OUT_DIR, sprintf("Fig_4_roc_%dyear.pdf", yr)), width = 5, height = 5)
  PlotsurvROC_nature(yr, Train, Test, stepcox)
  dev.off()
}
message("  ✓ Time-dependent ROC curves (1/3/5 year)")

# ── 3.10 AUC over time ──────────────────────────────────────────────────────
CALAUC_train <- function(time) {
  newdata <- Train[, colnames(Train) %in% c(cox_vip_list, "status", "stage", "OS", "id")]
  riskScore <- predict(stepcox, type = "risk", newdata = newdata)
  risk_df <- cbind.data.frame(newdata[, c("id", "OS", "status")], riskScore = riskScore)
  roc <- survivalROC(Stime = risk_df$OS / 365, status = risk_df$status,
                      marker = risk_df$riskScore, predict.time = time, method = "KM")
  roc$AUC
}
CALAUC_test <- function(time) {
  newdata <- Test[, colnames(Test) %in% c(cox_vip_list, "status", "stage", "OS", "id")]
  riskScore <- predict(stepcox, type = "risk", newdata = newdata)
  risk_df <- cbind.data.frame(newdata[, c("id", "OS", "status")], riskScore = riskScore)
  roc <- survivalROC(Stime = risk_df$OS / 365, status = risk_df$status,
                      marker = risk_df$riskScore, predict.time = time, method = "KM")
  roc$AUC
}

AUC <- data.frame(
  time = rep(seq(1, 15, 0.05), 2),
  group = rep(c("Training", "Testing"), each = length(seq(1, 15, 0.05))),
  AUC = c(sapply(seq(1, 15, 0.05), CALAUC_train),
          sapply(seq(1, 15, 0.05), CALAUC_test)))

cairo_pdf(file.path(OUT_DIR, "Fig_5_auc_over_time.pdf"), width = 6, height = 5)
print(
  ggplot(AUC, aes(x = time, y = AUC, colour = group)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = pal_train_test) +
    labs(title = "CoxPH Risk Prediction Performance Over Time",
         y = "Area Under Curve (AUC)", x = "Observation Time (years)") +
    geom_vline(xintercept = c(1.5, 3, 5, 7, 9),
               linetype = "dashed", color = nature_grey, size = 0.3) +
    theme_pubr() +
    theme(plot.title = element_text(size = 14, hjust = 0.5))
)
dev.off()
message("  ✓ AUC over time plot")

# ── 3.11 Risk score distribution ────────────────────────────────────────────
GGsurvplot_nature <- function(lncrna_list, dataset, stepcox) {
  newdata <- dataset[, colnames(dataset) %in% c(lncrna_list, "status", "stage", "OS", "id")]
  riskScore <- predict(stepcox, type = "risk", newdata = newdata)
  risk <- as.factor(ifelse(riskScore > median(riskScore), "High", "Low"))
  risk_df <- cbind.data.frame(newdata[, c("id", "OS", "status")],
                               riskScore = riskScore, risk = risk)
  risk_df <- risk_df %>% mutate(status = as.numeric(as.vector(status)))
  fit <- survfit(Surv(OS / 365, status) ~ risk, data = risk_df)
  ggsurvplot(fit, data = risk_df, risk.table = FALSE, pval = TRUE,
             ggtheme = theme_survminer(),
             title = "OS probability by predicted risk group",
             xlab = "Time (years)", legend = "top",
             surv.median.line = "hv",
             palette = pal_risk)
}

for (label in c("Train", "Test")) {
  ds <- if (label == "Train") Train else Test
  cairo_pdf(file.path(OUT_DIR, sprintf("Fig_6_survival_%s.pdf", tolower(label))),
            width = 5, height = 5)
  print(GGsurvplot_nature(cox_vip_list, ds, stepcox))
  dev.off()
}
message("  ✓ Risk-stratified KM curves")

# ── 3.12 Risk score distributions ───────────────────────────────────────────
risk_fun <- function(dataset) {
  newdata <- dataset[, colnames(dataset) %in% c(cox_vip_list, "status", "stage", "OS", "id")]
  riskScore <- predict(stepcox, type = "risk", newdata = newdata, reference = "strata")
  risk <- as.factor(ifelse(riskScore > median(riskScore), "High", "Low"))
  risk_df <- cbind.data.frame(newdata[, c("id", "OS", "status")],
                               riskScore = riskScore, risk = risk)
  risk_df <- risk_df %>% mutate(OS = as.numeric(OS), riskScore = as.numeric(riskScore)) %>%
    arrange(desc(riskScore)) %>%
    mutate(log10.riskScore = log10(riskScore),
           zscore = (riskScore - mean(riskScore)) / sd(riskScore))
  risk_df$survival.status <- ifelse(risk_df$status == 0, "Censored", "Dead")
  risk_df$survival.status <- as.factor(risk_df$survival.status)
  return(risk_df)
}

risk_df <- rbind(
  data.frame(risk_fun(Train), group = "Training"),
  data.frame(risk_fun(Test),  group = "Testing"))

cairo_pdf(file.path(OUT_DIR, "Fig_7_riskscore_comparison.pdf"), width = 6, height = 5)
print(
  ggboxplot(risk_df, x = "group", y = "log10.riskScore",
            color = "survival.status", palette = pal_status,
            add = "jitter") +
    stat_compare_means(aes(group = survival.status),
                       method = "t.test", label = "p.signif") +
    labs(title = "Risk score comparison by group")
)
dev.off()

cairo_pdf(file.path(OUT_DIR, "Fig_8_riskscore_distribution.pdf"), width = 8, height = 6)
print(
  ggdotchart(risk_df, x = "id", y = "log10.riskScore",
             shape = "group", color = "survival.status",
             title = "Risk score distribution of TCGA-COAD patients",
             repel = TRUE, ylab = "log10(RiskScore)",
             size = 1, palette = pal_status,
             ggtheme = theme_pubr()) +
    coord_flip() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.title = element_text(size = 14, hjust = 0.5)) +
    geom_hline(yintercept = median(risk_df$log10.riskScore),
               linetype = "dashed", color = nature_grey) +
    geom_rug(aes(color = survival.status), data = risk_df,
             size = 0.1, alpha = 1, position = "jitter") +
    scale_shape_manual(values = c(18, 20))
)
dev.off()
message("  ✓ Risk score distribution plots")

###############################################################################
# PART 4: Random Forest Survival Analysis
###############################################################################

message("\n═══ PART 4: Random Forest Survival ═══")

# ── 4.1 Prepare data ────────────────────────────────────────────────────────
rfdata_train <- Train %>% dplyr::select(-id) %>% mutate(OS = OS / 365)
rfdata_train <- rfdata_train[, colnames(rfdata_train) %in%
  c("stage", "OS", "status", cox_vip_list)]

rfdata_test <- Test %>% dplyr::select(-id) %>% mutate(OS = OS / 365)
rfdata_test <- rfdata_test[, colnames(rfdata_test) %in%
  c("stage", "OS", "status", cox_vip_list)]

# ── 4.2 Hyperparameter tuning ───────────────────────────────────────────────
set.seed(12345)
options(rf.cores = 4)
rf_tune <- tune(Surv(OS, status) ~ ., data = rfdata_train,
  mtryStart = floor(sqrt(ncol(rfdata_train))),
  nodesizeTry = c(1:9, seq(10, 100, by = 5)), ntreeTry = 500,
  stepFactor = 1.25, improve = 1e-3, strikeout = 3, maxIter = 50,
  trace = FALSE, doBest = TRUE)
message(sprintf("  Optimal RF params: mtry=%d nodesize=%d",
                rf_tune$optimal["mtry"], rf_tune$optimal["nodesize"]))

# OOB error surface
rftun_res <- as.data.frame(rf_tune[["results"]])
OOBerror  <- as.matrix(as.data.frame.matrix(
  xtabs(err ~ nodesize + mtry, data = rftun_res)))

OOBerror_plot <- plotly::plot_ly(
  x = ~rftun_res$mtry, y = ~rftun_res$nodesize, z = ~OOBerror) %>%
  add_surface(colorscale = list(
    c(0, nature_blue), c(0.5, nature_orange), c(1, nature_red)))
htmlwidgets::saveWidget(OOBerror_plot,
  file.path(OUT_DIR, "OOBerror_plot.html"))
message("  ✓ OOB error surface")

# ── 4.3 Fit Random Forest Survival ──────────────────────────────────────────
set.seed(12345)
options(rf.cores = 4)
rfsrc <- rfsrc(Surv(OS, status) ~ ., data = rfdata_train,
  nsplit = 3, mtry = rf_tune$optimal["mtry"],
  nodesize = rf_tune$optimal["nodesize"],
  bootstrap = "by.root", samptype = "swor",
  na.action = "na.impute", nimpute = 1,
  tree.err = TRUE, split.depth = "by.tree",
  importance = "permute", ntree = 500,
  do.trace = FALSE, statistics = TRUE)

rf.pred <- predict(rfsrc, newdata = rfdata_test, na.action = "na.impute")

# ── 4.4 RF variable importance ──────────────────────────────────────────────
cairo_pdf(file.path(OUT_DIR, "Fig_S2_rf_vimp.pdf"), width = 7, height = 5)
plot(rfsrc, main = "Random Forest OOB Error and Variable Importance")
dev.off()

cairo_pdf(file.path(OUT_DIR, "Fig_S3_rf_vimp_bar.pdf"), width = 6, height = 5)
if (has.ggRF) {
  print(
    plot(gg_vimp(rfsrc)) +
      theme(legend.position = "top") +
      labs(fill = "VIMP > 0", title = "Variable Importance for OS Prediction") +
      theme(plot.title = element_text(size = 14, hjust = 0.5))
  )
} else {
  # Base R fallback: variable importance from rfsrc
  vimp_data <- data.frame(
    Variable = names(rfsrc$importance),
    Importance = as.numeric(rfsrc$importance))
  vimp_data <- vimp_data[order(vimp_data$Importance, decreasing = TRUE), ]
  vimp_data <- head(vimp_data, 15)
  par(las = 2, mar = c(6, 8, 4, 2))
  barplot(vimp_data$Importance, names.arg = vimp_data$Variable,
          horiz = TRUE, col = nature_blue, border = NA,
          main = "Variable Importance for OS Prediction",
          xlab = "Variable Importance (VIMP)")
}
dev.off()

# ── 4.5 RF predicted survival ───────────────────────────────────────────────
for (label in c("training", "testing")) {
  obj  <- if (label == "training") rfsrc else rf.pred
  oob  <- round(mean(na.omit(obj$err.rate)) * 100, 3)
  cairo_pdf(file.path(OUT_DIR,
    sprintf("Fig_9_rf_survival_%s.pdf", label)), width = 5, height = 5)
  if (has.ggRF) {
    p <- plot(gg_rfsrc(obj)) +
      theme(legend.position = c(0.1, 0.2)) +
      labs(title = paste("RF predicted survival —", label),
           y = "Survival Probability", x = "Time (years)") +
      theme(plot.title = element_text(size = 14, hjust = 0)) +
      geom_vline(xintercept = c(1, 3), linetype = "dashed", color = nature_grey) +
      coord_cartesian(x = c(0, 4)) + theme_pubr() + geom_rug() +
      annotate("text", x = 0.85, y = 0.5, color = "black",
               label = paste("OOB error:", oob, "%")) +
      scale_color_manual(values = pal_status)
    print(p)
  } else {
    # Base R fallback: use built-in plot.rfsrc / plot.predict.rfsrc
    plot(obj, main = paste("RF predicted survival —", label))
    abline(v = c(1, 3), lty = 2, col = nature_grey)
    text(0.85, 0.5, paste("OOB error:", oob, "%"))
  }
  dev.off()
}
message("  ✓ RF survival prediction plots")

# ── 4.6 RF variable dependence ──────────────────────────────────────────────
if (has.ggRF) {
  gg_v <- gg_variable(rfsrc, time = c(1, 3, 5),
                       time.labels = c("1 Year", "3 Years", "5 Years"))

  cairo_pdf(file.path(OUT_DIR, "Fig_10_rf_stage_survival.pdf"), width = 5, height = 6)
  print(
    plot(gg_v, xvar = "stage", alpha = 0.6) +
      labs(title = "Survival Probability by Tumor Stage") +
      theme(legend.position = "top") +
      labs(y = "Survival", x = "Tumor Stage") +
      scale_color_manual(values = pal_status) +
      theme(plot.title = element_text(size = 14, hjust = 0.5))
  )
  dev.off()

  # lncRNA-dependent survival
  colnames(gg_v)[2:8] <- lookup(colnames(gg_v)[2:8], Annotation23[, -2])
  cairo_pdf(file.path(OUT_DIR, "Fig_11_rf_lncrna_survival.pdf"), width = 7, height = 6)
  print(
    plot(gg_v, xvar = colnames(gg_v)[2:8], panel = TRUE, alpha = 0.1) +
      theme(legend.position = "top") +
      labs(title = "LncRNA Expression-Dependent Survival Probability",
           y = "Survival Probability",
           x = "lncRNA Expression log2(FPKM+0.001)") +
      scale_color_manual(values = pal_status) +
      coord_cartesian(ylim = c(0.4, 1)) +
      theme(plot.title = element_text(size = 12, hjust = 0.5))
  )
  dev.off()
} else {
  # Base R fallback: plot.variable with single time point (v3.6+ restriction)
  for (t in c(1, 3, 5)) {
    cairo_pdf(file.path(OUT_DIR, sprintf("Fig_10_rf_stage_survival_%dy.pdf", t)),
              width = 6, height = 5)
    plot.variable(rfsrc, xvar.names = "stage", time = t,
                  surv.type = "surv", main = paste("Survival by Stage at", t, "Year(s)"))
    dev.off()
  }

  cairo_pdf(file.path(OUT_DIR, "Fig_11_rf_lncrna_survival.pdf"), width = 8, height = 6)
  lnc_vars <- intersect(cox_vip_list, colnames(rfsrc$xvar))
  if (length(lnc_vars) > 0) {
    plot.variable(rfsrc, xvar.names = lnc_vars[1], time = 3,
                  surv.type = "surv", partial = TRUE,
                  main = "LncRNA-Dependent Survival (3 Year)")
  }
  dev.off()
}
message("  ✓ RF variable dependence plots")

# ── 4.7 Brier Score & C-index ───────────────────────────────────────────────
# NOTE: pec package API changed; wrap in tryCatch to continue if unavailable
tryCatch({
  pecdata_train <- rfdata_train
  pecdata_test  <- rfdata_test

  Models_train <- list(
    "CoxPH (stage)" = coxph(Surv(OS, status) ~ stage, data = pecdata_train, x = TRUE, y = TRUE),
    "CoxPH (+lncRNA)" = coxph(Surv(OS, status) ~ ., data = pecdata_train, x = TRUE, y = TRUE),
    "Random Forest" = rfsrc)

  Bier_train <- pec::pec(Models_train,
    formula = Surv(OS, status) ~ ., data = pecdata_train,
    cens.model = "marginal", splitMethod = "bootcv",
    M = round(nrow(pecdata_train) * 0.6), B = 100,
    keep.index = TRUE, multiSplitTest = FALSE, confInt = FALSE,
    exact = TRUE, verbose = FALSE, maxtime = 3000,
    eval.times = seq(0, 3650, 30))

  Models_test <- list(
    "CoxPH (stage)" = coxph(Surv(OS, status) ~ stage, data = pecdata_test, x = TRUE, y = TRUE),
    "CoxPH (+lncRNA)" = coxph(Surv(OS, status) ~ ., data = pecdata_test, x = TRUE, y = TRUE),
    "Random Forest" = rfsrc)

  Bier_test <- pec::pec(Models_test,
    formula = Surv(OS, status) ~ ., data = pecdata_test,
    cens.model = "marginal", splitMethod = "bootcv",
    M = round(nrow(pecdata_test) * 0.6), B = 100,
    keep.index = TRUE, multiSplitTest = FALSE, confInt = FALSE,
    exact = TRUE, verbose = FALSE, maxtime = 3000,
    eval.times = seq(0, 3650, 30))

  model_colors <- c(nature_grey, nature_blue, nature_red)

  for (label in c("train", "test")) {
    obj <- if (label == "train") Bier_train else Bier_test
    cairo_pdf(file.path(OUT_DIR, sprintf("Fig_12_brier_%s.pdf", label)),
              width = 6, height = 5)
    plot(obj, smooth = TRUE, lwd = 1.5, legend.cex = 0.8,
         type = "s", add.refline = TRUE,
         xlim = c(0, 3650), ylim = c(0, 0.5),
         xlab = "Time (days)", ylab = "Prediction error (Brier score)",
         col = model_colors)
    dev.off()
  }

  # C-index
  Cindex_train <- pec::cindex(Models_train,
    formula = Surv(OS, status) ~ ., data = pecdata_train,
    eval.times = seq(0, 3650, 30),
    splitMethod = "bootcv", M = round(nrow(pecdata_train) * 0.6),
    cens.model = "marginal", B = 100, maxtime = 3000)

  Cindex_test <- pec::cindex(Models_test,
    formula = Surv(OS, status) ~ ., data = pecdata_test,
    eval.times = seq(0, 3650, 30),
    splitMethod = "bootcv", M = round(nrow(pecdata_test) * 0.6),
    cens.model = "marginal", B = 100, maxtime = 3000)

  for (label in c("train", "test")) {
    obj <- if (label == "train") Cindex_train else Cindex_test
    cairo_pdf(file.path(OUT_DIR, sprintf("Fig_13_cindex_%s.pdf", label)),
              width = 6, height = 5)
    plot(obj, smooth = TRUE, lwd = 1.5, legend.cex = 0.8,
         type = "s", add.refline = TRUE,
         xlim = c(0, 3650), ylim = c(0, 1),
         xlab = "Time (days)", ylab = "Concordance Index (C-index)",
         col = model_colors)
    dev.off()
  }
  message("  ✓ Brier score & C-index plots")
}, error = function(e) {
  message(sprintf("  ⚠ Brier/C-index skipped (pec error): %s", e$message))
})

###############################################################################
# PART 5: Classification Models
# SVM · Random Forest · Neural Network · Elastic Net · Logistic Regression
###############################################################################

message("\n═══ PART 5: Machine Learning Classification ═══")

# ── 5.1 Prepare classification data ─────────────────────────────────────────
classification_data <- expr[rownames(expr) %in% cox_vip_list, ]
classification_data <- as.data.frame(t(classification_data))
colnames(classification_data) <- lookup(colnames(classification_data),
                                         lncRNA_v23[, c(1, 2)])

classification_data <- data.frame(
  group = ifelse(grepl(".11$|GTEX", rownames(classification_data)),
                 "NonTumor", "Tumor"),
  classification_data)
classification_data$group <- factor(classification_data$group,
                                     levels = c("NonTumor", "Tumor"))

# Box-Cox normalization
normalization <- preProcess(classification_data, method = "BoxCox")
classification_data <- predict(normalization, newdata = classification_data)
write.csv(classification_data, file.path(OUT_DIR, "Machine_learning_data.csv"),
          row.names = FALSE)

# ── 5.2 Feature distribution plots (Nature style) ───────────────────────────
cairo_pdf(file.path(OUT_DIR, "Fig_S4_feature_density.pdf"), width = 10, height = 8)
featurePlot(x = classification_data[, -1], y = classification_data$group,
            plot = "density", center = FALSE, scale = FALSE,
            main = "Feature density distribution",
            scales = list(x = list(relation = "free"), y = list(relation = "free")),
            adjust = 0.5, pch = "|", layout = c(4, 2),
            auto.key = list(columns = 2))
dev.off()

cairo_pdf(file.path(OUT_DIR, "Fig_S5_feature_boxplot.pdf"), width = 10, height = 8)
featurePlot(x = classification_data[, -1], y = classification_data$group,
            plot = "box", center = FALSE, scale = FALSE,
            main = "Feature Box-Whisker Plot",
            scales = list(y = list(relation = "free"), x = list(rot = 90)),
            layout = c(4, 2), auto.key = list(columns = 2))
dev.off()
message("  ✓ Feature distribution plots")

# ── 5.3 Train/Test split ────────────────────────────────────────────────────
set.seed(12345)
ind <- createDataPartition(classification_data$group, p = 0.7, list = FALSE)
train <- classification_data[ind, ]
test  <- classification_data[-ind, ]

custom <- trainControl(method = "repeatedcv", number = 10, repeats = 5,
  verboseIter = FALSE, classProbs = TRUE, savePredictions = "final",
  summaryFunction = twoClassSummary)

# ── 5.4 Logistic Regression ─────────────────────────────────────────────────
set.seed(12345)
logisreg <- train(group ~ ., data = train, method = "glmStepAIC",
  family = "binomial", tuneLength = 10, trControl = custom)
message(sprintf("  Logistic Regression: Accuracy=%.3f",
                max(logisreg$results$Accuracy)))

# ── 5.5 Elastic Net ─────────────────────────────────────────────────────────
set.seed(12345)
elastnet <- train(group ~ ., data = train, method = "glmnet",
  tuneLength = 10, trControl = custom, metric = "ROC",
  tuneGrid = expand.grid(alpha = seq(0, 1, length = 10),
                          lambda = seq(0.0001, 0.1, length = 5)))
message(sprintf("  Elastic Net: best α=%.2f λ=%.4f",
                elastnet$bestTune$alpha, elastnet$bestTune$lambda))

# ── 5.6 Random Forest ───────────────────────────────────────────────────────
cl <- makePSOCKcluster(3)
registerDoParallel(cl)
set.seed(12345)
RF <- train(group ~ ., data = train, method = "ranger", tuneLength = 10,
  trControl = custom, metric = "ROC",
  tuneGrid = expand.grid(mtry = 1:(ncol(train) - 1),
                          splitrule = "gini", min.node.size = seq(1, 10)))
stopCluster(cl)
message(sprintf("  Random Forest: best mtry=%d", RF$bestTune$mtry))

# ── 5.7 SVM ─────────────────────────────────────────────────────────────────
cl <- makePSOCKcluster(3)
registerDoParallel(cl)
set.seed(12345)
SVM <- train(group ~ ., data = train, method = "svmRadial",
  tuneLength = 10, trControl = custom, metric = "ROC",
  tuneGrid = expand.grid(sigma = seq(0, 1, length = 10),
                          C = seq(0.1, 2, length = 10)))
stopCluster(cl)
message(sprintf("  SVM: best sigma=%.2f C=%.2f",
                SVM$bestTune$sigma, SVM$bestTune$C))

# ── 5.8 Neural Network ─────────────────────────────────────────────────────
cl <- makePSOCKcluster(3)
registerDoParallel(cl)
set.seed(12345)
NNET <- train(group ~ ., data = train, method = "nnet", tuneLength = 10,
  trControl = custom, metric = "ROC",
  tuneGrid = expand.grid(size = 1:15, decay = seq(0.1, 1, length = 5)))
stopCluster(cl)
message(sprintf("  Neural Network: best size=%d decay=%.2f",
                NNET$bestTune$size, NNET$bestTune$decay))

###############################################################################
# PART 6: Model Evaluation & Comparison
###############################################################################

message("\n═══ PART 6: Model Evaluation ═══")

# ── 6.1 Model performance metrics ───────────────────────────────────────────
model_evaluation <- function(models) {
  model_pred <- function(models, dataset) {
    predsclas <- predict(models, newdata = dataset)
    res <- confusionMatrix(data = predsclas, reference = dataset$group,
                           positive = "Tumor", mode = "prec_recall")
    data.frame(Accuracy  = res[["overall"]][["Accuracy"]],
               Lower     = res[["overall"]][["AccuracyLower"]],
               Upper     = res[["overall"]][["AccuracyUpper"]],
               AccPVal   = res[["overall"]][["AccuracyPValue"]],
               Kappa     = res[["overall"]][["Kappa"]],
               Precision = res[["byClass"]][["Precision"]],
               Recall    = res[["byClass"]][["Recall"]],
               F1        = res[["byClass"]][["F1"]])
  }
  model_ROC <- function(models, dataset) {
    pred <- predict(models, newdata = dataset, type = "prob")
    roc  <- roc(dataset$group, pred$Tumor, quiet = TRUE)
    as.numeric(auc(roc))
  }
  temp <- rbind(model_pred(models, train), model_pred(models, test))
  temp$DataSets <- c("Training", "Testing")
  temp$AUC <- c(model_ROC(models, train), model_ROC(models, test))
  temp <- dplyr::select(temp, DataSets, AUC, Accuracy, Lower, Upper,
                        AccPVal, Precision, Recall, F1, Kappa)
  return(temp)
}

Modelperformance <- rbind(
  model_evaluation(logisreg),
  model_evaluation(RF),
  model_evaluation(SVM),
  model_evaluation(elastnet),
  model_evaluation(NNET))
Modelperformance <- data.frame(
  Models = rep(c("Logistic", "RandomForest", "SVM", "ElasticNet", "NeuralNet"),
               each = 2),
  Modelperformance)
write.csv(Modelperformance, file.path(OUT_DIR, "model_performance.csv"),
          row.names = FALSE)
message("  ✓ Model performance table saved")

# ── 6.2 ROC curves (Nature style) ──────────────────────────────────────────
model_palette <- c("Logistic" = nature_grey, "ElasticNet" = nature_orange,
                   "RandomForest" = nature_blue, "SVM" = nature_red,
                   "NeuralNet" = nature_green)

Plot_model_ROC_nature <- function(model, dataset, color, add = FALSE, lty = 1) {
  pred <- predict(model, newdata = dataset, type = "prob")
  roc_obj <- roc(dataset$group, pred$Tumor, quiet = TRUE)
  if (add) {
    lines(1 - roc_obj$specificities, roc_obj$sensitivities,
          col = color, lwd = 2, lty = lty)
  } else {
    plot(roc_obj, col = color, lwd = 2, main = "ROC Curves — Classification Models")
  }
  return(as.numeric(auc(roc_obj)))
}

for (ds_label in c("Train", "Test")) {
  ds <- if (ds_label == "Train") train else test
  cairo_pdf(file.path(OUT_DIR, sprintf("Fig_14_roc_%s.pdf", tolower(ds_label))),
            width = 6, height = 6)
  aucs <- c(
    Plot_model_ROC_nature(logisreg, ds, model_palette["Logistic"]),
    Plot_model_ROC_nature(elastnet, ds, model_palette["ElasticNet"], add = TRUE, lty = 2),
    Plot_model_ROC_nature(RF,       ds, model_palette["RandomForest"], add = TRUE, lty = 3),
    Plot_model_ROC_nature(SVM,      ds, model_palette["SVM"], add = TRUE, lty = 4),
    Plot_model_ROC_nature(NNET,     ds, model_palette["NeuralNet"], add = TRUE, lty = 5))
  legend("bottomright",
         legend = paste(names(model_palette), "AUC=", round(aucs, 3)),
         col = model_palette, lwd = 2, lty = 1:5,
         bty = "n", cex = 0.8, title = paste(ds_label, "Data"))
  dev.off()
}
message("  ✓ ROC curves")

# ── 6.3 Model comparison (resampling profiles) ──────────────────────────────
p1 <- ggplot(elastnet) +
  geom_vline(xintercept = elastnet$bestTune$alpha, col = nature_grey, lty = 2) +
  theme_bw() + theme(legend.position = "top") +
  labs(title = "Elastic Net Regularized Regression",
       x = expression("Mixing Percentage (" * alpha * ")")) +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

p2 <- ggplot(RF) +
  geom_vline(xintercept = RF$bestTune$mtry, col = nature_grey, lty = 2) +
  theme_bw() + theme(legend.position = "top") +
  labs(title = "Random Forest", x = "Randomly Selected Predictors (mtry)") +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

p3 <- ggplot(SVM) +
  geom_vline(xintercept = SVM$bestTune$sigma, col = nature_grey, lty = 2) +
  theme_bw() + theme(legend.position = "top") +
  labs(title = "Support Vector Machine", x = "Sigma") +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

p4 <- ggplot(NNET) +
  geom_vline(xintercept = NNET$bestTune$size, col = nature_grey, lty = 2) +
  theme_bw() + theme(legend.position = "top") +
  labs(title = "Neural Network", x = "Hidden Units") +
  theme(plot.title = element_text(size = 11, hjust = 0.5))

cairo_pdf(file.path(OUT_DIR, "Fig_S6_model_tuning.pdf"), width = 10, height = 8)
gridExtra::grid.arrange(p1, p4, p3, p2, ncol = 2)
dev.off()
message("  ✓ Model tuning profiles")

# ── 6.4 Lift curves ─────────────────────────────────────────────────────────
PlotLift_nature <- function(dataset) {
  pred_list <- list(
    Logistic     = predict(logisreg, newdata = dataset, type = "prob")$Tumor,
    ElasticNet   = predict(elastnet, newdata = dataset, type = "prob")$Tumor,
    RandomForest = predict(RF, newdata = dataset, type = "prob")$Tumor,
    SVM          = predict(SVM, newdata = dataset, type = "prob")$Tumor,
    NeuralNet    = predict(NNET, newdata = dataset, type = "prob")$Tumor)

  actual <- ifelse(dataset$group == "Tumor", 1, 0)

  lift_data <- do.call(rbind, lapply(names(pred_list), function(m) {
    pred <- pred_list[[m]]
    idx  <- order(pred, decreasing = TRUE)
    cum_actual <- cumsum(actual[idx]) / sum(actual)
    data.frame(Model = m, Percentile = (1:length(idx)) / length(idx),
               Lift = cum_actual / ((1:length(idx)) / length(idx)))
  }))

  ggplot(lift_data, aes(x = Percentile, y = Lift, color = Model)) +
    geom_line(size = 1) +
    scale_color_manual(values = model_palette) +
    geom_hline(yintercept = 1, linetype = "dashed", color = nature_grey) +
    labs(title = "Cumulative Lift Curve", y = "Lift", x = "Percentile") +
    theme_pubr() + theme(legend.position = "bottom")
}

for (ds_label in c("Train", "Test")) {
  ds <- if (ds_label == "Train") train else test
  cairo_pdf(file.path(OUT_DIR,
    sprintf("Fig_15_lift_%s.pdf", tolower(ds_label))), width = 6, height = 5)
  print(PlotLift_nature(ds))
  dev.off()
}
message("  ✓ Lift curves")

# ── 6.5 Model comparison ────────────────────────────────────────────────────
model.list <- list(
  Logistic_Regression = logisreg,
  Random_Forest       = RF,
  Elastnet            = elastnet,
  SVM                 = SVM,
  Neural_Network      = NNET)
res <- resamples(model.list)

cairo_pdf(file.path(OUT_DIR, "Fig_16_model_comparison.pdf"), width = 12, height = 8)
par(mfrow = c(1, 4))
bwplot(res, layout = c(4, 1), main = "Model Performance Comparison")
dev.off()
message("  ✓ Model comparison plots")

###############################################################################
# PART 7: Summary & Export
###############################################################################

message("\n═══ PART 7: Summary & Export ═══")

# ── 7.1 Session info ────────────────────────────────────────────────────────
sink(file.path(OUT_DIR, "session_info.txt"))
cat("CRC lncRNA Pipeline — Nature-style\n")
cat("Completed:", as.character(Sys.time()), "\n\n")
sessionInfo()
sink()
message("  ✓ Session info saved")

# ── 7.2 Figure inventory ────────────────────────────────────────────────────
figs <- sort(list.files(OUT_DIR, pattern = "\\.pdf$"))
cat("\n═══ Generated Figures (", length(figs), "total) ═══\n", sep = "")
for (f in figs) cat("  ", f, "\n")

message("\n═══ Pipeline Complete ═══")
message("All outputs saved to: ", normalizePath(OUT_DIR))
