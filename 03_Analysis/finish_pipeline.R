###############################################################################
# finish_pipeline.R — Generate missing Cox prediction & evaluation figures
# Run in RStudio: setwd("03_Analysis"); source("finish_pipeline.R")
###############################################################################

# ── Setup ──
library(survival); library(survminer); library(dplyr); library(ggplot2)
library(ggpubr); library(ggsci); library(survivalROC); library(caret)
library(doParallel); library(parallel)
library(pec); library(riskRegression); library(prodlim)

OUT_DIR <- "outputs"
DATA_DIR_RAW <- "../02_Data/raw"

# Nature palette
nature_blue  <- "#1F4E79"; nature_red   <- "#C0392B"
nature_green <- "#27AE60"; nature_grey  <- "#95A5A6"
nature_orange <- "#E67E22"; nature_purple <- "#8E44AD"
pal_train_test <- c(nature_blue, nature_red)
pal_risk       <- c(nature_blue, nature_red)
pal_status     <- c(nature_blue, nature_red)
model_colors <- c("CoxPH (stage)" = nature_grey,
                  "CoxPH (+lncRNA)" = nature_blue,
                  "Random Forest" = nature_red)

# Helper
lookup <- function(terms, key.match, missing = NA) {
  if (is.data.frame(key.match) && ncol(key.match) >= 2) {
    idx <- match(terms, key.match[[1]])
    out <- key.match[[2]][idx]
    out[is.na(idx)] <- missing
    return(out)
  }
  return(rep(missing, length(terms)))
}

# Load saved data
Clinic_lncRNA_Exprs <- read.table(file.path(OUT_DIR, "Clinic_lncRNA_Exprs.txt"),
                                   header=TRUE, stringsAsFactors=FALSE)
DEA_OS_analysis <- read.csv(file.path(OUT_DIR, "DEA_OS_analysis.csv"), stringsAsFactors=FALSE)
coxph_os_res <- read.csv(file.path(OUT_DIR, "coxph_os.csv"), stringsAsFactors=FALSE)

message(sprintf("Loaded: %d patients, %d DEA-OS lncRNAs, %d Cox lncRNAs",
        nrow(Clinic_lncRNA_Exprs), nrow(DEA_OS_analysis), nrow(coxph_os_res)))

# Use top 20 most significant Cox lncRNAs (avoid stepwise timeout)
top_genes <- coxph_os_res %>% arrange(HR.pvalue) %>% head(20)
top_ensg <- top_genes$ENSG_id

# Build survdata
survdata <- Clinic_lncRNA_Exprs[, colnames(Clinic_lncRNA_Exprs) %in%
  c("id", "OS", "stage", "status", top_ensg)]
survdata <- survdata %>%
  filter(OS > 0) %>%
  mutate(status = as.numeric(as.vector(status)))
message(sprintf("Survival matrix: %d × %d", nrow(survdata), ncol(survdata)))

# Train/Test split
survdata$group <- ifelse(survdata$status == 0, "alive", "dead")
survdata$group <- factor(survdata$group)
set.seed(12345)
trainIndex <- createDataPartition(survdata$group, p = 0.7, list = FALSE, times = 1)
Train <- survdata[trainIndex, ] %>% dplyr::select(-group)
Test  <- survdata[-trainIndex, ] %>% dplyr::select(-group)
message(sprintf("Train=%d  Test=%d", nrow(Train), nrow(Test)))

# ── Multivariate Cox (no stepwise — use top 20) ──
cox_data_train <- Train %>% dplyr::select(-id) %>% mutate(status = as.numeric(status))
cox_data_test  <- Test  %>% dplyr::select(-id) %>% mutate(status = as.numeric(status))

cox_fit <- coxph(Surv(OS, status) ~ ., data = cox_data_train)
cox_vip_list <- names(coef(cox_fit))
message(sprintf("Cox model: %d features, AIC=%.2f", length(cox_vip_list), AIC(cox_fit)))

# ── Forest plot ──
cairo_pdf(file.path(OUT_DIR, "Fig_3_forestplot_cox.pdf"), width = 7, height = 5)
print(ggforest(cox_fit, main = "Feature importance in OS (TCGA-COAD)",
               data = cox_data_train))
dev.off()
message("✓ Forest plot")

# ── Time-dependent ROC (1/3/5 year) ──
PlotsurvROC <- function(predict.time, Train, Test, model = cox_fit, vip_list) {
  # Training
  newdata_t <- Train[, colnames(Train) %in% c(vip_list, "status", "stage", "OS", "id")]
  riskScore_t <- predict(model, type = "risk", newdata = newdata_t)
  risk_df_t <- cbind.data.frame(newdata_t[, c("id", "OS", "status")],
                                 riskScore = riskScore_t)
  roc_t <- survivalROC(Stime = risk_df_t$OS / 365, status = risk_df_t$status,
                        marker = risk_df_t$riskScore, predict.time = predict.time,
                        method = "KM")
  # Testing
  newdata_v <- Test[, colnames(Test) %in% c(vip_list, "status", "stage", "OS", "id")]
  riskScore_v <- predict(model, type = "risk", newdata = newdata_v)
  risk_df_v <- cbind.data.frame(newdata_v[, c("id", "OS", "status")],
                                 riskScore = riskScore_v)
  roc_v <- survivalROC(Stime = risk_df_v$OS / 365, status = risk_df_v$status,
                        marker = risk_df_v$riskScore, predict.time = predict.time,
                        method = "KM")

  plot(roc_t$FP, roc_t$TP, type = "l", xlim = c(0, 1), ylim = c(0, 1),
       col = nature_blue, lwd = 2, xlab = "False positive rate",
       ylab = "True positive rate",
       main = paste("ROC of CoxPH Risk Prediction (Year =", predict.time, ")"))
  abline(0, 1, col = nature_grey, lty = 2)
  lines(roc_v$FP, roc_v$TP, type = "l", col = nature_red, lwd = 2)
  legend("bottomright", col = c(nature_blue, nature_red), bty = "n",
         cex = 0.8, lwd = 2, title = "ROC AUC",
         legend = c(paste("Training:", round(roc_t$AUC, 3)),
                    paste("Testing:",  round(roc_v$AUC, 3))))
}

for (yr in c(1, 3, 5)) {
  cairo_pdf(file.path(OUT_DIR, sprintf("Fig_4_roc_%dyear.pdf", yr)), width = 5, height = 5)
  PlotsurvROC(yr, Train, Test, cox_fit, cox_vip_list)
  dev.off()
}
message("✓ ROC curves (1/3/5 year)")

# ── AUC over time ──
risk_train <- {
  newdata <- Train[, colnames(Train) %in% c(cox_vip_list, "status", "stage", "OS", "id")]
  cbind.data.frame(newdata[, c("id", "OS", "status")],
                   riskScore = predict(cox_fit, type = "risk", newdata = newdata))
}
risk_test <- {
  newdata <- Test[, colnames(Test) %in% c(cox_vip_list, "status", "stage", "OS", "id")]
  cbind.data.frame(newdata[, c("id", "OS", "status")],
                   riskScore = predict(cox_fit, type = "risk", newdata = newdata))
}

CALAUC <- function(time, risk_df) {
  roc <- survivalROC(Stime = risk_df$OS / 365, status = risk_df$status,
                      marker = risk_df$riskScore, predict.time = time, method = "KM")
  roc$AUC
}

time_grid <- seq(1, 15, 0.05)
AUC <- data.frame(
  time = rep(time_grid, 2),
  group = rep(c("Training", "Testing"), each = length(time_grid)),
  AUC = c(sapply(time_grid, CALAUC, risk_df = risk_train),
          sapply(time_grid, CALAUC, risk_df = risk_test)))

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
message("✓ AUC over time")

# ── Risk-stratified KM curves ──
GGsurvplot_nature <- function(lncrna_list, dataset, cox_model) {
  newdata <- dataset[, colnames(dataset) %in% c(lncrna_list, "status", "stage", "OS", "id")]
  riskScore <- predict(cox_model, type = "risk", newdata = newdata)
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
  print(GGsurvplot_nature(cox_vip_list, ds, cox_fit))
  dev.off()
}
message("✓ Risk-stratified KM curves")

# ── Risk score distributions ──
risk_fun <- function(dataset) {
  newdata <- dataset[, colnames(dataset) %in% c(cox_vip_list, "status", "stage", "OS", "id")]
  riskScore <- predict(cox_fit, type = "risk", newdata = newdata)
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
message("✓ Risk score distribution plots")

# ── Brier Score & C-index ──
rfdata_train <- cox_data_train %>% mutate(OS = OS / 365)
rfdata_test  <- cox_data_test  %>% mutate(OS = OS / 365)

# Collapse stages to avoid NA in small test sets
rfdata_train$stage_simple <- ifelse(rfdata_train$stage %in% c("i","ii"), "early", "late")
rfdata_test$stage_simple  <- ifelse(rfdata_test$stage %in% c("i","ii"), "early", "late")

# Use simple Cox models that won't have NA coefficients
simple_train <- rfdata_train[, c("OS","status","stage_simple")]
simple_test  <- rfdata_test[, c("OS","status","stage_simple")]

# Only keep lncRNA columns that exist
lnc_cols <- intersect(cox_vip_list, colnames(rfdata_train))
rf_train_cox <- cbind(simple_train, rfdata_train[, lnc_cols])
rf_test_cox  <- cbind(simple_test,  rfdata_test[, lnc_cols])

# Simple models
Models_train <- list(
  "CoxPH (stage)" = coxph(Surv(OS, status) ~ stage_simple, data = rf_train_cox, x = TRUE, y = TRUE),
  "CoxPH (+lncRNA)" = coxph(Surv(OS, status) ~ ., data = rf_train_cox, x = TRUE, y = TRUE))

Models_test <- list(
  "CoxPH (stage)" = coxph(Surv(OS, status) ~ stage_simple, data = rf_test_cox, x = TRUE, y = TRUE),
  "CoxPH (+lncRNA)" = coxph(Surv(OS, status) ~ ., data = rf_test_cox, x = TRUE, y = TRUE))

eval_days <- seq(30, 3650, 30)

# C-index
message("Computing C-index...")
tryCatch({
  ci_train <- pec::cindex(Models_train,
    formula = Surv(OS, status) ~ 1, data = rf_train_cox,
    eval.times = eval_days / 365, cens.model = "marginal",
    splitMethod = "bootcv", B = 50,
    M = round(nrow(rf_train_cox) * 0.6),
    verbose = FALSE, maxtime = max(eval_days / 365))

  ci_test <- pec::cindex(Models_test,
    formula = Surv(OS, status) ~ 1, data = rf_test_cox,
    eval.times = eval_days / 365, cens.model = "marginal",
    splitMethod = "bootcv", B = 50,
    M = round(nrow(rf_test_cox) * 0.6),
    verbose = FALSE, maxtime = max(eval_days / 365))

  model_colors_ci <- c("CoxPH (stage)" = nature_grey, "CoxPH (+lncRNA)" = nature_blue)
  for (label in c("train", "test")) {
    obj <- if (label == "train") ci_train else ci_test
    cairo_pdf(file.path(OUT_DIR, sprintf("Fig_13_cindex_%s.pdf", label)),
              width = 6, height = 5)
    plot(obj, smooth = TRUE, lwd = 1.5, legend.cex = 0.8,
         type = "s", add.refline = TRUE,
         xlim = c(0, 3650), ylim = c(0, 1),
         xlab = "Time (days)", ylab = "Concordance Index (C-index)",
         col = model_colors_ci)
    dev.off()
  }
  message("✓ C-index plots")
}, error = function(e) {
  message(sprintf("⚠ C-index skipped: %s", e$message))
})

# Brier Score
message("Computing Brier scores...")
tryCatch({
  score_train <- riskRegression::Score(Models_train,
    formula = Surv(OS, status) ~ 1, data = rf_train_cox,
    metrics = "brier", times = eval_days / 365,
    cens.model = "marginal", splitMethod = "bootcv",
    B = 50, M = round(nrow(rf_train_cox) * 0.6), verbose = FALSE)

  score_test <- riskRegression::Score(Models_test,
    formula = Surv(OS, status) ~ 1, data = rf_test_cox,
    metrics = "brier", times = eval_days / 365,
    cens.model = "marginal", splitMethod = "bootcv",
    B = 50, M = round(nrow(rf_test_cox) * 0.6), verbose = FALSE)

  model_colors_brier <- c("CoxPH (stage)" = nature_grey, "CoxPH (+lncRNA)" = nature_blue)
  for (label in c("train", "test")) {
    sc <- if (label == "train") score_train else score_test
    bd <- as.data.frame(sc$Brier$score)
    bd$times_days <- bd$times * 365
    cairo_pdf(file.path(OUT_DIR, sprintf("Fig_12_brier_%s.pdf", label)),
              width = 6, height = 5)
    p <- ggplot(bd, aes(x = times_days, y = Brier, color = model, fill = model)) +
      geom_line(size = 1.1) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.08, color = NA) +
      scale_color_manual(values = model_colors_brier) +
      scale_fill_manual(values = model_colors_brier) +
      labs(title = paste("Prediction Error (Brier Score) —", label),
           y = "Brier Score", x = "Time (days)") +
      ylim(0, 0.5) + xlim(0, 3650) +
      theme_pubr() + theme(legend.position = "bottom",
                            plot.title = element_text(size = 13, hjust = 0.5))
    print(p)
    dev.off()
  }
  message("✓ Brier score plots")
}, error = function(e) {
  message(sprintf("⚠ Brier score skipped: %s", e$message))
})

# ── Session info ──
sink(file.path(OUT_DIR, "session_info.txt"))
cat("CRC lncRNA Pipeline — finish_pipeline.R\n")
cat("Completed:", as.character(Sys.time()), "\n\n")
sessionInfo()
sink()

# Final inventory
figs <- sort(list.files(OUT_DIR, pattern = "\\.pdf$"))
message(sprintf("\n═══ All Figures (%d total) ═══", length(figs)))
for (f in figs) message("  ", f)
message("\n═══ finish_pipeline.R complete ═══")
