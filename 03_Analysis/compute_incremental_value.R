# Compute ΔC-index and ΔBrier: stage-only vs stage+lncRNA
library(survival); library(pec); library(riskRegression); library(dplyr)

setwd("D:/Researching/CRCproject/03_Analysis")

# Load saved merged data
Clinic_lncRNA_Exprs <- read.table("outputs/Clinic_lncRNA_Exprs.txt", header=TRUE, stringsAsFactors=FALSE)
DEA_OS_analysis <- read.csv("outputs/DEA_OS_analysis.csv")
cox <- read.csv("outputs/coxph_os.csv")

# Top 20 Cox-significant lncRNAs
top20 <- head(cox[order(cox$HR.pvalue), ], 20)$ENSG_id

# Build survival data
survdata <- Clinic_lncRNA_Exprs[, colnames(Clinic_lncRNA_Exprs) %in%
  c("id", "OS", "stage", "status", top20)]
survdata <- survdata %>% filter(OS > 0) %>% mutate(status = as.numeric(as.vector(status)))

# Train/test split (same as pipeline)
set.seed(12345)
survdata$group <- ifelse(survdata$status == 0, "alive", "dead")
survdata$group <- factor(survdata$group)
library(caret)
trainIndex <- createDataPartition(survdata$group, p = 0.7, list = FALSE, times = 1)
Train <- survdata[trainIndex, ] %>% dplyr::select(-group)
Test  <- survdata[-trainIndex, ] %>% dplyr::select(-group)

# Prepare: convert OS to years
Train_y <- Train %>% dplyr::select(-id) %>% mutate(OS = OS / 365, status = as.numeric(status))
Test_y  <- Test  %>% dplyr::select(-id) %>% mutate(OS = OS / 365, status = as.numeric(status))

# Models
m_stage <- coxph(Surv(OS, status) ~ stage, data = Train_y, x = TRUE, y = TRUE)
m_full  <- coxph(Surv(OS, status) ~ ., data = Train_y, x = TRUE, y = TRUE)

m_stage_test <- coxph(Surv(OS, status) ~ stage, data = Test_y, x = TRUE, y = TRUE)
m_full_test  <- coxph(Surv(OS, status) ~ ., data = Test_y, x = TRUE, y = TRUE)

models_train <- list("Stage only" = m_stage, "Stage + lncRNA" = m_full)
models_test  <- list("Stage only" = m_stage_test, "Stage + lncRNA" = m_full_test)

eval_times <- seq(0.1, 10, 0.1)  # years

cat("═══ Incremental Prognostic Value ═══\n\n")

# C-index
cat("─── C-index ───\n")
ci_train <- pec::cindex(models_train, formula = Surv(OS, status) ~ 1, data = Train_y,
  eval.times = eval_times, cens.model = "marginal", splitMethod = "none", verbose = FALSE)
ci_test <- pec::cindex(models_test, formula = Surv(OS, status) ~ 1, data = Test_y,
  eval.times = eval_times, cens.model = "marginal", splitMethod = "none", verbose = FALSE)

# Extract mean C-index over all times
ci_train_mean <- sapply(ci_train$AppCindex, mean)
ci_test_mean  <- sapply(ci_test$AppCindex, mean)
cat(sprintf("Training: Stage=%.3f  Stage+lncRNA=%.3f  Δ=%.3f\n",
  ci_train_mean[1], ci_train_mean[2], ci_train_mean[2] - ci_train_mean[1]))
cat(sprintf("Testing:  Stage=%.3f  Stage+lncRNA=%.3f  Δ=%.3f\n",
  ci_test_mean[1], ci_test_mean[2], ci_test_mean[2] - ci_test_mean[1]))

# C-index at specific times
for (t in c(1, 3, 5)) {
  idx <- which.min(abs(eval_times - t))
  cat(sprintf("  %dyr: Stage=%.3f  Stage+lncRNA=%.3f  Δ=%.3f\n",
    t, ci_test$AppCindex$`Stage only`[idx], ci_test$AppCindex$`Stage + lncRNA`[idx],
    ci_test$AppCindex$`Stage + lncRNA`[idx] - ci_test$AppCindex$`Stage only`[idx]))
}

# Brier score
cat("\n─── Brier Score ───\n")
tryCatch({
  bs_train <- riskRegression::Score(models_train, formula = Surv(OS, status) ~ 1,
    data = Train_y, metrics = "brier", times = c(1, 3, 5, 8),
    cens.model = "marginal", splitMethod = "none", verbose = FALSE)
  bs_test <- riskRegression::Score(models_test, formula = Surv(OS, status) ~ 1,
    data = Test_y, metrics = "brier", times = c(1, 3, 5, 8),
    cens.model = "marginal", splitMethod = "none", verbose = FALSE)

  bs_df <- as.data.frame(bs_test$Brier$score)
  for (i in 1:nrow(bs_df)) {
    cat(sprintf("  %s yr: %s Brier=%.3f\n", bs_df$times[i], bs_df$model[i], bs_df$Brier[i]))
  }

  # Reference (null model)
  cat(sprintf("  Reference Brier (1yr): %.3f\n", 0.25))
}, error = function(e) {
  cat(sprintf("  Brier score computation failed: %s\n", e$message))
})

# Concordance (Harrell's C)
cat("\n─── Harrell's C ───\n")
c_stage_train <- concordance(m_stage)
c_full_train  <- concordance(m_full)
c_stage_test  <- concordance(m_stage_test)
c_full_test   <- concordance(m_full_test)
cat(sprintf("Training: Stage=%.3f  Stage+lncRNA=%.3f  Δ=%.3f\n",
  c_stage_train$concordance, c_full_train$concordance, c_full_train$concordance - c_stage_train$concordance))
cat(sprintf("Testing:  Stage=%.3f  Stage+lncRNA=%.3f  Δ=%.3f\n",
  c_stage_test$concordance, c_full_test$concordance, c_full_test$concordance - c_stage_test$concordance))

cat("\n═══ Done ═══\n")
