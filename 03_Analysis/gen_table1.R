dea_os <- read.csv("outputs/DEA_OS_analysis.csv")
cox <- read.csv("outputs/coxph_os.csv")
logr <- read.table("outputs/log_rank_os.txt", header = TRUE)

# Top 17 by Cox significance
idx <- order(cox$HR.pvalue)[1:17]
top <- cox[idx, ]

# Match DEA and logrank stats by ENSG
m_dea <- match(top$ENSG_id, dea_os$ENSG_id)
m_lr  <- match(top$ENSG_id, logr$ENSG_id)

cat("% Auto-generated Table 1 rows - top 17 by Cox significance\n")
for (i in 1:nrow(top)) {
  fc <- dea_os$logFC[m_dea[i]]
  fdr <- dea_os$adj.P.Val[m_dea[i]]
  hr <- top$HR[i]
  pv <- top$HR.pvalue[i]
  lr <- logr$logrank.pvalue[m_lr[i]]

  fc_str <- if (is.na(fc)) "---" else sprintf("%.3f", fc)
  fdr_str <- if (is.na(fdr)) "---" else sprintf("%.2e", fdr)
  hr_str <- sprintf("%.2f", hr)
  pv_str <- sprintf("%.2e", pv)
  lr_str <- sprintf("%.2e", lr)

  prefix <- if (i %% 2 == 0) "\\rowcolor{LightBlue} " else ""

  cat(sprintf("%s%s & %s & %s & %s & %s & %s & %s & %s \\\\\n",
    prefix, top$Symbol[i], top$ENSG_id[i], fc_str, fdr_str,
    hr_str, top$CI[i], pv_str, lr_str))
}
