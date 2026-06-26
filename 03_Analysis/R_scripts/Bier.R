Bier = function(dataset, cl = NULL) {
	require(prodlim)
	require(pec)
	require(riskRegression)
	require(doParallel)

	pecdata = dataset[, colnames(dataset) %in% c(cox_vip_list, "stage", "OS", "status")]
	Models <- list(
		"CoxPH.stage" = coxph(Surv(OS, status) ~ stage, data = pecdata, x = T, y = T),
		"CoxPH" = coxph(Surv(OS, status) ~ ., data = pecdata, x = T, y = T),
		"RandomForest" = rfsrc)

	if (!is.null(cl)) {
		registerDoParallel(cl)
		on.exit(registerDoSEQ())
	}

	Bier.imp <- pec::pec(Models, formula = Hist(OS, status) ~ ., data = pecdata,
		cens.model = "marginal", splitMethod = "bootcv", M = round(nrow(pecdata)*0.6), B = 500,
		keep.index = T, multiSplitTest = T, confInt = T, confLevel = 0.95,
		exact = T, verbose = T, maxtime = 3000, eval.times = seq(0, 3650, 30))
	return(Bier.imp)
}
