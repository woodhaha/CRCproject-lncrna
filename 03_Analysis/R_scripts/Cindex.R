Cindex = function(dataset, cl = NULL) {
	require(prodlim)
	require(pec)
	require(riskRegression)
	require(doParallel)

	cdata = dataset[, colnames(dataset) %in% c(cox_vip_list, "stage", "OS", "status")]
	Models <- list(
		"CoxPH.stage" = coxph(Surv(OS, status) ~ stage, data = cdata, x = T, y = T),
		"CoxPH" = coxph(Surv(OS, status) ~ ., data = cdata, x = T, y = T),
		"RandomForest" = rfsrc)

	if (!is.null(cl)) {
		registerDoParallel(cl)
		on.exit(registerDoSEQ())
	}

	Cimp = pec::cindex(Models, formula = Hist(OS, status) ~ ., data = cdata,
		eval.times = seq(0, 3650, 30), splitMethod = "bootcv",
		M = round(nrow(cdata)*0.6), cens.model = "marginal", B = 500, maxtime = 3000)
	return(Cimp)
}
