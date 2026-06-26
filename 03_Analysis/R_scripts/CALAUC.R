CALAUC_precompute = function(dataset, cox_vip_list, stepcox) {
	require(survivalROC)
	newdata = dataset[, colnames(dataset) %in% c(cox_vip_list, "status", "stage", "OS", "id")]
	riskScore = predict(stepcox, type = "risk", newdata = newdata)
	risk_df = cbind.data.frame(newdata[, c("id", "OS", "status")], riskScore = riskScore)
	list(risk_df = risk_df)
}

CALAUC_train = function(time) {
	auc = survivalROC(Stime = pre_train$risk_df$OS / 365,
		status = pre_train$risk_df$status,
		marker = pre_train$risk_df$riskScore,
		predict.time = time, method = "KM")$AUC
	auc
}

CALAUC_test = function(time) {
	auc = survivalROC(Stime = pre_test$risk_df$OS / 365,
		status = pre_test$risk_df$status,
		marker = pre_test$risk_df$riskScore,
		predict.time = time, method = "KM")$AUC
	auc
}
