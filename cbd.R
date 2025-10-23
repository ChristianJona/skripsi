# 2.2. Fungsi CBD (Model & Prediksi) - Sama
fit_cbd_model <- function(ages_matrix, qx_matrix, forecast_periods = 7, verbose = TRUE) {
  # ... (kode fit_cbd_model asli) ...
  if(verbose) cat("Fitting CBD model...\n"); hx_matrix <- -log(1-pmax(qx_matrix, 1e-10)); hx_matrix[is.infinite(hx_matrix)] <- NA; n_ages <- nrow(hx_matrix); n_periods <- ncol(hx_matrix); all_ages <- as.vector(ages_matrix); mean_age <- mean(all_ages, na.rm=TRUE); if(verbose) cat("  Mean age:", round(mean_age,2), "\n"); kappa0 <- numeric(n_periods); kappa1 <- numeric(n_periods)
  for(t in 1:n_periods) { valid_idx <- !is.na(ages_matrix[,t]) & !is.na(hx_matrix[,t]); if(sum(valid_idx)<2){kappa0[t]<-NA; kappa1[t]<-NA; next}; ages_t<-ages_matrix[valid_idx,t]; log_hx_t<-log(hx_matrix[valid_idx,t]); centered_ages<-ages_t-mean_age; tryCatch({ lm_fit<-lm(log_hx_t~centered_ages); kappa0[t]<-coef(lm_fit)[1]; kappa1[t]<-coef(lm_fit)[2]}, error=function(e){kappa0[t]<-NA; kappa1[t]<-NA})}
  time_index<-1:n_periods; valid_kappa<-!is.na(kappa0)&!is.na(kappa1); kappa0_forecast<-rep(NA, forecast_periods); kappa1_forecast<-rep(NA, forecast_periods)
  if(sum(valid_kappa)>=3){ kappa0_trend<-lm(kappa0[valid_kappa]~time_index[valid_kappa]); kappa1_trend<-lm(kappa1[valid_kappa]~time_index[valid_kappa]); future_time<-(n_periods+1):(n_periods+forecast_periods); kappa0_forecast<-predict(kappa0_trend, newdata=data.frame(time_index=future_time)); kappa1_forecast<-predict(kappa1_trend, newdata=data.frame(time_index=future_time)); if(verbose) cat("  Forecast kappa OK.\n")} else { if(verbose) cat("  Warning: Kappa < 3, cannot forecast.\n")}
  cbd_model<-list(kappa0=kappa0, kappa1=kappa1, kappa0_forecast=kappa0_forecast, kappa1_forecast=kappa1_forecast, mean_age=mean_age, ages_matrix=ages_matrix, qx_matrix=qx_matrix, n_periods=n_periods, forecast_periods=forecast_periods); class(cbd_model)<-"cbd"; return(cbd_model)
}
predict.cbd <- function(object, newdata_ages, period = NULL, forecast_period = NULL, ...) {
  # ... (kode predict.cbd asli dengan as.numeric) ...
  if(is.null(period)==is.null(forecast_period)) stop("Specify either 'period' or 'forecast_period'."); k0<-NA; k1<-NA
  if(!is.null(period)){if(period<1||period>object$n_periods) stop("Invalid period."); k0<-object$kappa0[period]; k1<-object$kappa1[period]} else {if(forecast_period<1||forecast_period>object$forecast_periods) stop("Invalid forecast_period."); k0<-object$kappa0_forecast[forecast_period]; k1<-object$kappa1_forecast[forecast_period]}
  if(is.na(k0)||is.na(k1)) return(rep(NA,length(newdata_ages))); centered_ages<-as.numeric(newdata_ages)-object$mean_age; log_hx_pred<-k0+k1*centered_ages; hx_pred<-exp(log_hx_pred); qx_pred<-1-exp(-hx_pred); qx_pred<-pmax(0, pmin(1, qx_pred, na.rm=T), na.rm=T); return(qx_pred)
}