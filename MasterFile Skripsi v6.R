library(tidyverse)
library(knitr)
library(MortalityLaws)
library(scales)
library(xtable)

file_sources <- list.files("R", pattern="*.R$", full.names=TRUE)
sapply(file_sources, source)
cat("Fungsi model dimuat\n")

data_bersih <- readRDS("C:/Users/rapha/Documents/1 PROJECT/0 SKRIPSI/Raw Data dan R/data bersih/data_final_siap_pakai.rds")
data_females <- data_bersih %>% filter(Sex == "Female")
data_males   <- data_bersih %>% filter(Sex == "Male")

cat("Data dimuat - Kohor wanita:", unique(data_females$Cohort), "\n")

# FUNGSI COALE-KISKER
fit_coale_kisker <- function(ages, qx, mx = NULL, x0 = 84, x1 = 110, mx1 = 1.0) {
  # PENTING: pred_ages harus lebih luas dari ages untuk extrapolation
  pred_ages <- min(ages):x1  # Extrapolate sampai x1
  
  if(is.null(mx)) mx <- qx / (1 - 0.5 * qx)
  
  # Hitung mx_x0 (starting point extrapolation)
  if(x0 %in% ages) {
    mx_x0 <- mx[ages == x0]
  } else {
    mx_x0 <- approx(ages, mx, xout = x0, rule = 2)$y
  }
  
  # Hitung parameter R
  if(x0 %in% ages && (x0+1) %in% ages) {
    k_x0 <- log(mx[ages == x0] / mx[ages == (x0+1)])
  } else {
    mx_x0_plus1 <- approx(ages, mx, xout = x0+1, rule = 2)$y
    k_x0 <- log(mx_x0_plus1/ mx_x0)
  }
  
  sum_i <- sum(1:(x1 - x0))
  R <- (k_x0 * ((x1 - x0 + 1)/2) + log(mx_x0) - log(mx1)) / sum_i
  
  
  # Prediksi mx untuk semua ages
  pred_mx <- numeric(length(pred_ages))
  
  for(i in seq_along(pred_ages)) {
    age <- pred_ages[i]
    
    if(age < x0) {
      # Gunakan observed atau interpolasi
      if(age %in% ages) {
        pred_mx[i] <- mx[ages == age]
      } else {
        pred_mx[i] <- approx(ages, mx, xout = age, rule = 2)$y
      }
    } else {
      # Coale-Kisker extrapolation
      k_age <- k_x0 - R * (age - x0)
      
      if(age == x0) {
        pred_mx[i] <- mx_x0
      } else {
        # Recursive calculation
        pred_mx[i] <- pred_mx[i-1] * exp(k_age)
      }
    }
  }
  
  # Convert mx to qx
  pred_qx <- pred_mx / (1 + 0.5 * pred_mx)
  pred_qx <- pmin(pred_qx, 1)
  
  return(pred_qx)
}

# FUNGSI CBD
fit_cbd_model <- function(ages_matrix, qx_matrix, forecast_periods = 7, verbose = TRUE) {
  if(verbose) cat("Fitting CBD model...\n")
  hx_matrix <- -log(1 - pmax(qx_matrix, 1e-10))
  n_ages <- nrow(hx_matrix)
  n_periods <- ncol(hx_matrix)
  all_ages <- as.vector(ages_matrix)
  mean_age <- mean(all_ages, na.rm = TRUE)
  if(verbose) cat("  Mean age:", round(mean_age, 2), "\n")
  
  kappa0 <- numeric(n_periods)
  kappa1 <- numeric(n_periods)
  
  for(t in 1:n_periods) {
    valid_idx <- !is.na(ages_matrix[, t]) & !is.na(hx_matrix[, t])
    if(sum(valid_idx) < 2) {
      kappa0[t] <- NA
      kappa1[t] <- NA
      next
    }
    ages_t <- ages_matrix[valid_idx, t]
    log_hx_t <- log(hx_matrix[valid_idx, t])
    centered_ages <- ages_t - mean_age
    
    tryCatch({
      lm_fit <- lm(log_hx_t ~ centered_ages)
      kappa0[t] <- coef(lm_fit)[1]
      kappa1[t] <- coef(lm_fit)[2]
    }, error = function(e) {
      kappa0[t] <- NA
      kappa1[t] <- NA
    })
  }
  
  time_index <- 1:n_periods
  valid_kappa <- !is.na(kappa0) & !is.na(kappa1)
  
  if(sum(valid_kappa) < 3) {
    kappa0_forecast <- rep(NA, forecast_periods)
    kappa1_forecast <- rep(NA, forecast_periods)
  } else {
    kappa0_trend <- lm(kappa0[valid_kappa] ~ time_index[valid_kappa])
    kappa1_trend <- lm(kappa1[valid_kappa] ~ time_index[valid_kappa])
    future_time <- (n_periods + 1):(n_periods + forecast_periods)
    kappa0_forecast <- predict(kappa0_trend, newdata = data.frame(time_index = future_time))
    kappa1_forecast <- predict(kappa1_trend, newdata = data.frame(time_index = future_time))
  }
  
  cbd_model <- list(
    kappa0 = kappa0, kappa1 = kappa1,
    kappa0_forecast = kappa0_forecast, kappa1_forecast = kappa1_forecast,
    mean_age = mean_age, ages_matrix = ages_matrix, qx_matrix = qx_matrix,
    n_periods = n_periods, forecast_periods = forecast_periods
  )
  class(cbd_model) <- "cbd"
  return(cbd_model)
}

predict.cbd <- function(object, newdata, period = NULL, forecast_period = NULL, ...) {
  if(is.null(period) && is.null(forecast_period)) {
    stop("Must specify either 'period' or 'forecast_period'")
  }
  
  if(!is.null(period)) {
    kappa0 <- object$kappa0[period]
    kappa1 <- object$kappa1[period]
  } else {
    kappa0 <- object$kappa0_forecast[forecast_period]
    kappa1 <- object$kappa1_forecast[forecast_period]
  }
  
  centered_ages <- newdata - object$mean_age
  log_hx_pred <- kappa0 + kappa1 * centered_ages
  hx_pred <- exp(log_hx_pred)
  qx_pred <- 1 - exp(-hx_pred)
  qx_pred <- pmax(0, pmin(1, qx_pred))
  return(qx_pred)
}

# FITTING SEMUA MODEL STATIS
fit_static_models_updated <- function(cohort_data, verbose = TRUE) {
  if(verbose) cat("Fitting models untuk kohor", unique(cohort_data$Cohort), "...\n")
  models <- list()
  
  tryCatch({
    models$stlt <- stlt(ages = cohort_data$Age, qx = cohort_data$qx)
    if(verbose) cat("  STLT: OK\n")
  }, error = function(e) {
    if(verbose) cat("  STLT: GAGAL\n")
    models$stlt <- NULL
  })
  
  tryCatch({
    models$tlt <- tlt(ages = cohort_data$Age, qx = cohort_data$qx, 
                      lx = cohort_data$lx, dx = cohort_data$dx)
    if(verbose) cat("  TLT: OK\n")
  }, error = function(e) {
    if(verbose) cat("  TLT: GAGAL\n")
    models$tlt <- NULL
  })
  
  target_models <- c("gompertz", "makeham", "HP2")
  for(law in target_models) {
    tryCatch({
      models[[law]] <- get_qx(x = cohort_data$Age, qx = cohort_data$qx, 
                              law = law, pred_ages = cohort_data$Age)
      if(verbose) cat("  ", law, ": OK\n")
    }, error = function(e) {
      if(verbose) cat("  ", law, ": GAGAL\n")
      models[[law]] <- NULL
    })
  }
  
  # --- COALE-KISKER ---
  tryCatch({
    ck_full <- fit_coale_kisker(
      ages = cohort_data$Age, 
      qx = cohort_data$qx,
      x0 = 95,   # Start extrapolation age
      x1 = 110   # End extrapolation age
    )
    
    # Subset to match cohort_data length
    models$coale_kisker <- ck_full[1:nrow(cohort_data)]
    
    if(verbose) cat("  Coale-Kisker: OK\n")
    
  }, error = function(e) {
    if(verbose) cat("  Coale-Kisker: GAGAL -", e$message, "\n")
    models$coale_kisker <- NULL
  })
  
  return(models)
}

# TABEL PERBANDINGAN
generate_comparison_table_fixed <- function(cohort_data, models, output_format = "kable") {
  
  pred_list <- list()
  
  # --- COLLECT PREDICTIONS FROM ALL MODELS ---
  
  # STLT
  if(!is.null(models$stlt)) {
    tryCatch({
      pred_list$STLT <- predict(models$stlt, newdata = cohort_data$Age)
    }, error = function(e) {
      warning("STLT prediction failed: ", e$message)
    })
  }
  
  # TLT
  if(!is.null(models$tlt)) {
    tryCatch({
      pred_list$TLT <- predict.tlt(models$tlt, newdata = cohort_data$Age)
    }, error = function(e) {
      warning("TLT prediction failed: ", e$message)
    })
  }
  
  # Gompertz, Makeham, HP2
  mortality_laws <- c("gompertz", "makeham", "HP2")
  law_names <- c("GOMPERTZ", "MAKEHAM", "HP2")
  
  for(j in seq_along(mortality_laws)) {
    law <- mortality_laws[j]
    if(!is.null(models[[law]])) {
      pred_list[[law_names[j]]] <- models[[law]]
    }
  }
  
  # Coale-Kisker
  if(!is.null(models$coale_kisker)) {
    tryCatch({
      # Coale-Kisker may return extended ages, subset to match cohort_data
      ck_pred <- models$coale_kisker
      
      # Ensure length matches
      if(length(ck_pred) >= nrow(cohort_data)) {
        pred_list$`COALE-KISKER` <- ck_pred[1:nrow(cohort_data)]
      } else {
        warning("Coale-Kisker prediction too short")
      }
      
    }, error = function(e) {
      warning("Coale-Kisker prediction failed: ", e$message)
    })
  }
  
  # --- CHECK IF ANY PREDICTIONS AVAILABLE ---
  if(length(pred_list) == 0) {
    warning("No model predictions available for comparison")
    return(NULL)
  }
  
  # --- VALIDATE LENGTHS ---
  pred_lengths <- sapply(pred_list, length)
  expected_length <- nrow(cohort_data)
  
  # Remove models with incorrect length
  valid_models <- names(pred_list)[pred_lengths == expected_length]
  
  if(length(valid_models) == 0) {
    warning("No models have correct prediction length")
    return(NULL)
  }
  
  pred_list <- pred_list[valid_models]
  
  # --- CREATE PREDICTION MATRIX ---
  pred_matrix <- do.call(cbind, pred_list)
  colnames(pred_matrix) <- names(pred_list)

  
  # --- COMPUTE SSE FOR EACH MODEL ---
  sse_values <- compute_sse(
    qx = cohort_data$qx,
    pred_qx = pred_matrix
  )
  
  # --- CREATE RESULTS TABLE ---
  results <- data.frame(
    Method = names(pred_list),
    SSE = sse_values,
    stringsAsFactors = FALSE
  )
  
  # Sort by SSE (best model first)
  results <- results[order(results$SSE), ]
  rownames(results) <- NULL
  
  # --- FORMAT OUTPUT ---
  if(output_format == "raw") {
    return(results)
  } else {
    return(knitr::kable(results, 
                        format = "pipe", 
                        digits = 4,
                        caption = paste("Model Comparison (SSE) - Cohort", unique(cohort_data$Cohort))))
  }
}

# ANALISIS KOHOR TUNGGAL (1901)
cat("\nANALISIS KOHOR WANITA 1901\n")
target_cohort <- 1901
kohor_target <- data_females %>% filter(Cohort == target_cohort)
cat("Rentang usia:", min(kohor_target$Age), "-", max(kohor_target$Age), "dengan", nrow(kohor_target), "observasi\n")

models_1901 <- fit_static_models_updated(kohor_target, verbose = TRUE)

# PLOT TLT vs STLT
create_tlt_stlt_fitted_plot <- function(cohort_data, models, title_suffix = "") {
  if(is.null(models$stlt) || is.null(models$tlt)) return(NULL)
  
  fitted_df <- tibble(
    Age = cohort_data$Age,
    qx_stlt = predict(models$stlt, newdata = cohort_data$Age),
    qx_tlt = predict.tlt(models$tlt, newdata = cohort_data$Age)
  ) %>%
    pivot_longer(cols = c(qx_stlt, qx_tlt), names_to = "Model", values_to = "qx_fitted") %>%
    mutate(Model = case_when(Model == "qx_stlt" ~ "STLT", Model == "qx_tlt" ~ "TLT", TRUE ~ Model))
  
  p <- ggplot() +
    geom_point(data = cohort_data, aes(x = Age, y = qx, size = lx), alpha = 0.7, shape = 1, color = "black") +
    geom_line(data = fitted_df, aes(x = Age, y = qx_fitted, color = Model), linewidth = 1.1) +
    scale_color_manual(name = NULL, values = c("STLT" = "red", "TLT" = "blue")) +
    scale_size_continuous(name = "Population (lx)", guide = "none") +
    labs(title = paste("Perbandingan TLT dan STLT", title_suffix),
         subtitle = paste("Kohor Wanita Belanda", unique(cohort_data$Cohort)),
         x = "Usia", y = "Tingkat Mortalitas (qx)") +
    theme_bw() +
    theme(legend.position = c(0.2, 0.8),
          legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12)) +
    coord_cartesian(xlim = c(min(cohort_data$Age), max(cohort_data$Age) + 2), ylim = c(0, 1.05))
  
  return(p)
}

# PLOT SEMUA MODEL
create_all_models_fitted_plot <- function(cohort_data, models, title_suffix = "") {
  if(is.null(models$stlt)) return(NULL)
  
  fitted_list <- list()
  fitted_list$STLT <- predict(models$stlt, newdata = cohort_data$Age)
  if(!is.null(models$tlt)) fitted_list$TLT <- predict.tlt(models$tlt, newdata = cohort_data$Age)
  
  model_names <- c("gompertz", "makeham", "HP2", "coale_kisker")
  model_labels <- c("Gompertz", "Makeham", "Heligman-Pollard", "Coale-Kisker")
  
  for(i in seq_along(model_names)) {
    if(!is.null(models[[model_names[i]]])) {
      fitted_list[[model_labels[i]]] <- models[[model_names[i]]]
    }
  }
  
  fitted_df <- map_dfr(names(fitted_list), function(model_name) {
    tibble(Age = cohort_data$Age, qx_fitted = fitted_list[[model_name]], Model = model_name)
  })
  
  p <- ggplot() +
    geom_point(data = cohort_data, aes(x = Age, y = qx, size = lx), alpha = 0.6, shape = 1, color = "black") +
    geom_line(data = fitted_df, aes(x = Age, y = qx_fitted, color = Model), linewidth = 1) +
    scale_color_manual(name = "Model",
                       values = c("STLT" = "red", "TLT" = "darkblue", 
                                  "Gompertz" = "blue", "Makeham" = "green", 
                                  "Heligman-Pollard" = "purple", "Coale-Kisker" = "orange")) +
    scale_size_continuous(name = "Population (lx)", guide = "none") +
    labs(title = paste("Perbandingan Berbagai Model", title_suffix),
         subtitle = paste("Kohor Wanita Belanda", unique(cohort_data$Cohort)),
         x = "Usia", y = "Tingkat Mortalitas (qx)") +
    theme_bw() +
    theme(legend.position = "right",
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12)) +
    coord_cartesian(xlim = c(min(cohort_data$Age), max(cohort_data$Age) + 2), ylim = c(0, 1.05))
  
  return(p)
}

plot_tlt_stlt <- create_tlt_stlt_fitted_plot(kohor_target, models_1901)
plot_all_models <- create_all_models_fitted_plot(kohor_target, models_1901)
if(!is.null(plot_tlt_stlt)) print(plot_tlt_stlt)
if(!is.null(plot_all_models)) print(plot_all_models)

# ANALISIS MULTI-KOHOR UNTUK DSTLT
cat("\nANALISIS MULTI-KOHOR UNTUK DSTLT\n")

prepare_dstlt_data <- function(data, cohorts, sex = "Female") {
  filtered_data <- data %>% filter(Sex == sex, Cohort %in% cohorts) %>% arrange(Cohort, Age)
  cohort_list <- sort(unique(filtered_data$Cohort))
  age_range <- sort(unique(filtered_data$Age))
  
  ages_matrix <- matrix(NA, nrow = length(age_range), ncol = length(cohort_list))
  qx_matrix <- matrix(NA, nrow = length(age_range), ncol = length(cohort_list))
  
  for(i in seq_along(cohort_list)) {
    cohort_data <- filtered_data %>% filter(Cohort == cohort_list[i])
    for(j in seq_along(age_range)) {
      age_idx <- which(age_range == age_range[j])
      cohort_age_data <- cohort_data %>% filter(Age == age_range[j])
      if(nrow(cohort_age_data) > 0) {
        ages_matrix[age_idx, i] <- cohort_age_data$Age[1]
        qx_matrix[age_idx, i] <- cohort_age_data$qx[1]
      }
    }
  }
  
  return(list(ages = ages_matrix, qx = qx_matrix, cohorts = cohort_list, sex = sex))
}

training_cohorts <- 1893:1901
test_cohorts <- 1902:1908
cat("Menyiapkan data DSTLT untuk kohor", min(training_cohorts), "-", max(training_cohorts), "\n")

dstlt_data_female <- prepare_dstlt_data(data_females, training_cohorts, "Female")

cat("Fitting DSTLT...\n")
dstlt_fit_female <- tryCatch({
  dstlt(ages = dstlt_data_female$ages, qxs = dstlt_data_female$qx, startN = 85, endN = 105, hessian = FALSE)
}, error = function(e) {
  cat("DSTLT gagal:", e$message, "\n")
  NULL
})

if(!is.null(dstlt_fit_female)) {
  cat("DSTLT berhasil - N:", dstlt_fit_female$coefficients$N, "Omega:", round(dstlt_fit_female$Omega, 2), "\n")
}

cat("Fitting CBD...\n")
cbd_fit_female <- tryCatch({
  fit_cbd_model(ages_matrix = dstlt_data_female$ages, qx_matrix = dstlt_data_female$qx, 
                forecast_periods = length(test_cohorts), verbose = TRUE)
}, error = function(e) {
  cat("CBD gagal:", e$message, "\n")
  NULL
})

# TABEL PARAMETER
create_parameter_table <- function(models, cohort_name) {
  param_df <- data.frame(Parameter = c("B", "C", "Gamma", "N", "Omega"), STLT = NA, TLT = NA, stringsAsFactors = FALSE)
  
  if(!is.null(models$stlt)) {
    param_df$STLT <- c(models$stlt$coefficients$B, models$stlt$coefficients$C, 
                       models$stlt$coefficients$gamma, models$stlt$coefficients$N, models$stlt$Omega)
  }
  
  if(!is.null(models$tlt)) {
    tlt_coeffs <- models$tlt$coefficients
    param_df$TLT <- c(as.numeric(tlt_coeffs["B"]), as.numeric(tlt_coeffs["C"]), 
                      as.numeric(tlt_coeffs["gamma"]), as.numeric(tlt_coeffs[5]),
                      as.numeric(tlt_coeffs[5]) - as.numeric(tlt_coeffs["theta"])/as.numeric(tlt_coeffs["gamma"]))
  }
  
  return(param_df)
}

# TREN PARAMETER
create_parameter_trends_plot <- function(data, cohorts, sex = "Female") {
  param_trends <- data.frame(Cohort = integer(), B = numeric(), C = numeric(), 
                             Gamma = numeric(), N = integer(), Omega = numeric())
  
  for(cohort in cohorts) {
    cohort_data <- data %>% filter(Sex == sex, Cohort == cohort)
    if(nrow(cohort_data) > 10) {
      tryCatch({
        stlt_fit <- stlt(ages = cohort_data$Age, qx = cohort_data$qx)
        param_trends <- rbind(param_trends, data.frame(
          Cohort = cohort, B = stlt_fit$coefficients$B, C = stlt_fit$coefficients$C,
          Gamma = stlt_fit$coefficients$gamma, N = stlt_fit$coefficients$N, Omega = stlt_fit$Omega
        ))
      }, error = function(e) {
        cat("Gagal fitting kohor", cohort, "\n")
      })
    }
  }
  
  if(nrow(param_trends) < 3) return(NULL)
  
  param_long <- param_trends %>%
    select(Cohort, B, C, Gamma, N) %>%
    pivot_longer(cols = c(B, C, Gamma, N), names_to = "Parameter", values_to = "Value")
  
  p <- ggplot(param_long, aes(x = Cohort, y = Value)) +
    geom_line(color = "blue", linewidth = 1) +
    geom_point(color = "red", size = 2) +
    facet_wrap(~Parameter, scales = "free_y", nrow = 2) +
    labs(title = "Tren Parameter STLT Antar Kohor", subtitle = paste("Kohor", sex, "Belanda"),
         x = "Kohor Kelahiran", y = "Nilai Parameter") +
    theme_bw() +
    theme(plot.title = element_text(size = 14, face = "bold"), plot.subtitle = element_text(size = 12),
          strip.text = element_text(face = "bold"))
  
  return(list(plot = p, data = param_trends))
}

# PERBANDINGAN DSTLT vs CBD
# NOTE: Fungsi ini didefinisikan ulang di bawah (line ~639) dengan error handling yang lebih baik

# GENERATE PLOTS DAN TABEL
if(length(unique(data_females$Cohort)) >= 3) {
  param_trends_result <- create_parameter_trends_plot(data_females, unique(data_females$Cohort), "Female")
  if(!is.null(param_trends_result)) {
    plot_param_trends <- param_trends_result$plot
    param_trends_data <- param_trends_result$data
    print(plot_param_trends)
  }
}

if(!is.null(dstlt_fit_female) && !is.null(cbd_fit_female)) {
  test_data_female <- data_females %>% filter(Cohort %in% test_cohorts)
  dstlt_cbd_result <- create_dstlt_cbd_comparison_plot(dstlt_fit_female, cbd_fit_female, test_data_female, test_cohorts, training_cohorts)
  
  if(!is.null(dstlt_cbd_result)) {
    plot_dstlt_cbd <- dstlt_cbd_result$plot
    print(plot_dstlt_cbd)
    
    pred_error_data <- dstlt_cbd_result$data
    error_summary <- pred_error_data %>%
      group_by(Cohort) %>%
      summarise(DSTLT_MAE = mean(abs(qx_dstlt - qx_observed), na.rm = TRUE),
                DSTLT_RMSE = sqrt(mean((qx_dstlt - qx_observed)^2, na.rm = TRUE)),
                CBD_MAE = mean(abs(qx_cbd - qx_observed), na.rm = TRUE),
                CBD_RMSE = sqrt(mean((qx_cbd - qx_observed)^2, na.rm = TRUE)), .groups = 'drop')
  }
}

# CETAK SEMUA TABEL
cat("\n=== TABEL UNTUK SKRIPSI ===\n")

cat("\nTabel 1 - Estimasi Parameter:\n")
param_table <- create_parameter_table(models_1901, "1901")
print(param_table)

cat("\n=== Tabel 2 - Perbandingan Model Statis (SSE) ===\n")
comparison_sse <- generate_comparison_table_fixed(kohor_target, models_1901, output_format = "raw")

if(!is.null(comparison_sse)) {
  print(comparison_sse)
  
  # Tampilkan model terbaik
  best_model <- comparison_sse$Method[which.min(comparison_sse$SSE)]
  cat(sprintf("\n✓ Model terbaik untuk kohor %d: %s (SSE = %.4f)\n\n", 
              target_cohort, best_model, min(comparison_sse$SSE, na.rm = TRUE)))}

if(exists("param_trends_data")) {
  cat("\nTabel 3 - Tren Parameter:\n")
  print(param_trends_data)
}

if(!is.null(dstlt_fit_female)) {
  cat("\nTabel 5 - Parameter DSTLT:\n")
  dstlt_param_table <- data.frame(
    Parameter = c("a", "b", "theta", "gamma", "N", "Omega"),
    Estimate = c(dstlt_fit_female$coefficients$a, dstlt_fit_female$coefficients$b,
                 dstlt_fit_female$coefficients$theta, dstlt_fit_female$coefficients$gamma,
                 dstlt_fit_female$coefficients$N, dstlt_fit_female$Omega)
  )
  print(dstlt_param_table)
}

if(exists("error_summary")) {
  cat("\nTabel 8 - Error Prediksi:\n")
  print(error_summary)
}




create_dstlt_cbd_comparison_plot <- function(dstlt_model, cbd_model, test_data, test_cohorts, training_cohorts) {
  if(is.null(dstlt_model) || is.null(cbd_model)) return(NULL)

  pred_results <- list()
  n_training <- length(training_cohorts)  # 9 untuk kohor 1893-1901
  
  cat("Memprediksi kohor test:\n")
  
  for(i in seq_along(test_cohorts)) {
    cohort <- test_cohorts[i]
    cohort_data <- test_data %>% filter(Cohort == cohort)
    if(nrow(cohort_data) == 0) next

    cat(sprintf("  Kohor %d (i=%d, t=%d)...", cohort, i, n_training + i))

    # Parameter t untuk predict.dstlt:
    # t=1 → kohor training pertama (1893)
    # t=n_training → kohor training terakhir (1901)
    # t=n_training+1 → kohor test pertama (1902)
    # t=n_training+i → kohor test ke-i
    dstlt_pred <- tryCatch({
      result <- predict.dstlt(dstlt_model, newdata = cohort_data$Age, t = n_training + i)
      cat(" DSTLT OK")
      result
    }, error = function(e) {
      cat(" DSTLT GAGAL - Error:", e$message)
      rep(NA, nrow(cohort_data))
    })
    
    # Prediksi CBD
    # forecast_period=i → periode forecast ke-i (untuk kohor test ke-i)
    cbd_pred <- tryCatch({
      result <- predict.cbd(cbd_model, newdata = cohort_data$Age, forecast_period = i)
      cat(", CBD OK\n")
      result
    }, error = function(e) {
      cat(", CBD GAGAL - Error:", e$message, "\n")
      rep(NA, nrow(cohort_data))
    })
    
    pred_results[[paste0("cohort_", cohort)]] <- data.frame(
      Age = cohort_data$Age, 
      qx_observed = cohort_data$qx, 
      qx_dstlt = dstlt_pred,
      qx_cbd = cbd_pred, 
      Cohort = cohort, 
      lx = cohort_data$lx
    )
  }
  
  if(length(pred_results) == 0) return(NULL)
  
  all_pred <- do.call(rbind, pred_results)
  
  # Debug: Cek apakah ada nilai yang aneh
  cat("\nDEBUG - Statistik prediksi:\n")
  cat("  DSTLT range:", min(all_pred$qx_dstlt, na.rm=T), "-", max(all_pred$qx_dstlt, na.rm=T), "\n")
  cat("  CBD range:", min(all_pred$qx_cbd, na.rm=T), "-", max(all_pred$qx_cbd, na.rm=T), "\n")
  cat("  Observed range:", min(all_pred$qx_observed, na.rm=T), "-", max(all_pred$qx_observed, na.rm=T), "\n")
  
  pred_long <- all_pred %>%
    select(Age, qx_observed, qx_dstlt, qx_cbd, Cohort, lx) %>%
    pivot_longer(cols = c(qx_dstlt, qx_cbd), names_to = "Model", values_to = "qx_pred") %>%
    mutate(Model = case_when(Model == "qx_dstlt" ~ "DSTLT", Model == "qx_cbd" ~ "CBD", TRUE ~ Model)) %>%
    filter(!is.na(qx_pred))
  
  p <- ggplot() +
    geom_point(data = all_pred, aes(x = Age, y = qx_observed, size = log(lx + 1)), 
               alpha = 0.6, shape = 1, color = "grey50") +
    geom_line(data = pred_long, aes(x = Age, y = qx_pred, color = Model), linewidth = 1) +
    facet_wrap(~Cohort, scales = "free_y", nrow = 2) +
    scale_color_manual(name = "Model", values = c("DSTLT" = "red", "CBD" = "blue")) +
    scale_size_continuous(guide = "none") +
    labs(title = "Perbandingan Prediksi DSTLT vs CBD", 
         subtitle = paste("Kohor Test Wanita Belanda", min(test_cohorts), "-", max(test_cohorts)),
         x = "Usia", y = "qx") +
    theme_bw() +
    theme(legend.position = "bottom")
  
  return(list(plot = p, data = all_pred))
}

cat("\n=== INSTRUKSI PENGGUNAAN ===\n")
cat("1. Source file ini SETELAH source kode utama Anda\n")
cat("2. Fungsi create_dstlt_cbd_comparison_plot sudah diperbaiki dengan parameter training_cohorts\n")
cat("3. Jalankan ulang bagian prediksi:\n\n")
cat("   result <- create_dstlt_cbd_comparison_plot(dstlt_fit_female, cbd_fit_female, \n")
cat("                                               data_females %>% filter(Cohort %in% test_cohorts), \n")
cat("                                               test_cohorts, training_cohorts)\n")
cat("   print(result$plot)\n\n")
cat("4. Lihat output DEBUG untuk memastikan prediksi masuk akal\n")
cat("5. Jika masih error, salin error message lengkapnya\n\n")










