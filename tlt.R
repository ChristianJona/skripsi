tlt <- function(ages, qx, lx, dx, start_N = 85, end_N = 100) {
  
  # --- Fungsi Internal: Log-Likelihood untuk bagian Gompertz (tetap sama) ---
  loglik_gompertz <- function(params, ages_gomp, dx_gomp, lx_gomp) {
    B <- exp(params[1])
    C <- 1 + exp(params[2])
    lnC <- log(C)
    S_x <- exp(-(B / lnC) * (C^ages_gomp - 1))
    S_x1 <- exp(-(B / lnC) * (C^(ages_gomp + 1) - 1))
    prob <- S_x - S_x1
    prob[prob <= 0] <- 1e-9
    ll <- sum(dx_gomp * log(prob)) + lx_gomp[length(lx_gomp)] * log(S_x1[length(S_x1)])
    return(-ll)
  }
  
  # --- Fungsi Internal: Log-Likelihood untuk bagian GPD (DENGAN PERBAIKAN) ---
  loglik_gpd <- function(params, ages_gpd, dx_gpd, lx_gpd, S_at_N) {
    theta <- exp(params[1])
    gamma <- params[2]
    
    # ---- PERBAIKAN DIMULAI DI SINI ----
    # 1. PENGECEKAN KEAMANAN PARAMETER
    # Untuk gamma < 0, pastikan data tidak melebihi batas atas teoretis GPD.
    # Jika melebihi, parameter ini tidak valid, beri penalti besar.
    if (gamma < 0) {
      max_excess_age <- max(ages_gpd) - min(ages_gpd) + 1
      if (max_excess_age >= -theta / gamma) {
        return(1e10) # Kembalikan nilai penalti yang sangat besar
      }
    }
    # ---- AKHIR DARI PERBAIKAN UTAMA ----
    
    excess_ages <- ages_gpd - min(ages_gpd)
    excess_ages_1 <- (ages_gpd + 1) - min(ages_gpd)
    
    if (gamma != 0) {
      surv_factor <- (1 + gamma * excess_ages / theta)^(-1 / gamma)
      surv_factor_1 <- (1 + gamma * excess_ages_1 / theta)^(-1 / gamma)
    } else {
      surv_factor <- exp(-excess_ages / theta)
      surv_factor_1 <- exp(-excess_ages_1 / theta)
    }
    
    S_x <- S_at_N * surv_factor
    S_x1 <- S_at_N * surv_factor_1
    
    prob <- S_x - S_x1
    
    # Pengecekan keamanan tambahan untuk probabilitas
    if (any(is.na(prob)) || any(prob <= 0)) {
      return(1e10)
    }
    
    ll <- sum(dx_gpd * log(prob)) + lx_gpd[length(lx_gpd)] * log(S_x1[length(S_x1)])
    
    if (is.infinite(ll) || is.na(ll)) {
      return(1e10)
    }
    
    return(-ll)
  }
  
  # --- Proses Utama: Mencari N Optimal (tetap sama) ---
  possible_N <- start_N:end_N
  total_logliks <- numeric(length(possible_N))
  
  for (i in seq_along(possible_N)) {
    N <- possible_N[i]
    data_gomp <- data.frame(ages, dx, lx) %>% filter(ages < N)
    data_gpd  <- data.frame(ages, dx, lx) %>% filter(ages >= N)
    
    if (nrow(data_gomp) < 2 || nrow(data_gpd) < 2) {
      total_logliks[i] <- -Inf
      next
    }
    
    fit_gomp <- optim(par = c(log(0.0001), log(0.1)), fn = loglik_gompertz,
                      ages_gomp = data_gomp$ages, dx_gomp = data_gomp$dx, lx_gomp = data_gomp$lx,
                      method = "Nelder-Mead")
    
    B_hat <- exp(fit_gomp$par[1])
    C_hat <- 1 + exp(fit_gomp$par[2])
    S_at_N <- exp(-(B_hat / log(C_hat)) * (C_hat^N - 1))
    
    fit_gpd <- optim(par = c(log(2), -0.1), fn = loglik_gpd,
                     ages_gpd = data_gpd$ages, dx_gpd = data_gpd$dx, lx_gpd = data_gpd$lx,
                     S_at_N = S_at_N, method = "Nelder-Mead")
    
    if (is.infinite(fit_gomp$value) || is.infinite(fit_gpd$value)) {
      total_logliks[i] <- -Inf
    } else {
      total_logliks[i] <- -(fit_gomp$value + fit_gpd$value)
    }
  }
  
  # ... (Sisa kode untuk fit ulang tetap sama) ...
  optimal_N <- possible_N[which.max(total_logliks)]
  data_gomp_final <- data.frame(ages, dx, lx) %>% filter(ages < optimal_N)
  data_gpd_final  <- data.frame(ages, dx, lx) %>% filter(ages >= optimal_N)
  fit_gomp_final <- optim(par = c(log(0.0001), log(0.1)), fn = loglik_gompertz,
                          ages_gomp = data_gomp_final$ages, dx_gomp = data_gomp_final$dx, lx_gomp = data_gomp_final$lx)
  B_final <- exp(fit_gomp_final$par[1])
  C_final <- 1 + exp(fit_gomp_final$par[2])
  S_at_N_final <- exp(-(B_final / log(C_final)) * (C_final^optimal_N - 1))
  fit_gpd_final <- optim(par = c(log(2), -0.1), fn = loglik_gpd,
                         ages_gpd = data_gpd_final$ages, dx_gpd = data_gpd_final$dx, lx_gpd = data_gpd_final$lx,
                         S_at_N = S_at_N_final)
  theta_final <- exp(fit_gpd_final$par[1])
  gamma_final <- fit_gpd_final$par[2]
  
  return(list(
    coefficients = c(B = B_final, C = C_final, theta = theta_final, gamma = gamma_final, N = optimal_N),
    logLik = max(total_logliks)
  ))
}