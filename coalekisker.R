fit_coale_kisker <- function(ages, qx, mx = NULL, x0 = 85, x1 = 110, mx1 = 1.0) {
  # Validasi input
  if (length(ages) != length(qx)) stop("Panjang ages dan qx tidak sama.")
  if (x1 <= x0) stop("x1 harus > x0.")
  if (!any(ages >= x0) || !any(ages >= x0 + 1)) stop("Data tidak cukup di sekitar x0.")
  
  # Bersihkan duplikat dan NA
  valid_idx <- !is.na(ages) & !is.na(qx)
  ages <- ages[valid_idx]
  qx <- qx[valid_idx]
  ages <- unique(ages)
  qx <- qx[match(ages, ages)]
  
  pred_ages <- seq(min(ages), x1, by = 1)
  
  if (is.null(mx)) mx <- qx / (1 - 0.5 * qx)
  mx[mx <= 0 | is.na(mx)] <- 1e-10
  mx <- pmin(mx, 10)
  
  # Interpolasi di x0 dan x0+1
  mx_x0_val <- approx(ages, mx, xout = x0, rule = 2)$y
  mx_x0_plus1_val <- approx(ages, mx, xout = x0 + 1, rule = 2)$y
  if (any(is.na(c(mx_x0_val, mx_x0_plus1_val))) || mx_x0_val <= 0 || mx_x0_plus1_val <= 0) {
    warning("Nilai mx tidak valid di x0/x0+1.")
    return(rep(NA_real_, length(pred_ages)))
  }
  
  k_x0 <- log(mx_x0_val / mx_x0_plus1_val)
  denominator_R <- (x1 - x0) * (x1 - x0 + 1) / 2
  numerator_R <- (x1 - x0) * k_x0 + log(mx_x0_val) - log(mx1)
  R <- numerator_R / denominator_R
  
  pred_mx <- numeric(length(pred_ages))
  k_vals <- numeric(length(pred_ages))
  idx_x0 <- which(pred_ages == x0)
  
  for (i in 1:(idx_x0 - 1)) {
    pred_mx[i] <- approx(ages, mx, xout = pred_ages[i], rule = 2)$y
    k_vals[i] <- NA
  }
  
  pred_mx[idx_x0] <- mx_x0_val
  k_vals[idx_x0] <- k_x0
  
  if (x1 > x0) {
    for (i in (idx_x0 + 1):length(pred_ages)) {
      k_vals[i] <- k_vals[i - 1] - R
      prev_mx <- pred_mx[i - 1]
      prev_k <- k_vals[i - 1]
      if (is.na(prev_mx) || prev_mx <= 0) {
        pred_mx[i:length(pred_ages)] <- NA
        break
      }
      pred_mx[i] <- prev_mx * exp(-prev_k)
      pred_mx[i] <- min(pred_mx[i], 10)
    }
  }
  
  pred_qx <- pred_mx / (1 + 0.5 * pred_mx)
  pred_qx <- pmin(pmax(pred_qx, 0), 1)
  
  if (pred_ages[length(pred_qx)] == x1) pred_qx[length(pred_qx)] <- 1.0
  
  return(pred_qx)
}
