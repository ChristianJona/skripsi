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
    k_x0 <- log(mx_x0 / mx_x0_plus1)
  }
  
  sum_i <- sum(1:(x1 - x0))
  R <- ((x1 - x0) * k_x0 + log(mx_x0) - log(mx1)) / sum_i
  
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
        pred_mx[i] <- pred_mx[i-1] / exp(k_age)
      }
    }
  }
  
  # Convert mx to qx
  pred_qx <- pred_mx / (1 + 0.5 * pred_mx)
  pred_qx <- pmin(pred_qx, 1)
  
  return(pred_qx)
}