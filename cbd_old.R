#' @title Fit Coale-Kisker Model
#'
#' @description Fits the Coale-Kisker extrapolation method for mortality rates
#'
#' @param ages Vector of ages
#' @param qx Vector of observed mortality rates
#' @param mx Vector of central death rates (optional, will be estimated if not provided)
#' @param x0 Starting age for extrapolation (default 84)
#' @param x1 Ending age for extrapolation (default 110)
#' @param mx1 Central death rate at ending age (default 1.0)
#' @param pred_ages Ages for which to predict qx
#'
#' @return Vector of predicted qx values
#'
#' @export
fit_coale_kisker <- function(ages, qx, mx = NULL, x0 = 84, x1 = 110, mx1 = 1.0, pred_ages = ages) {
  
  # Convert qx to mx if mx not provided
  if(is.null(mx)) {
    # Approximate conversion: mx ≈ qx / (1 - 0.5*qx) for most ages
    mx <- qx / (1 - 0.5 * qx)
  }
  
  # Create data frame for easier manipulation
  mortality_data <- data.frame(
    age = ages,
    qx = qx,
    mx = mx
  )
  
  # Predict mx for all requested ages
  pred_mx <- numeric(length(pred_ages))
  
  for(i in seq_along(pred_ages)) {
    age <- pred_ages[i]
    
    if(age <= max(ages)) {
      # Use observed data for ages within range
      if(age %in% ages) {
        pred_mx[i] <- mx[ages == age]
      } else {
        # Interpolate for missing ages within observed range
        pred_mx[i] <- approx(ages, mx, xout = age, rule = 2)$y
      }
    } else if(age <= x0) {
      # Extrapolate using last observed mx
      pred_mx[i] <- mx[length(mx)]
    } else {
      # Apply Coale-Kisker method for ages > x0
      
      # Calculate k(x) values for observed data near x0
      mx_x0 <- approx(ages, mx, xout = x0, rule = 2)$y
      
      # Calculate R parameter
      if(x0 %in% ages) {
        k_x0 <- log(mx[ages == x0] / mx[ages == (x0+1)])
      } else {
        # Approximate k(x0) from nearby values
        mx_x0_plus1 <- approx(ages, mx, xout = x0+1, rule = 2)$y
        k_x0 <- log(mx_x0 / mx_x0_plus1)
      }
      
      # Calculate R using Coale-Kisker formula
      sum_i <- sum(1:(x1 - x0))
      R <- ((x1 - x0) * k_x0 + log(mx_x0) - log(mx1)) / sum_i
      
      # Calculate k(age) for current age
      k_age <- k_x0 - R * (age - x0)
      
      # Calculate mx for previous age
      if(age == x0 + 1) {
        mx_prev <- mx_x0
      } else {
        # Use recursion: mx(age-1) based on previous calculation
        mx_prev <- pred_mx[pred_ages == (age - 1)]
      }
      
      # Calculate mx(age) = mx(age-1) / exp(k(age))
      pred_mx[i] <- mx_prev / exp(k_age)
    }
  }
  
  # Convert mx back to qx: qx = mx / (1 + 0.5*mx)
  pred_qx <- pred_mx / (1 + 0.5 * pred_mx)
  
  # Ensure qx doesn't exceed 1
  pred_qx <- pmin(pred_qx, 1)
  
  return(pred_qx)
}