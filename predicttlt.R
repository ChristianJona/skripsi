#' @title Predict using a fitted Threshold Life Table (TLT)
#'
#' @description For a given TLT fit, this function computes the qx associated with any
#' particular age. Unlike STLT, TLT may have discontinuity at the threshold age N.
#'
#' @param object a TLT model object
#'
#' @param newdata Vector of ages for which qx are to be predicted
#'
#' @param ... Additional arguments passed to methods
#'
#' @return Vector of predicted qx under the fitted TLT model
#'
#' @export
predict.tlt <- function(object, newdata, ...) {
  
  # Extract parameters from named vector
  coeffs <- object$coefficients
  B <- coeffs["B"]
  C <- coeffs["C"]
  theta <- coeffs["theta"]
  gam <- coeffs["gamma"]
  N <- coeffs[5]  # The 5th element (unnamed, likely N)
  
  # Remove names to avoid issues
  names(B) <- NULL
  names(C) <- NULL
  names(theta) <- NULL
  names(gam) <- NULL
  names(N) <- NULL
  
  # Calculate omega (highest attained age)
  if (gam < 0) {
    omega <- N - theta/gam  # Since gam is negative, this adds theta/|gam|
  } else {
    omega <- Inf  # No upper limit if gamma >= 0
  }
  
  # Initialize output vector
  qx <- rep(NA, length(newdata))
  
  # Loop through each age to predict
  for (i in 1:length(newdata)) {
    x <- newdata[i]
    
    if (x < N) {
      # Gompertz part (before threshold)
      # qx = (S(x) - S(x+1)) / S(x)
      # where S(x) = exp(-B/ln(C) * (C^x - 1))
      S_x <- exp(-B/log(C) * (C^x - 1))
      S_x1 <- exp(-B/log(C) * (C^(x+1) - 1))
      qx[i] <- (S_x - S_x1) / S_x
      
    } else {
      # GPD part (after threshold)
      if (abs(gam) > 1e-6) {  # gam != 0
        if (is.finite(omega) && x < omega - 1) {
          # Normal GPD case
          # S(x) = S(N) * (1 + gam*(x-N)/theta)^(-1/gam)
          S_N <- exp(-B/log(C) * (C^N - 1))  # Survival at threshold from Gompertz
          
          # Check if the GPD terms are valid
          term_x <- 1 + gam*(x - N)/theta
          term_x1 <- 1 + gam*(x + 1 - N)/theta
          
          if (term_x > 0 && term_x1 > 0) {
            S_x <- S_N * term_x^(-1/gam)
            S_x1 <- S_N * term_x1^(-1/gam)
            qx[i] <- (S_x - S_x1) / S_x
          } else {
            qx[i] <- 1  # Beyond valid range
          }
          
        } else if (is.finite(omega) && x < omega) {
          # Near the upper bound
          qx[i] <- 1
        } else if (is.infinite(omega)) {
          # No upper bound case (gamma > 0)
          S_N <- exp(-B/log(C) * (C^N - 1))
          term_x <- 1 + gam*(x - N)/theta
          term_x1 <- 1 + gam*(x + 1 - N)/theta
          
          if (term_x > 0 && term_x1 > 0) {
            S_x <- S_N * term_x^(-1/gam)
            S_x1 <- S_N * term_x1^(-1/gam)
            qx[i] <- (S_x - S_x1) / S_x
          } else {
            qx[i] <- NA
          }
        } else {
          # Beyond omega
          qx[i] <- NA
        }
        
      } else {
        # gam ≈ 0, exponential case
        S_N <- exp(-B/log(C) * (C^N - 1))
        S_x <- S_N * exp(-(x - N)/theta)
        S_x1 <- S_N * exp(-(x + 1 - N)/theta)
        qx[i] <- (S_x - S_x1) / S_x
      }
    }
  }
  
  return(qx)
}