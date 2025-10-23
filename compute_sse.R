#' @title Compute Sum of Squared Errors (SSE)
#' @description Calculate SSE for model comparison as in Huang et al. (2020)
#' 
#' @param qx Vector of observed mortality rates
#' @param pred_qx Vector or matrix of predicted mortality rates
#'                If matrix, each column is one model
#' 
#' @return Numeric vector of SSE values (one per model)
#' @export

compute_sse <- function(qx, pred_qx) {
  
  # Convert to matrix if vector
  if(is.vector(pred_qx)) {
    pred_qx <- matrix(pred_qx, ncol = 1)
  }
  
  # Validate dimensions
  if(nrow(pred_qx) != length(qx)) {
    stop(sprintf("Length mismatch: qx has %d values but pred_qx has %d rows", 
                 length(qx), nrow(pred_qx)))
  }
  
  # Compute SSE for each model
  n_models <- ncol(pred_qx)
  sse_values <- numeric(n_models)
  
  for(i in 1:n_models) {
    pred <- pred_qx[, i]
    
    # Remove NA and infinite values from BOTH observations and predictions
    valid_data <- !is.na(qx) & !is.na(pred) & 
      is.finite(qx) & is.finite(pred)
    
    n_valid <- sum(valid_data)
    
    if(n_valid < 2) {
      sse_values[i] <- NA
      warning(sprintf("Model %d (%s): only %d valid data points", 
                      i, 
                      ifelse(is.null(colnames(pred_qx)), "Unknown", colnames(pred_qx)[i]),
                      n_valid))
      next
    }
    
    # Calculate SSE: Σ(qx - q̂x)²
    errors <- qx[valid_data] - pred[valid_data]
    sse_values[i] <- sum(errors^2)
  }
  
  return(sse_values)
}