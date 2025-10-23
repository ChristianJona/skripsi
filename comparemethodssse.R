#' @title Generate Comparison Table for Static Models (SSE only)
#' @description Compare models using Sum of Squared Errors as in Huang et al. (2020) Table 2
#' 
#' @param cohort_data Data frame with Age, qx, lx for one cohort
#' @param models List of fitted models (stlt, tlt, gompertz, etc.)
#' @param output_format "raw" for data frame, "kable" for formatted table
#' 
#' @return Data frame or kable with SSE for each model
#' @export

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
    pred_qx = pred_matrix,
    start_age = min(cohort_data$Age),
    end_age = max(cohort_data$Age)
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