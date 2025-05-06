#' Predict method for lassovar objects
#'
#' @param object a lassovar object.
#' @param newdata a tibble or matrix with latest observations for prediction.
#' @param ... not used.
#' @return A tibble with predictions.
#' @export
predict.lassovar <- function(object, newdata, ...){
  # Validate input
  if(missing(newdata)){
    stop("newdata must be provided")
  }
  
  # Convert to tibble
  newdata <- tibble::as_tibble(newdata)
  
  # Extract model parameters
  p <- object$p
  intercept <- object$intercept
  beta <- object$coefficients
  
  # Prepare design matrix from newdata
  X <- prepare_predict_design(newdata, p, intercept)
  
  # Generate predictions
  predictions <- tibble::as_tibble(X %*% beta)
  
  return(predictions)
}

#' Create design matrix for prediction
#'
#' @param newdata Data for generating prediction
#' @param p Lag order
#' @param intercept Whether to include intercept
#' @return Design matrix
#' @keywords internal
prepare_predict_design <- function(newdata, p, intercept){
  # Convert to matrix for easier handling
  newdata_mat <- as.matrix(newdata)
  
  # Get dimensions
  n <- nrow(newdata_mat)
  k <- ncol(newdata_mat)
  
  # Check if we have enough observations
  if(n < p){
    stop("newdata must have at least p observations")
  }
  
  # Create lagged matrices
  lagged_data <- purrr::map(1:p, function(lag) {
    # Create lagged data
    lag_mat <- newdata_mat[(p-lag+1):(n-lag), , drop = FALSE]
    return(lag_mat)
  })
  
  # Combine all lags into design matrix
  X <- do.call(cbind, lagged_data)
  
  # Add intercept if requested
  if(intercept) {
    X <- cbind(1, X)
  }
  
  return(X)
}