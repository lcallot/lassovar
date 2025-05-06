#' Forecasting from a lassovar object
#'
#' @param object a lassovar object.
#' @param h horizon of the forecast, default is 1.
#' @param newdata a tibble or matrix with latest observations for forecasting.
#' @param ... not used.
#' @return A tibble with forecasts.
#' @export
forecast.lassovar <- function(object, h = 1, newdata = NULL, ...){
  # Extract model parameters
  p <- object$p
  k <- object$k
  intercept <- object$intercept
  
  # Extract coefficients
  beta <- object$coefficients
  
  # If newdata is provided, use it for forecasting
  if(!is.null(newdata)){
    # Convert to tibble
    newdata <- tibble::as_tibble(newdata)
    
    # Check if newdata has enough observations
    if(nrow(newdata) < p){
      stop("newdata must contain at least p observations")
    }
    
    # Use the most recent p observations for forecasting
    xnew <- tibble::as_tibble(newdata[(nrow(newdata)-p+1):nrow(newdata),])
  } else {
    # Use the fitted data from the model
    xnew <- object$fitted[(nrow(object$fitted)-p+1):nrow(object$fitted),]
  }
  
  # Prepare storage for forecasts
  forecasts <- tibble::tibble(
    .rows = h,
    !!!rlang::set_names(rep(list(numeric(h)), k), colnames(xnew))
  )
  
  # Generate forecasts for each horizon
  for(i in 1:h){
    # Create design matrix for the current forecast step
    X <- create_forecast_design(xnew, p, intercept)
    
    # Generate forecast
    pred <- X %*% beta
    
    # Store forecast
    forecasts[i, ] <- pred
    
    # Update xnew for the next iteration
    if(i < h){
      xnew <- rbind(xnew[-1,], pred)
    }
  }
  
  return(forecasts)
}

#' Create design matrix for forecasting
#'
#' @param x Data for generating forecast
#' @param p Lag order
#' @param intercept Whether to include intercept
#' @return Design matrix
#' @keywords internal
create_forecast_design <- function(x, p, intercept){
  # Convert to matrix for easier handling
  x_mat <- as.matrix(x)
  
  # Create design matrix with lagged values
  design <- numeric()
  
  # Add lags
  for(lag in 1:p){
    design <- cbind(design, x_mat[nrow(x_mat)-lag+1,])
  }
  
  # Add intercept if needed
  if(intercept){
    design <- c(1, design)
  }
  
  return(matrix(design, nrow = 1))
}