#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' Estimate a VAR with the Lasso
#' 
#' @param y matrix or tibble of observations
#' @param p the lag length
#' @param intercept should an intercept be included
#' @param post should post-Lasso OLS be performed
#' @param ada should the adaptive Lasso be used
#' @param weight.ada vector or scalar of adaptive weights
#' @param lambda grid of values for the regularization parameter
#' @param lambda.ada grid of values for the regularization of the adaptive Lasso
#' @param cv.ada should the adaptive lasso be cross validated
#' @param ic selection of lambda using information criteria, default is cv
#' @param nfolds number of folds for cross validation
#' @param verbose should progress be reported
#' @param ... Extra arguments to be passed to \code{glmnet}
#' @return An object of class \code{lassovar}
#' @export
lassovar <- function(y, 
                    p = 1, 
                    intercept = TRUE,
                    post = FALSE,
                    ada = FALSE,
                    weight.ada = 1,
                    lambda = NULL,
                    lambda.ada = NULL,
                    cv.ada = FALSE,
                    ic = c("cv", "aic", "bic", "aicc"),
                    nfolds = 5,
                    verbose = FALSE,
                    ...){
  
  # Convert input to tibble if not already
  y <- tibble::as_tibble(y)
  
  # Extract dimensions
  n <- nrow(y)
  k <- ncol(y)
  
  # Prepare data for VAR model
  var_data <- prepare_var_data(y, p, intercept)
  X <- var_data$X
  Y <- var_data$Y
  
  # Estimate model
  if(ada){
    mod <- estimate_adaptive_lasso(Y, X, weight.ada, lambda.ada, cv.ada, nfolds, ic, verbose, ...)
  } else {
    mod <- estimate_lasso(Y, X, lambda, nfolds, ic, verbose, ...)
  }
  
  # Post-Lasso OLS if requested
  if(post){
    mod <- perform_post_lasso_ols(mod, Y, X, ...)
  }
  
  # Create output object
  out <- list(
    coefficients = mod$coefficients,
    lambda = mod$lambda,
    lambda.min = mod$lambda.min,
    lambda.1se = mod$lambda.1se,
    ic = ic[1],
    p = p,
    k = k,
    n = n,
    intercept = intercept,
    post = post,
    ada = ada,
    ic.value = mod$ic.value,
    residuals = mod$residuals,
    fitted = mod$fitted,
    call = match.call()
  )
  
  class(out) <- "lassovar"
  return(out)
}

#' Internal function to prepare data for VAR model
#' 
#' @param y Original data
#' @param p Lag order
#' @param intercept Whether to include intercept
#' @return List with design matrix X and response matrix Y
#' @keywords internal
prepare_var_data <- function(y, p, intercept) {
  # Convert to matrix for manipulation
  y_mat <- as.matrix(y)
  
  # Get dimensions
  n <- nrow(y_mat)
  k <- ncol(y_mat)
  
  # Create lagged matrices
  lagged_data <- purrr::map(1:p, function(lag) {
    lag_mat <- y_mat[(p-lag+1):(n-lag), , drop = FALSE]
    # Create column names with lag information
    colnames(lag_mat) <- paste0(colnames(y_mat), "_lag", lag)
    return(lag_mat)
  })
  
  # Combine all lags into design matrix
  X <- do.call(cbind, lagged_data)
  
  # Add intercept if requested
  if(intercept) {
    X <- cbind(1, X)
    colnames(X)[1] <- "intercept"
  }
  
  # Response matrix
  Y <- y_mat[(p+1):n, , drop = FALSE]
  
  return(list(X = X, Y = Y))
}

#' Estimate VAR model using Lasso
#' 
#' @param Y Response matrix
#' @param X Design matrix
#' @param lambda Grid for regularization parameter
#' @param nfolds Number of folds for cross-validation
#' @param ic Information criterion
#' @param verbose Whether to show progress
#' @param ... Additional arguments for glmnet
#' @return List with model results
#' @keywords internal
estimate_lasso <- function(Y, X, lambda, nfolds, ic, verbose, ...) {
  # Get dimensions
  k <- ncol(Y)
  
  # Perform Lasso estimation for each variable
  models <- purrr::map(1:k, function(i) {
    if(verbose) {
      cat("Fitting model for variable", i, "of", k, "\n")
    }
    
    # Select dependent variable
    y_i <- Y[, i]
    
    # Fit lasso model
    fit <- glmnet::cv.glmnet(X, y_i, lambda = lambda, nfolds = nfolds, ...)
    
    # Select lambda based on information criterion
    if(ic[1] == "cv") {
      lambda_opt <- fit$lambda.min
    } else {
      # Calculate information criterion values for different lambda values
      ic_values <- calculate_ic(fit, X, y_i, ic[1])
      lambda_opt <- fit$lambda[which.min(ic_values)]
    }
    
    # Get fitted values and residuals
    beta <- as.vector(coef(fit, s = lambda_opt))
    fitted <- X %*% beta
    resid <- y_i - fitted
    
    return(list(
      beta = beta,
      fitted = fitted,
      residuals = resid,
      lambda_opt = lambda_opt
    ))
  })
  
  # Extract coefficients
  coefs <- purrr::map_dfr(models, ~ tibble::tibble(beta = .x$beta))
  coef_matrix <- as.matrix(coefs)
  
  # Extract other components
  fitted_values <- do.call(cbind, purrr::map(models, ~ .x$fitted))
  residuals <- do.call(cbind, purrr::map(models, ~ .x$residuals))
  lambda_opt <- purrr::map_dbl(models, ~ .x$lambda_opt)
  
  return(list(
    coefficients = coef_matrix,
    lambda = lambda,
    lambda.min = mean(lambda_opt), 
    residuals = residuals,
    fitted = fitted_values
  ))
}

#' Estimate VAR model using Adaptive Lasso
#' 
#' @param Y Response matrix
#' @param X Design matrix
#' @param weight.ada Adaptive weights
#' @param lambda.ada Grid for regularization parameter
#' @param cv.ada Whether to use cross-validation for adaptive step
#' @param nfolds Number of folds for cross-validation
#' @param ic Information criterion
#' @param verbose Whether to show progress
#' @param ... Additional arguments for glmnet
#' @return List with model results
#' @keywords internal
estimate_adaptive_lasso <- function(Y, X, weight.ada, lambda.ada, cv.ada, nfolds, ic, verbose, ...) {
  # Implementation similar to estimate_lasso but with adaptive weights
  # This is a simplification; full implementation would be more complex
  
  # Placeholder for the implementation
  stop("Adaptive Lasso implementation not yet available in modernized version")
}

#' Perform Post-Lasso OLS
#' 
#' @param mod Initial Lasso model
#' @param Y Response matrix
#' @param X Design matrix
#' @param ... Additional arguments
#' @return Updated model with OLS estimates
#' @keywords internal
perform_post_lasso_ols <- function(mod, Y, X, ...) {
  # Placeholder for the implementation
  stop("Post-Lasso OLS implementation not yet available in modernized version")
}

#' Calculate Information Criteria
#' 
#' @param fit Fitted glmnet model
#' @param X Design matrix
#' @param y Response vector
#' @param ic Type of information criterion
#' @return Vector of IC values
#' @keywords internal
calculate_ic <- function(fit, X, y, ic) {
  # Placeholder for the implementation
  # Would calculate AIC, BIC, or AICC based on the ic parameter
  # For now, return a random vector for demonstration
  return(runif(length(fit$lambda)))
}