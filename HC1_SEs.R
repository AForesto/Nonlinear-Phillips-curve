############################ FUNCTION TO ESTIMATE HETEROSKEDASTICITY AND AUTOCORRELATION ADJUSTED STANDARD ERRORS ############

#' @X matrix or vector of independent variables
#' @residuals vector of error values from the regression
#' @bandwith the number of lags considered in the covariance calculation, default to Newey-West
HAC_se <- function(X, residuals, bandwidth = NULL) {
  n <- length(residuals)
  k <- ncol(X)
  
  # Default bandwidth using Newey-West optimal bandwidth
  if (is.null(bandwidth)) {
    bandwidth <- floor(4 * (n/100)^(2/9))
  }
  
  # Initialize the middle term of the sandwich estimator
  Omega <- matrix(0, k, k)
  
  # Compute X'ee'X accounting for autocorrelation
  for (j in 0:bandwidth) {
    # Compute autocorrelation at lag j
    Gamma_j <- matrix(0, k, k)
    for (t in (j+1):n) {
      xee <- X[t,] * residuals[t]
      if (j > 0) {
        xee_lag <- X[t-j,] * residuals[t-j]
      } else {
        xee_lag <- xee
      }
      Gamma_j <- Gamma_j + (xee %*% t(xee_lag))
    }
    
    # Apply Bartlett kernel
    if (j > 0) {
      weight <- 1 - j/(bandwidth + 1)
      Omega <- Omega + weight * (Gamma_j + t(Gamma_j))
    } else {
      Omega <- Omega + Gamma_j
    }
  }
  
  # Compute bread of sandwich
  bread <- solve(t(X) %*% X)
  
  # Compute HAC variance-covariance matrix
  V_HAC <- bread %*% Omega %*% bread
  
  # Return HAC standard errors
  return(sqrt(diag(V_HAC)))
}

  
