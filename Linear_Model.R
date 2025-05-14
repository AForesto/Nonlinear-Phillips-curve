############################# FUNCTION TO ESTIMATE A LINEAR MODEL ################################

#' @Y vector for the dependent variable
#' @X vector or matrix of independent variables
#' @const boolean that appends a value of 1's in the independent variable matrix
#' @conf_level numerical value between 0 and 0.99 that sets the level of confidence 

Linear_Model <- function(Y, X, conf_level = 0.95, HC1 = TRUE, const = TRUE){
  
  # Input validation
  if (!is.numeric(Y) || !is.matrix(X)){
    stop("Y must be a vector and X must be a matrix")
  }
  if (length(Y) != nrow(X)){
    stop("Number of observations in Y and X must be the same")
  }
  
  # Add constant term if requested and not already present
  if (const) {
    if (!all(X[,1] == 1)) {
      X <- cbind(1, X)
    }
  } else {
    if (all(X[,1] == 1)) {
      X <- X[,-1]  # Remove constant if present when const = FALSE
    }
  }
  # Create lags of Y and store in X
  L1 <- lag(Y, 1)
  L2 <- lag(Y, 2)
  Lags <- c(L1, L2)

  X <- cbind(Lags,X)
  
  # Calculate the coefficients (Betas)
  # Beta = (X'X)(-1)X'Y
  Beta <- solve(t(X) %*% X) %*% t(X) %*% Y
    
  # Fitted values
  # Y = XB
  Y_hat <- X %*% Beta
  # Residuals
  # e = Y - (BX)
  epsilon <- Y - Y_hat
    
  # Calculate standard errors
  if (HAC) {
    Std_error <- HAC_se(X, residuals = epsilon, bandwidth)
  } else {
    # Conventional standard errors
    N <- length(Y)                              # Number of observations
    K <- ncol(X)                                # Number of independent variables 
    df <- N - K                                 # Degrees of freedom
    MSE <- sum(epsilon^2)/df                    # Mean squared error
    VCOV <- as.numeric(MSE) * solve(t(X) %*% X) # Variance-covariance matrix: (Σe^2/df)* X'X
    Std_error <- sqrt(diag(VCOV))               # Standard error is the square root on the diagonal of the variance-covariance matrix
  }
    
  # Calculate t values
  # t = Beta/Standard error
  t_val <- as.vector(Beta)/Std_error
    
  # Calculate p values
  N <- length(Y)
  K <- ncol(X)
  df <- N - K
  p_val <- 2 * pt(abs(as.numeric(t_val)), df = df, lower.tail = F)
    
  # Calculate confidence bands
  # (Beta - critical value * standard error; Beta; Beta + critical value * standard error)
  crit_value <- qt(1 - (1 - conf_level)/2, df = df)
  conf_low <- Beta - crit_value * Std_error
  conf_high <- Beta + crit_value * Std_error
    
  # Diagnostic calculations 
  # Total sum of squares: TSS = Σ(y_i - y_bar)^2
  Y_bar <- mean(Y)
  TSS <- sum((Y - Y_bar)^2)
    
  # Residual sum of squares: RSS = Σ(y_i - y_hat_i)^2
  RSS <- sum(epsilon^2)
    
  # Explanatory sum of squares: ESS = TSS - RSS
  ESS <- TSS - RSS
    
  # R-squared
  R <- 1 - RSS/TSS
  # Adjusted R-Squared: R^2_adj = 1 - (RSS/(n-k))/(TSS/(n-1))
  R_adj <- 1 - (RSS/(N-K))/(TSS/(N-1))
    
  # F-statistic: F = (ESS/(k-1))/(RSS/(n-k))
  F_stat <- (ESS/(K-1))/(RSS/(N-K))
  F_p_val <- pf(F_stat, df1 = K-1, df2 = df, lower.tail = F)
    
  # return results
  diagnostics <- list(
    ` No. Observations` = N,
    `R squared` = R,
    `Adjusted R squared` = R_adj,
    `F statistic` = F_stat,
    `F p-value` = F_p_val
  )
    
  results <- list(
    Coefficients = as.vector(Beta),
    `Standard error` = as.vector(Std_error),
    `t statistic` = as.vector(t_val),
    `p value` =  as.vector(p_val),
    residuals = as.vector(epsilon),
    conf_low = as.vector(conf_low),
    conf_high = as.vector(conf_high),
    fitted = as.vector(Y_hat)
  )
    
return(list(results = results,
            diagnostics = diagnostics))
}
