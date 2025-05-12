############################# FUNCTION TO RUN MARKOV SWITCHING MODEL ################################3

#' @theta an object (list or matrix) of initial coefficients for the variables
#' @series_ts vector of the dependent variable
#' @controls matrix of exogenous control variables

Conditional_Log_Likelihood <- function(theta, series_ts, controls){
  
  # Calculate conditional log-likelihood for AR(2) MSM
  # 2 states: High and Low Inflation
  # Parameters that will change:
  # p11 probability to stay in state 1 (High Inflation)
  # p22 probability to stay in state 2 (Low Inflation)
  # c1 and c2, constant terms of the AR in states 1 and 2
  # Control variable parameters (beta's)
  # sigma, the standard deviation of disturbances
  # Function returns the log-likelihood of the model
  
  # Insert checks or print statements
  if(any(is.infinite(theta)) || any(is.nan(theta))) {
    print("Invalid theta values detected")
    print(theta)
    return(Inf)
  }
  
  # Ensure all inputs are numeric
  series <- as.numeric(series_ts)
  controls <- as.matrix(controls)
  theta <- as.numeric(theta)
  years <- floor(time(series_ts))
  
  p11 = theta[1]     # Transition probability: State 1
  p22 = theta[2]     # Transition probability: State 2
  c1 = theta[3]      # Constant in State 1
  c2 = theta[4]      # Constant in State 2
  phi11 = theta[5]   # First Lag in State 1
  phi12 = theta[6]   # First Lag in State 2
  phi21 = theta[7]   # Second Lag in State 1
  phi22 = theta[8]   # Second Lag in State 2
  beta11 = theta[9]  # Expectations in State 1   
  beta12 = theta[10] # Expectations in State 2  
  beta21 = theta[11] # Output Gap in State 1
  beta22 = theta[12] # Output Gap in state 2
  beta31 = theta[13] # REER in State 1
  beta32 = theta[14] # REER in State 2
  beta41 = theta[15] # ULC in State 1
  beta42 = theta[16] # ULC in state 2
  sigma = theta[17]  # Standard deviation (not state-dependent)
  
  # Check for invalid parameter values
  if (sigma <= 0) {
    print("Invalid sigma value detected")
    print(theta)
    return(Inf)
  }
  
  if (p11 < 0 || p11 > 1 || p22 < 0 || p22 > 1) {
    print("Invalid transition probabilities detected")
    print(theta)
    return(Inf)
  }
  
  n = length(series) # Number of observations
  trans_mat = matrix(c(p11, 1-p11, 1-p22, p22), nrow = 2, byrow = TRUE) # Transition probability matrix
  zeta = matrix(0, n, 2) # State matrix
  f = rep(0, n) # log likelihood
  eta = matrix(0, n, 2) # residual
  
  # Initial state probabilities
  zeta_Initial = c((1-trans_mat[2,2])/(2-trans_mat[1,1]-trans_mat[2,2]), 
                   (1-trans_mat[1,1])/(2-trans_mat[2,2]-trans_mat[1,1]))
  
  # Estimating the mean's and densities of the residuals
  for (i in 3:n) {
    # Evaluate the densities under both regimes
    control_vars = as.numeric(controls[i,])
    mu_1 = c1 + phi11*series[i-1] + phi21*series[i-2] + beta11*control_vars[1] + beta21*control_vars[2] + beta31*control_vars[3] + beta41*control_vars[4]
    mu_2 = c2 + phi12*series[i-1] + phi22*series[i-2] + beta12*control_vars[1] + beta22*control_vars[2] + beta32*control_vars[3] + beta42*control_vars[4]
    
    eta[i, 1] = dnorm(x = series[i], mean = mu_1, sd = sigma) + 1e-10
    eta[i, 2] = dnorm(x = series[i], mean = mu_2, sd = sigma) + 1e-10
    
    # Condition for bounds of density, cannot return zero or NaN
    if (eta[i, 1] == 0 || eta[i, 2] == 0) {
      print("Zero density detected")
      print(i)
      print(eta[i, ])
      return(Inf)
    }
    
    # Estimating the likelihood of the function
    if (i == 3) {
      f[i] = sum((trans_mat %*% c(eta[i,1], eta[i,2])) * zeta_Initial)
    } else {
      f[i] = sum((trans_mat %*% c(eta[i,1], eta[i,2])) * c(zeta[i-1,1], zeta[i-1,2]))
    }
    
    # Condition for negative values of the likelihood
    if (f[i] <= 0) {
      print("Non-positive f detected")
      print(i)
      print(f[i])
      return(Inf)
    }
    
    # Probabilities of states 1 or 2 at each point in time
    if (i == 3) {
      zeta[i,] = zeta_Initial
    } else {
      zeta[i,1] = (eta[i,1] * sum(trans_mat[,1] * c(zeta[i-1,1], zeta[i-1,2]))) / f[i]
      zeta[i,2] = (eta[i,2] * sum(trans_mat[,2] * c(zeta[i-1,1], zeta[i-1,2]))) / f[i]
    }
  }
  
  # Smooth probabilities 
  smooth_likelihood <- matrix(0, nrow = n, ncol = 2)
  smooth_prob <- matrix(0, nrow = n, ncol = 2)
  
  for(i in 3:n){
    smooth_likelihood[i, ] <- zeta[i, ] / sum(zeta[i, ])
  }
  
  smooth_prob[n, ] <- smooth_likelihood[n, ]
  
  for(i in (n-1):3){
    hk <- smooth_likelihood[i, ] %*% trans_mat
    for(j in 1:2){
      if(hk[j] > 1e-150){
        hk[j] <- smooth_prob[i+1, j] / hk[j]
      } else {
        hk[j] <- 0
      }
    }
    smooth_prob[i, ] <- smooth_likelihood[i, ] * (hk %*% trans_mat)
  }
  
  logf = sum(log(f[3:n]))
  if (is.nan(logf)) {
    logf = 99999999
    cat("Error: Not a number", "\n")
  }
  output = cbind(eta, f, zeta, smooth_prob)
  
  results <- list(-logf, zeta, smooth_prob, eta)
  
  return(-logf)
}
