############################ FUNCTION TO ESTIMATE HETEROSKEDASTICITY-CONSISTENT STANDARD ERRORS #######################

#' @model Object that contains the model diagnostics - can be pulled from the lm function or other model objects
#' @data  Data frame or matrix that contains both the dependent (Y) and independent variables (X)

HC1_Estimation <- function(model, data){

  # Convert the model object to matrix form
  X <- model.matrix(model)
  # Extract the residuals of the model
  epsilon <- residuals(model) 
  # Extract the fitted values form the model
  H <- hatvalues(model)

  # The HC1 adjustment factor
  D <- diag((1-H)^2)

  # Estimate the robust covariance matrix
  covariance_matrix <- solve(t(X) %*% X) %*% (t(X) %*% D %*% X) %*% solve(t(X) %*% X)
  # Estimate robust standard errors
  standard_errors_robust <- sqrt(diag(covariance_matrix))

  # Estimate t-statistics and p-values
  # model coefficients
  coefficients <- coef(model)
  # t-values
  t_values_robust <- coefficients/standard_errors_robust
  # p-values
  p_values_robust <- 2*pt(-abs(t_values_robust), df = df.residual(model))

  # Return a list of the outputs
  return(list(standard_errors = standard_errors_robust,
              t_values        = t_values_robust,
              p_values        = p_values_robust))
  
