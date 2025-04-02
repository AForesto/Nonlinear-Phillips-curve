######################################### FUNCTION TO ESTIMATE A TIME-VARYING SETAR MODEL ########################################

#' @data Object that is typically a data frame and contains the dependent (Y) and independent (X) variables
#' @y_name String that indicates the dependent variables
#' @control_vars List of strings that contains our independent variables, lags and exogenous controls. Default is set to NULL
#' @start_size Numeric value that indicates what the starting sample size is for the intitial regression. Default is set to 8 quarters.
#' @min_ratio Lower bound for the grid search, default is 0.2
#' @max_ratio Upper bound for the grid search, default is 0.8

TV-SETAR <- function (data, y_name, control_vars, start_size = 8, min_ratio = 0.2, max_ratio = 0.8){

   # Initial data preparation
  data <- data %>%
    mutate(L1 = lag(!!rlang::sym(y_name), 1),
           L2 = lag(!!rlang::sym(y_name), 2)) %>%
    na.omit()
  
  # Check data after preparation
  if (nrow(data) <= start_size) {
    stop("Not enough data after preparation")
  }
  
  # Storage for results
  results <- list()
  
  # Matrix-based model estimation function
  OLS_Estim <- function(X, Y, criterion = "BIC") {
    # Beta = (X'X)^(-1)X'Y
    B <- solve(t(X) %*% X) %*% t(X) %*% Y
    
    # Compute residuals
    epsilon <- Y - X %*% B
    
    # Compute degrees of freedom
    n <- nrow(X)
    p <- ncol(X)
    df <- n - p
    
    # Compute variance estimate
    var <- sum(epsilon^2) / df
    
    # Compute covariance matrix of coefficients
    cov <- var * solve(t(X) %*% X)
    
    # Compute standard errors
    se <- sqrt(diag(cov))
    
    # Compute t-statistics
    t_stats <- B / se
    
    # Compute p-values
    p_values <- 2 * (1 - pt(abs(t_stats), df = df))
    
    # Compute information criteria
    rss <- sum(epsilon^2)
    aic <- n * log(rss/n) + 2 * p
    bic <- n * log(rss/n) + log(n) * p
    
    # Select criterion
    if (criterion == "AIC") {
      selected_criterion <- aic
    } else if (criterion == "BIC") {
      selected_criterion <- bic
    } else {
      stop("Invalid criterion specified. Use 'AIC' or 'BIC'.")
    }
    
    return(list(
      coefficients = B,
      std_errors = se,
      t_stats = t_stats,
      p_values = p_values,
      residuals = epsilon,
      df = df,
      X = X,
      rss = rss,
      aic = aic,
      bic = bic,
      criterion = selected_criterion,
      fitted_values = X %*% B
    ))
  }
   
  # Recursive expanding window estimation
  for (end in (start_size + 1):nrow(data)) {
    # Subset data from start to the current endpoint
    subset_data <- data[1:end, ]
    Y <- subset_data[[y_name]]
    min_ta <- round(min_ratio * length(Y), 0)
    max_ta <- round(max_ratio * length(Y), 0)
    
    # Check indices
    if (min_ta > length(subset_data$L1) || max_ta > length(subset_data$L1)) {
      next
    }
    
    tau <- sort(abs(subset_data$L1))[min_ta:max_ta]
    
    # Initialize the best model with high initial criterion value
    best_model <- list(criterion = Inf)
    
    # Compute threshold statistics
    mean_tau <- mean(tau)
    median_tau <- median(tau)
    lower_quartile <- quantile(tau, 0.25)
    upper_quartile <- quantile(tau, 0.75)
    min_tau <- min(tau)
    max_tau <- max(tau)
    
    # Loop through possible thresholds
    for (tau in tau) {
      # Create threshold indicator
      df <- subset_data %>%
        mutate(d = if_else(L1 <= tau, 1, 0),
               i.t = d * L1)
      
      # Prepare predictors
      predictors <- c("L1", "L2")
      if (!is.null(control_vars)) {
        predictors <- c(predictors, control_vars)
      }
      
      # Prepare design matrix
      X <- as.matrix(df[, c(predictors, "d", "i.t")])
      X <- cbind(1, X)  # Add intercept
      
      # Estimate model using matrix algebra
      model_results <- OLS_Estim(X, as.matrix(df[[y_name]]), criterion)
      
      # Check if the current model is better based on the selected criterion
      if (model_results$criterion < best_model$criterion) {
        # Convert matrix results to list format
        best_model <- list(
          tau = tau,
          criterion = model_results$criterion,
          rss = model_results$rss,
          aic = model_results$aic,
          bic = model_results$bic,
          coefficients = as.list(model_results$coefficients),
          t_stats = as.list(model_results$t_stats),
          p_values = as.list(model_results$p_values),
          std_errors = as.list(model_results$std_errors),
          fitted_values = as.vector(model_results$fitted_values)
        )
      }
    }
    
    # Add threshold statistics to the best model
    best_model$mean_tau <- mean_tau
    best_model$median_tau <- median_tau
    best_model$lower_quartile <- lower_quartile
    best_model$upper_quartile <- upper_quartile
    best_model$min_tau <- min_tau
    best_model$max_tau <- max_tau
    
    # Store the best model from this iteration
    results[[length(results) + 1]] <- best_model
  }
  
  return(results)
}
