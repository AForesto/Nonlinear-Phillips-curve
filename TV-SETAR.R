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
           L2 = lag(!!rlang::sym(y_name), 2),
           Avg_inflation = zoo::rollmeanr(!!rlang::sym(y_name), k = 8, fill = NA)) %>%
    na.omit()
  
  # Check data after preparation
  if (nrow(data) <= start_size) {
    stop("Not enough data after preparation")
  }
  
  # Validate control variables if provided
  if (!is.null(control_vars)) {
    missing_cols <- control_vars[!control_vars %in% colnames(data)]
    if (length(missing_cols) > 0) {
      stop(paste("The following control variables are not in the dataset:", 
                 paste(missing_cols, collapse = ", ")))
    }
  }
  
  # Storage for results
  results <- list()
  
  # Matrix-based model estimation function
  OLS_Estim <- function(X, Y, crit = "BIC") {
    # Check for sufficient observations
    if (nrow(X) <= ncol(X)) {
      return(NULL)  # Not enough observations for estimation
    }
    
    # Beta = (X'X)^(-1)X'Y
    B <- tryCatch({
      solve(t(X) %*% X) %*% t(X) %*% Y
    }, error = function(e) {
      # Handle singular matrix errors
      return(NULL)
    })
    
    # If matrix inversion failed, return NULL
    if (is.null(B)) {
      return(NULL)
    }
    
    # Compute residuals
    epsilon <- Y - X %*% B
    
    # Compute degrees of freedom
    n <- nrow(X)
    p <- ncol(X)
    df <- n - p
    
    # Compute variance estimate
    var <- sum(epsilon^2) / df
    
    # Compute covariance matrix of coefficients
    cov <- tryCatch({
      var * solve(t(X) %*% X)
    }, error = function(e) {
      # Handle singular matrix errors
      return(NULL)
    })
    
    # If covariance calculation failed, return limited results
    if (is.null(cov)) {
      return(list(
        coefficients = B,
        residuals = epsilon,
        df = df,
        X = X,
        rss = sum(epsilon^2),
        fitted_values = X %*% B
      ))
    }
    
    # Compute standard errors
    se <- sqrt(diag(cov))
    
    # Compute t-statistics
    t_stats <- B / se
    
    # Compute p-values (with NaN handling)
    p_values <- tryCatch({
      2 * (1 - pt(abs(t_stats), df = df))
    }, warning = function(w) {
      # Handle NaN warnings
      rep(NA, length(t_stats))
    })
    
    # Compute information criteria
    rss <- sum(epsilon^2)
    aic <- n * log(rss/n) + 2 * p
    bic <- n * log(rss/n) + log(n) * p
    
    # Select criterion
    if (crit == "AIC") {
      selected_criterion <- aic
    } else if (crit == "BIC") {
      selected_criterion <- bic
    } else {
      selected_criterion <- bic  # Default to BIC if invalid criterion
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
    
    # Ensure we have a valid range of indices
    if (min_ta >= max_ta || min_ta < 1 || max_ta > length(subset_data$L1)) {
      next
    }
    
    tau_values <- sort(abs(subset_data$L1))[min_ta:max_ta]
    
    # Skip if we don't have enough threshold values
    if (length(tau_values) < 1) {
      next
    }
    
    # Initialize the best model with high initial criterion value
    best_model <- list(criterion = Inf)
    best_tau <- NA
    
    # Compute threshold statistics
    mean_tau <- mean(tau_values, na.rm = TRUE)
    median_tau <- median(tau_values, na.rm = TRUE)
    lower_quartile <- quantile(tau_values, 0.25, na.rm = TRUE)
    upper_quartile <- quantile(tau_values, 0.75, na.rm = TRUE)
    min_tau <- min(tau_values, na.rm = TRUE)
    max_tau <- max(tau_values, na.rm = TRUE)
    
    # Loop through possible thresholds
    for (tau in tau_values) {
      # Create threshold indicator
      df <- subset_data %>%
        mutate(d = if_else(L1 <= tau, 1, 0),
               i.t = d * L1)
      
      # Prepare predictors
      predictors <- c("L1", "L2")
      if (!is.null(control_vars)) {
        # Make sure control variables exist in the dataset
        valid_controls <- control_vars[control_vars %in% colnames(df)]
        if (length(valid_controls) > 0) {
          predictors <- c(predictors, valid_controls)
        }
      }
      
      # Check if all columns exist in the dataframe
      all_vars <- c(predictors, "d", "i.t")
      if (!all(all_vars %in% colnames(df))) {
        next  # Skip this iteration if any variable is missing
      }
      
      # Prepare design matrix
      X <- as.matrix(df[, all_vars])
      X <- cbind(1, X)  # Add intercept
      
      # Estimate model using matrix algebra
      model_results <- OLS_Estim(X, as.matrix(df[[y_name]]), criterion)
      
      # Skip if estimation failed
      if (is.null(model_results)) {
        next
      }
      
      # Check if the current model is better based on the selected criterion
      if (!is.null(model_results$criterion) && 
          !is.na(model_results$criterion) && 
          model_results$criterion < best_model$criterion) {
        # Convert matrix results to list format
        best_tau <- tau
        best_model <- list(
          tau = tau,
          criterion = model_results$criterion,
          rss = model_results$rss,
          aic = model_results$aic,
          bic = model_results$bic,
          coefficients = as.list(model_results$coefficients),
          std_errors = if (!is.null(model_results$std_errors)) as.list(model_results$std_errors) else NULL,
          t_stats = if (!is.null(model_results$t_stats)) as.list(model_results$t_stats) else NULL,
          p_values = if (!is.null(model_results$p_values)) as.list(model_results$p_values) else NULL,
          fitted_values = as.vector(model_results$fitted_values)
        )
      }
    }
    
    # If we found a good model, store it with threshold statistics
    if (!is.infinite(best_model$criterion)) {
      # Calculate threshold dummies based on best threshold
      mean_threshold <- mean_tau  # Using mean threshold if no best model was found
      if (!is.na(best_tau)) {
        mean_threshold <- best_tau
      }
      
      # Create dummies for threshold
      dummy_avg <- if_else(subset_data$Avg_inflation > mean_threshold, 1, 0)  
      dummy_L1 <- if_else(subset_data$L1 > mean_threshold, 1, 0)
      
      # Add threshold statistics and dummy variables to the best model
      best_model$mean_tau <- mean_tau
      best_model$median_tau <- median_tau
      best_model$lower_quartile <- lower_quartile
      best_model$upper_quartile <- upper_quartile
      best_model$min_tau <- min_tau
      best_model$max_tau <- max_tau
      best_model$dummy_avg <- dummy_avg
      best_model$dummy_L1 <- dummy_L1
      best_model$date <- end
      best_model$n_obs <- nrow(subset_data)
      best_model$Avg_inflation <- subset_data$Avg_inflation
      
      # Store the best model from this iteration
      results[[length(results) + 1]] <- best_model
    }
  }
  
  return(results)
}
