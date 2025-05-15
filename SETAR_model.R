############################### FUNCTION TO RUN SETAR MODEL ################################################

#' @data Object is a data frame that contains the dependent (Y) and independent (X) variables
#' @y_name Dependent variable
#' @control_vars Independent variables that contain lags of Y and exogenous controls
#' @grid_search This is a boolean and the default is set to TRUE
#' @min_ratio Sets the lower bound for the grid search, default is 0.2
#' @max_ratio Sets the upper bound for the grid search, default is 0.8

SETAR <- function(data, y_name, control_vars, grid_search = T, min_ratio = 0.2, max_ratio = 0.8){

  # Preliminary check to determine if Y is in our data frame
  if (!y_name %in% colnames(data)){
    stop("Y column does not exist in the data")
  }
  
  # Check if Y is numeric
  if (!is.numeric(data[[y_name]])){
    data[[y_name]] <- as.numeric(data[[y_name]])
  }
  
  # Validate control variables
  if (!is.null(control_vars)) {
    missing_cols <- control_vars[!control_vars %in% colnames(data)]
    if (length(missing_cols) > 0) {
      stop(paste("The following control variables are not in the dataset:", 
                 paste(missing_cols, collapse = ", ")))
    }
  }
  
  # Create lags for Y 
  data <- data %>%
    mutate(L1 = lag(!!sym(y_name), 1),
           L2 = lag(!!sym(y_name), 2)) %>%
    na.omit()
  
  # Check if there are sufficient observations following lag creation
  if (nrow(data) < 10){
    stop("Not enough observations after creating lags")
  }
  
  # Create the dependent variable and the bounds for our thresholds
  Y <- data[[y_name]]
  min_ta <- round(min_ratio*length(Y), digits = 0)
  max_ta <- round(max_ratio*length(Y), digits = 0)
  
  # Prepare independent variables
  predictors <- c("L1", "L2")
  all_predictors <- predictors
  if (!is.null(control_vars)) {
    all_predictors <- c(predictors, control_vars)
  }
  
  # OLS Estimation 
  OLS_Estim <- function(X, Y){
    
    # Beta = (X'X)^(-1)X'Y
    B <- solve(t(X) %*% X) %*% t(X) %*% Y
    # Residuals
    epsilon <- Y - X %*% B  # Fixed this line
    # Degrees of freedom
    n <- nrow(X)
    p <- ncol(X)
    df <- n - p
    
    # Calculate the variance
    var <- sum(epsilon^2)/df
    # Covariance matrix
    cov <- var * solve(t(X) %*% X)
    # Standard errors
    se <- sqrt(diag(cov))
    # T-statistics
    t_stats <- B/se
    # P-values
    p_val <- 2 * (1 - pt(abs(t_stats), df = df))
    
    # Return the outputs from the estimation
    return(list(Beta = B,
                Std_errors = se,
                t_stats = t_stats,
                p_values = p_val,
                residuals = epsilon,
                df = df,
                X = X))
  }
  
  # Initialize results list
  results <- list()
  
  # Grid search for threshold values
  if (grid_search) {
    # Store residual sum of squares
    all_rss <- as.numeric()
    # Store all thresholds from grid search
    all_tau <- as.numeric()
    
    # Set the threshold variable, lag of our dependent variable
    tau_values <- sort(abs(data$L1))
    # Only search within the defined bounds
    search_indices <- min_ta:max_ta
    if (length(search_indices) > 0 && max(search_indices) <= length(tau_values)) {
      # Create the indicator functions
      for (i in search_indices) {
        tau <- tau_values[i]
        df <- data %>%
          mutate(d = if_else(L1 <= tau, 1, 0),
                 i.t = d * L1)
        
        # Prepare independent variable matrix with threshold variables
        model_vars <- c(all_predictors, "d", "i.t")
        
        # Check if all columns exist
        if (all(model_vars %in% colnames(df))) {
          # Prepare independent variable matrix
          X <- as.matrix(df[, model_vars])
          # Add a constant term
          X <- cbind(1, X)
          
          # Estimate the model using OLS
          model_results <- OLS_Estim(X, as.matrix(df[[y_name]]))
          
          # Store the residual sum of squares
          current_rss <- sum(model_results$residuals^2)
          # Store all residuals and thresholds
          all_rss <- c(all_rss, current_rss)
          all_tau <- c(all_tau, tau)
        }
      }
      
      if (length(all_rss) > 0) {
        # Find the threshold values with the lowest residual sum of squares
        # Take up to 100 best models (or fewer if we have fewer results)
        num_models <- min(100, length(all_rss))
        top_idx <- order(all_rss)[1:num_models]
        top_tau <- all_tau[top_idx]
        top_rss <- all_rss[top_idx]
        
        # Re-estimate the model using the best thresholds based on our grid search
        models <- list()
        for (i in seq_along(top_tau)) {
          tau <- top_tau[i]
          df <- data %>%
            mutate(d = if_else(L1 <= tau, 1, 0),
                   i.t = d * L1)
          
          # Prepare independent variable matrix
          if (!is.null(control_vars)) {
            X <- as.matrix(df[, c(predictors, "d", "i.t", control_vars)])
          } else {
            X <- as.matrix(df[, c(predictors, "d", "i.t")])
          }
          
          # Add intercept
          X <- cbind(1, X)
          
          # Estimate the model
          model_results <- OLS_Estim(X, as.matrix(df[[y_name]]))
          
          models[[i]] <- list(tau = tau,
                              model_results = model_results)
        }
        
        # Compute the mean threshold - divides data into two regimes
        mean_tau <- mean(top_tau, na.rm = TRUE)
        
        # Divide models into two inflation regimes, based on threshold
        # Data above the mean threshold (high regime)
        df_above <- data %>% 
          filter(L1 > mean_tau)
        # Data below the mean threshold (low regime)
        df_below <- data %>% 
          filter(L1 <= mean_tau)
        
        # Check if we have enough data in both regimes
        if (nrow(df_above) > 0 && nrow(df_below) > 0) {
          # Prepare the matrices for independent variables
          # High regime independent variables
          if (!is.null(control_vars)) {
            X_above <- as.matrix(df_above[, c(predictors, control_vars)])
          } else {
            X_above <- as.matrix(df_above[, predictors])
          }
          # Add constant
          X_above <- cbind(1, X_above)
          
          # Low regime independent variables - Fixed using df_below instead of df_above
          if (!is.null(control_vars)) {
            X_below <- as.matrix(df_below[, c(predictors, control_vars)])
          } else {
            X_below <- as.matrix(df_below[, predictors])
          }
          # Add constant
          X_below <- cbind(1, X_below)
          
          # Estimate High and Low inflation regime models
          High_regime <- OLS_Estim(X_above, as.matrix(df_above[[y_name]]))
          Low_regime <- OLS_Estim(X_below, as.matrix(df_below[[y_name]]))
          
          # List to store the submodel summaries
          submodels <- list(High = High_regime,
                            Low = Low_regime)
          
          # Store results
          results <- list(
            models = models,
            tau = top_tau,
            mean_tau = mean_tau,
            rss = top_rss,
            submodels = submodels
          )
        } else {
          warning("Not enough data in one or both regimes after threshold splitting")
          results <- list(
            models = models,
            tau = top_tau,
            mean_tau = mean_tau,
            rss = top_rss
          )
        }
      } else {
        warning("No valid threshold models were found")
      }
    } else {
      warning("Invalid search indices based on min_ratio and max_ratio")
    }
  }
  
  return(results)
}

      


  
