############################### FUNCTION TO RUN SETAR MODEL ################################################

#' @data Object is a data frame that contains the dependent (Y) and independent (X) variables
#' @y_name Dependent variable
#' @control_vars Independent variables that contain lags of Y and exogenous controls
#' @grid_search This is a boolean and the default is set to TRUE
#' @min_ratio Sets the lower bound for the grid search, default is 0.2
#' @max_ratio Sets the upper bound for the grid search, default is 0.8

SETAR <- function(data, y_name, control_vars, grid_search, min_ratio, max_ratio){

  # Preliminary check to determine if Y is in our data frame
  if (!y_name %in% colnames(data)){
    stop("Y column does not exist in the data")
    }

  # Check if Y is numeric
  if (!is.numeric(data[[y_name]]){
    data[[y_name]] <- as.numeric(data[[y_name]])
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
  if (!is.null(control_vars)) {
    predictors <- c(predictors, control_vars)
    }

  # OLS Estimation 
  OLS_Estim <- function(X, Y){
  
   # Beta = (X'X)^(-1)X'Y
   B <- solve(t(X) %*% X) %*% t(X) %*% Y
   # Residuals
   epsilon <- Y - B %*% X
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

  # Grid search for threshold values
  if (grid_search) {
    # Store residual sum of squares
    all_rss <- as.numeric()
    # Store all thresholds from grid search
    all_tau <- as.numeric()
    
    # Set the threshold variable, lag of our dependent variable
    tau <- sort(abs(data$L1))
    # Create the indicator functions
    for (tau in tau[min_ta:max_ta]) {
      df <- data %>%
        mutate(d = if_else(L1 <= tau, 1, 0),
               i.t = d * L1)

      # Prepare independent variable matrix
      X <- as.matrix(df[, c(predictors, "d", "i.t")])
      # Add a constant term to the independent variable matrix
      X <- cbind(1, X)

      # Estimate the model using OLS
      model_results <- OLS_Estim(X, as.matrix(df[[y_name]]))

      # Store the residual sum of squares
      current_rss <- sum(model_results$residuals^2)
      # Store all residuals and thresholds
      all_rss <- c(all_rss, current_rss)
      all_tau <- c(all_tau, tau)
    }

    # Find the threshold values with the lowest residual sum of squares
    top_idx <- order(all_rss)[1:100]
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
      X <- as.matrix(df[, c(predictors, control_vars)])
      # Add intercept
      X <- cbind(1, X)

      # Estimate the model
      model_results <- OLS_Estim(X, as.matrix(df[[y_name]]))

      models[i] <- list(tau = tau,
                        model_results = model_results)
    }

    # Compute the mean threshold - divides data into two reimes
    mean_tau <- mean(all_tau, na.rm = TRUE)

    # Divide models into two inflation regimes, based on threshold
    # Data above the mean threshold (high regime)
    df_above <- data %>% 
      filter(L1 > mean_tau)
    # Data below the mean threshold (low regime)
    df_below <- data %>% 
      filter(L1 <= mean_tau)
    # Prepare the matrices for independent variables
    # High regime independent variables
    X_above <- as.matrix(df_above[, predictors])
    # Add constant
    X_above <- cbind(1, X)
    # Low regime independent variables
    X_below <- as.matrix(df_above[, predictors])
    # Add constant
    X_below <- cbind(1, X)

    # Estimate High and Low inflation regime models
    High_regime <- OLS_Estim(X_above, as.matrix(df_above[[y_name]]))
    Low_regime <- OLS_Estim(X_below, as.matrix(df_below[[y_name]]))

    # List to store the submodel summaries
    submodels <- list(High = High_regime,
                      Low = Low_regime)

    # Return a list of all the outputs from the function
    results <- list(models = models,
                    tau = top_thresholds,
                    mean_tau = mean_tau,
                    rss = top_rss,
                    submodels = submodels)
  }
}

      


  
