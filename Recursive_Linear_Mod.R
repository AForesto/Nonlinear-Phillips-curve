############################# FUNCTION TO ESTIMATE A TIME-VARYING LINEAR MODEL ################################

#' @data a matrix or dataframe containing both the dependent and independent variables
#' @y_name a string that indicates the name of the dependent variable in our dataframe/matrix
#' @controls a list of strings indicating the independent variables
#' @start_size takes a numeric value that indicates the how many time steps are used in the initial regression

Recursive_Linear_Model <- function(data, y_name, controls = NULL, start_size) {

  # Check the data after creating lags
  print(paste("There are", nrow(data), "rows after preparation"))
  if (nrow(data) <= start_size) {
    stop("Not enough rows in our data set after preparation")
  }
  
  # Initialize storage of ratios in a vector with proper length
  ratios <- numeric(nrow(data))
  coef_L1_values <- numeric(nrow(data))
  coef_control_values <- numeric(nrow(data))
  
  # Handle controls parameter - can be single string or character vector
  if (is.null(controls)) {
    stop("Controls parameter cannot be NULL")
  }
  
  # Ensure controls is treated as a character vector
  if (!is.character(controls)) {
    stop("Controls must be a character string or vector of column names")
  }
  
  # Convert single string to vector if needed
  if (length(controls) == 1) {
    controls <- c(controls)
  }
  
  # OLS Estimation function
  OLS_Estim <- function(X, Y) {
    # Add intercept
    X <- cbind(1, X)
    
    # Check for potential issues before calculations
    if (ncol(X) <= 1 || nrow(X) <= ncol(X)) {
      return(list(
        Beta = rep(NA, ncol(X)),
        Std_errors = rep(NA, ncol(X)),
        t_stats = rep(NA, ncol(X)),
        p_values = rep(NA, ncol(X)),
        residuals = rep(NA, length(Y)),
        df = NA
      ))
    }
    
    # Try-catch to handle potential numerical issues
    tryCatch({
      # Beta = (X'X)^(-1)X'Y
      XtX <- t(X) %*% X
      if (det(XtX) < 1e-10) {
        # Near-singular matrix, use regularization or pseudo-inverse
        B <- MASS::ginv(XtX) %*% t(X) %*% Y
      } else {
        B <- solve(XtX) %*% t(X) %*% Y
      }
      
      # Residuals
      epsilon <- Y - X %*% B
      
      # Degrees of freedom
      n <- nrow(X)
      p <- ncol(X)
      df <- n - p
      
      if (df <= 0) {
        return(list(
          Beta = B,
          Std_errors = rep(NA, length(B)),
          t_stats = rep(NA, length(B)),
          p_values = rep(NA, length(B)),
          residuals = epsilon,
          df = df
        ))
      }
      
      # Calculate the variance
      var <- sum(epsilon^2) / df
      
      # Covariance matrix
      cov <- var * solve(XtX)
      
      # Standard errors
      se <- sqrt(diag(cov))
      
      # T-statistics
      t_stats <- B / se
      
      # P-values - safely calculated
      p_val <- rep(NA, length(t_stats))
      for (i in 1:length(t_stats)) {
        if (!is.na(t_stats[i]) && is.finite(t_stats[i]) && !is.na(df) && df > 0) {
          p_val[i] <- 2 * (1 - pt(abs(t_stats[i]), df = df))
        }
      }
      
      return(list(
        Beta = B,
        Std_errors = se,
        t_stats = t_stats,
        p_values = p_val,
        residuals = epsilon,
        df = df
      ))
    }, error = function(e) {
      # Return NA values if calculation fails
      return(list(
        Beta = rep(NA, ncol(X)),
        Std_errors = rep(NA, ncol(X)),
        t_stats = rep(NA, ncol(X)),
        p_values = rep(NA, ncol(X)),
        residuals = rep(NA, length(Y)),
        df = NA,
        error = e$message
      ))
    })
  }
  
  # For loop to iterate through linear model recursively 
  for (end in (start_size + 1):nrow(data)) {
    subset_data <- data[1:end, ]
    Y <- subset_data[[y_name]]
    
    # Create design matrix from the controls
    if (length(controls) == 1) {
      X <- matrix(subset_data[[controls]], ncol = 1)
    } else {
      X <- as.matrix(subset_data[, controls, drop = FALSE])
    }
    
    # Check if there's enough data
    if (nrow(X) <= ncol(X) + 1) {
      # Not enough data for this iteration, skip
      next
    }
    
    # Run OLS estimation
    result <- OLS_Estim(X, Y)
    
    # Extract coefficients - need to account for potential errors
    if (!is.null(result$Beta) && length(result$Beta) >= 2) {
      coef_L1 <- result$Beta[2]  # First control variable (after intercept)
      
      # Only try to get second coefficient if it exists
      if (length(result$Beta) >= 3) {
        coef_control <- result$Beta[3]  # Second control variable (if exists)
      } else {
        coef_control <- NA
      }
      
      # Store the coefficients
      coef_L1_values[end] <- coef_L1
      coef_control_values[end] <- coef_control
      
      # Compute the ratio of control/lag
      if (!is.na(coef_L1) && !is.na(coef_control) && abs(coef_L1) > 1e-10) {
        ratio <- coef_control / coef_L1
        ratios[end] <- ratio
      }
    }
  }
  
  # Remove the empty entries (for indices below start_size)
  valid_indices <- (start_size + 1):length(ratios)
  ratios <- ratios[valid_indices]
  coef_L1_values <- coef_L1_values[valid_indices]
  coef_control_values <- coef_control_values[valid_indices]
  
  # Check if we have a single control or multiple controls
  if (length(controls) == 1) {
    # Only one control variable - no ratio to calculate
    return(list(
      B1 = coef_L1_values,
      single_control = TRUE
    ))
  } else {
    # Multiple control variables - can calculate ratios
    return(list(
      ratios = ratios,
      B1 = coef_L1_values,
      B2 = coef_control_values,
      single_control = FALSE
    ))
  }
}

