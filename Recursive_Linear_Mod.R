############################# FUNCTION TO ESTIMATE A TIME-VARYING LINEAR MODEL ################################

#' @data a matrix or dataframe containing both the dependent and independent variables
#' @y_name a string that indicates the name of the dependent variable in our dataframe/matrix
#' @controls a list of strings indicating the independent variables
#' @start_size takes a numeric value that indicates the how many time steps are used in the initial regression

Recursive_Linear_Model <- function(data, y_name, controls = NULL, start_size) {

  # Create lags of our dependent variable
  data <- data %>%
    mutate(L1 = lag(!!rlang::sym(y_name), 1)) %>%
    na.omit()

 # Check the data after creating lags
 print(paste("There are ", nrow(data), " after preparation"))
 if (nrow(data) <= start_size) {
   stop("Not enough rows in our data set after preparation")
 }

 # Initialise storage of ratios in a vector
 ratios <- as.numeric()

 # For loop to iterate through linear model recursively 
 for (end in (start_size + 1):nrow(data)) {
   subset_data <- data[1:end,]
   Y <- subset_data[[y_name]]
   Lags <- subset_data[["L1"]]
   X <- c(Lags, controls)

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

coef_L1 <- Beta[,1]
coef_control <- Beta[,2]

# Compute the ratio of control/lag 
if (!is.na(coef_L1) && !is.na(coef_control)) {
  ratio <- coef_control/coef_L1
  ratios[[end]] <- c(ratio)
  } 
 }

object <- list(
ratios <- ratios
B1 <- coef_L1
B2 <- coef_control
)

object
}
