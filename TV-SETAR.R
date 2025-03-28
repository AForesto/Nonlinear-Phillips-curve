######################################### FUNCTION TO ESTIMATE A TIME-VARYING SETAR MODEL ########################################

#' @data Object that is typically a data frame and contains the dependent (Y) and independent (X) variables
#' @y_name String that indicates the dependent variables
#' @control_vars List of strings that contains our independent variables, lags and exogenous controls. Default is set to NULL
#' @start_size Numeric value that indicates what the starting sample size is for the intitial regression
#' @min_ratio Lower bound for the grid search, default is 0.2
#' @max_ratio Upper bound for the grid search, default is 0.8

TV-SETAR <- function (data, y_name, control_vars, start_size, min_ratio = 0.2, max_ratio = 0.8){

  # OLS Estimation
  OLS_Estim <- function (X, Y){
    # 
    
