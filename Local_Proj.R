####### FUNCTION THAT ESTIMATES SMOOTH LOCAL PROJECTION OF BARNICHON AND BRONWLEES (2019) ########

#' @y Dependent variable, stored as a vector 
#' @x Shock variable, stored as a vector
#' @z Control or independent variables that includes the shock variable, stored as a matrix,
#' @const This is a boolean and the default is set to TRUE
#' @type Determines the kind of local projection, options are "regular" or "smooth". Default is "regular"
#' @Horizon The number of h-steps in the local projection, takes on a numeric value
#' @h1 The initial point of the shock or impulse
#' @degree This is the order of the polynomial function used in the b-spline
#' @zero Set to true to exclude contemproaneous response
#' @lambda Used on the smooth local projection, numeric vector for the shrinking parameter
#' @trace Set to 1, to make the outputs more readable 
#' @shock_scale Applies a factor to the shock to scale it. e.g. based on weight of item in CPI basket

Local_Proj <- function(y, x, z = NULL, const = TRUE, type = "regular", Horizon, h1, degree = 0,
                       zero = FALSE, lambda = 0, trace = 1, shock_scale = 1){
  # Dimensions of the data are the same length as the dependent variable
  T = length(y)
  
  # Create the basis function for the smooth projection, based on spline interpolation
  if(type == "smooth"){
    horizon.range <- h1:Horizon # Horizon range from initial impact until period end
    basis.deg <- 3 # Degree of the polynomial basis function, cubic
    knots <- seq(-basis.deg + h1, Horizon + basis.deg, 1) # Points for the spline 
    basis.function <- spline.des(knots, horizon.range, basis.deg + 1, 0 * horizon.range, outer.ok = T)$design
  }
  
  # Create control matrix
  if(const == TRUE){
    z <- cbind(rep(1,T), z) # Row length the same as dependent var, columns same as controls chosen, add a constant term
  }
  
  # Shock size is determined based on regression of shock on control variables, else it is the std.dev of the shock
  if(!is.null(z)){
    delta <- summary(lm(x ~ lag(x)))$sigma # Shock size is variance based on the regression, can use controls
  }
  else{
    delta <- sd(x) # Or shock size is the standard deviation of shock variable
  }
  
  delta <- delta * shock_scale # Scale the shock if needed
  
  # Dimensions for impulse responses
  HR <- H + 1 - h1 # Number of periods projection spans
  TS <- T*HR # Total number of time steps which is the dimension multiplied by the horizon
  
  # Size of the XS matrix based on the projection type
  if (type == "regular"){
    XS <- HR # time steps of X are equivalent to Horizon if regular
  }
  else{
    XS <- ncol(basis.function) # For smooth projections the time steps are a function of the number of columns in the spline
  }
  ZS <- HR # Steps for the control variables are equivalent ot the horizon
  if(!is.null(z) | const == T){
    NZ <- ncol(z)
  }
  else{
    NZ <- 0
  }
  
  # Create matrices for the dependent variable, independent variable and control variables
  IDX <- matrix(0, TS, 2) # Index matrix to store time steps and horizons
  Y <- rep(NA, TS) # Shifted response vector used in Projections
  X.b <- matrix(0, TS, XS) # Matrix for the shock variable
  X.c <- array(0, dim = c(TS, HR, NZ)) # Array for the control variables
  
  # Construst dependent and independent for each point in the forecast
  for (t in 1:(T-h1)){
    idx.start <- (t-1)*HR + 1 # Start of the index for the forecast
    idx.end <- t*HR # End point for the index
    
    IDX[idx.start:idx.end, 1] <- t # Time index
    IDX[idx.start:idx.end, 2] <- h1:Horizon # Horizon index
    
    # y, dependent variable
    y.range <- (t+h1):min((t+Horizon), T) # Range of values for the dependent variable
    Y[idx.start:idx.end] <- c(y[y.range], rep(NA, HR - length(y.range))) # Full the range with NA's
    
    # x, shock variable
    if(type == "regular"){
      X.b[idx.start:idx.end, ] <- diag(HR) * x[t] # Shock at time t for each forecast step in the horizon if regular IRF is used
    }
    else{
      X.b[idx.start:idx.end, ] <- basis.function * x[t] # For the smooth projection we use the basis function
    }
    
    # z, control variables
    for (i in seq_len(NZ)){
      X.c[idx.start:idx.end, ,i] <- diag(HR) * z[t, i] # Controls for each forecast step in the horizon
    }
  }
  
  # Combine the shock variable (X.b) and the control variables (X.c) into a single matrix
  X <- cbind(X.b)
  for(i in seq_len(NZ)){
    X <- cbind(X, X.c[,,i]) # Add each control variable to main matrix X
  }
  X <- Matrix::Matrix(X, sparse = T) # Convert to a sparse matrix to optimise memory usage
  
  # Filter out rows where dependet variable (Y) is NA
  sel <- !is.na(Y)
  IDX <- IDX[sel,]
  Y  <- Y[sel]
  X  <- X[sel,]
  TS <- length(Y) # Update the time steps after omitting missing or NA values
  
  # Product matrices for the regression
  XX <- t(X)%*%X # Product of X'X
  XY <- t(X)%*%Y # Product of X'Y
  
  # Shrinkage or penalty matrix used for regularisation
  P <- matrix(0, ncol(X), ncol(X)) #Penalty Matrix
  P <- Matrix(P, sparse = T)
  
  # Construct the penalty matrix for smooth projections
  if(type=="smooth"){
    D <- diag(XS) # Identity matrix for shrinkage
    for(k in seq_len(degree)) D <- diff(D) # Differencing for smoothing
    
    if(zero){
      DP <- rep(0, XS)
      DP[XS] <- 1 # Penalty for the last parameter if needed
      D <- rbind(D, DP)
    }
    P[1:XS, 1:XS] <- t(D) %*% D # Updating the penalty matrix at each time step
  }
  
  # Matrices to store the impulse responses and coefficients of the regression
  impulse.response <- matrix(0, Horizon + 1, length(lambda)) # IRF's
  coeffs <- matrix(0, ncol(X), length(lambda)) # Coefficients
  optim <- matrix(0, HR, length(lambda)) # Optimised values based on lambda

  # Loop over each shrinkage parameter (lambda) to compute the coefficients and IRF's
  for(i in 1:length(lambda)){
    if (trace==1){
      cat('.')
    }
    
    # Solve for the coefficients (X'X + lambda*p)^(-1)
    A <- XX + lambda[i]*TS*P # Regularised X'X matrix
    b <- XY # Regularised X'Y matrix
    coeffs[, i] <- as.vector(Matrix::solve(A, b)) #Solve for the coefficients
    
    # Extract the coefficients of the shock
    if (type == "regular"){
      optim[, i] <- coeffs[1:XS, i] # For regular IRF's, use the first XS coefficients
    }
    else{
      beta <- coeffs[1:XS, i] # For smooth IRF's apply the basis function
      optim[, i] <- as.matrix(basis.function) %*% as.vector(beta)
    }
    
    # Compute the IRF's by multiplying by the shock size
    impulse.response[(h1+1):(Horizon+1), i] <- optim[, i] * delta
    
    Y_h <- X %*% coeffs[, i] # Predicted y based on coefficients
    block <- nrow(Y_h)/(Horizon + 1)
    
    Y_h <- colMeans(matrix(Y_h, block, Horizon + 1))

  }
  if (trace == 1){
    cat('\n')
  }

  # Create list to store results
  objects <- list()
  objects$type <- type # Store the type of projection used
  
  if(type == "smooth"){
    objects$basis <- basis.function # If smooth IRF used store the basis function
  }
  
  # Store all the outputs
  objects$h1               <- h1
  objects$Horizon          <- Horizon
  objects$XS               <- XS
  objects$HR               <- HR
  objects$T                <- T
  objects$Horizon          <- Horizon
  objects$TS               <- TS
  objects$IDX              <- IDX
  objects$y                <- y 
  objects$x                <- x
  objects$z                <- z
  objects$Y                <- Y
  objects$X                <- X
  objects$coeffs           <- coeffs
  objects$optim            <- optim
  objects$lambda           <- lambda
  objects$P                <- P
  objects$impulse.response <- impulse.response
  objects$delta            <- delta
  objects$forecast_y       <- Y_h
  
  objects
}
