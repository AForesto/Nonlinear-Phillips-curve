############### FUNCTION THAT ESTIMATES THE CONFIDENCE INTERVALS OF OUR LP ######################

#' @object The object returned following our Local_Proj and Cross_Validation functions
#' @l The index of the theta matrix with the optimal estimates

Confidence_Int <- function(object, l=1){
  
  u <- object$Y - object$X %*% object$coeffs[,l];
  S <- object$X * (u %*% t(rep(1, ncol(object$X))))
  
  # OUTER
  outer <- solve(t(object$X) %*% object$X + object$lambda[l]*object$TS*object$P)
  
  # INNER
  num_lag = min(floor(1.2*(object$T)**(1/3)), object$T)
  num_lag = object$Horizon
  weight = (num_lag + 1 -(0:num_lag))/(num_lag + 1)
  V <- t(S) %*% S
  for (i in 1:num_lag){
    gamma_i <- t(S[(i + 1):object$T,]) %*% S[1:(object$T - i),]
    g.plus.g.prime <- gamma_i + t(gamma_i)
    V <- V + weight[i + 1] * g.plus.g.prime
  }
  
  inner <- V
  
  V <- outer %*% inner %*% outer
  
  if (object$type == "regular"){
    standard_error <- sqrt(diag(V[1:(object$Horizon + 1 - object$h1), 1:(object$Horizon + 1 - object$h1)]))
    confidence <- matrix(0, length(standard_error), 2)
    confidence[, 1] <- object$optim[, l] + standard_error*qnorm(0.05)
    confidence[, 2] <- object$optim[, l] + standard_error*qnorm(0.95)
  }
  else{
    V <- as.matrix(object$basis) %*% V[1:object$XS, 1:object$XS] %*% t(as.matrix(object$basis))
    standard_error <- sqrt(diag(V))
    confidence <- matrix(0, length(standard_error), 2)
    confidence[, 1] <- object$optim[, l] + standard_error*qnorm(0.05)
    confidence[, 2] <- object$optim[, l] + standard_error*qnorm(0.95)
  }
  
  impulse.response.conf <- matrix(NA, object$Horizon + 1, 2)
  impulse.response.conf[(1 + object$h1):(object$Horizon + 1), ] <- confidence*object$delta
  
  object$standard_error <- standard_error
  object$impulse.response.conf <- impulse.response.conf
  
  object
}
