################### FUNCTION THAT DETERMINES OPTIMAL LAMBDA IN SMOOTHING - B&B (2019) #########################

#' @object Object created in the Local_Proj function, a list of outputs.
#' @size Size of the cross validation sample  
#' @trace Set to 1 to make the output more readable in R's console.

Cross_Validation <- function(object, size, trace = 1){
  
  T <- object$T
  Length <- length(object$lambda)
  
  index <- ceiling((object$IDX[,1]/T)* size)
  RSS <- rep(0, Length)
  
  for(l in 1:Length){
    if (trace == 1){
      cat('.')
    }
    RSS.l <- rep(0, size)
    
    for (i in 1:size){
      Y.in <- object$Y[index != i]
      X.in <- object$X[index != i,]
      Y.out <- object$Y[index == i]
      X.out <- object$X[index == i,]
      
      A <- t(X.in)%*%X.in + object$lambda[l] * object$TS * ((size-1)/size) * object$P
      b <- t(X.in)%*%Y.in
      beta <- Matrix::solve(A,b)
      RSS.l[i] <- mean((Y.out - X.out %*% beta)**2) 
    }
    RSS[l] <- mean(RSS.l)
  }
  if (trace == 1){
    cat('\n')
  }
  object$RSS <- RSS
  object$IDX.opt <- tail(which(min(RSS)==RSS), 1)
  object$impulse.response.opt <- object$impulse.response[, tail(which(min(RSS)==RSS), 1)]
  
  object
}
