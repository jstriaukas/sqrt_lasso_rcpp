rsqrtlasso <- function(X,Y,standardize=FALSE,Lambda=NULL,lambda.choice="bc",maxIter=10000, OptTolNorm=1e-6, OptTolObj=1e-8) {
  if (is.null(Lambda)){ # code take from RPtests and scalreg packages // ADD: AB choice.
    if (p == 1) {
      L = 0.5
    }
    else if(lambda.choice=="sh") {
      L <- 0.1
      Lold <- 0
      while (abs(L - Lold) > 0.001) {
        k <- (L^4 + 2 * L^2)
        Lold <- L
        L <- -qnorm(min(k/p, 0.99))
        L <- (L + Lold)/2
      }
      Lambda = sqrt(2/n) * L
    } else if (lambda.choice=="bc"){
      c <- 1.1
      alpha <- 0.05
      Lambda <- c*sqrt(nrow(X))*qnorm(1-(alpha/(2*ncol(X))))
    }
  }
  beta.min <- sqrtlasso(X,Y,Lambda,maxIter=maxIter,OptTolNorm=OptTolNorm,OptTolObj=OptTolObj,standardize=standardize,as.numeric(ncol(X)),1)
  obj <- list(beta.min=beta.min,lambda.min=Lambda)
  return(obj)
}