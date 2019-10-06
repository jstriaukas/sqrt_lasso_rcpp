setwd('~/Documents/GitHub/hachdr/rcpp/')
Rcpp::sourceCpp('~/Documents/GitHub/hachdr/rcpp/sqrt_lasso.cpp')
source("rsqrtlasso.R")
require(R.matlab)
p <- 200
n <- 100

# covariates:
X <- matrix(rnorm(n*p),nrow = n, ncol = p)
# slope:
beta <- matrix(runif(p),nrow = p, ncol = 1)
# randomly assign sparse beta:
ind <- sample(1:190, replace = FALSE)
beta[ind] <- 0#beta[ind]/100
# dependent variable:
Y <- X%*%beta + matrix(rnorm(n,1),nrow = n, ncol = 1)/10

beta0 <- beta*0
writeMat("data.mat", Y=Y,X=X)
beta <- rsqrtlasso(X,Y,standardize=FALSE,Lambda=NULL)$beta.min
writeMat("beta.mat", beta=beta)




