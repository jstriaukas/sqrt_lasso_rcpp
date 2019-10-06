#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec sqrtlasso(const arma::mat& X, const arma::vec& Y, double lambda1, double maxIter, double OptTolNorm, double OptTolObj, double standardize, arma::vec beta0, double init) {
  // double zeroThreshold
  int n = X.n_rows, p = X.n_cols;
  double Iter;
  double MaxErrorNorm = 1.0e-10;
  arma::mat X1 = X;
  arma::vec Y1 = Y;
  arma::vec tmp, error;
  
  if (standardize == 1){
    arma::mat M1 = mean(X,0);
    arma::mat S1 = stddev(X,0,0);
    X1 = (X-repmat(M1, n, 1))/repmat(S1, n, 1);
    double M2 = mean(Y);
    Y1 = Y - M2;
  } else {
    X1 = X;
    Y1 = Y;
  }

  // =============================================
  // Initial Point: start from the Ridge estimator
  // =============================================
  arma::mat XX = X1.t()*X1;
  arma::mat Xy = X1.t()*Y1;
  arma::mat M  = eye(p,p)*lambda1;
  if (init){
   beta0 = pinv((XX + M) * eye(p, p)) * Xy; 
  }
  arma::vec beta = beta0;
  XX = XX/n;    // Gram matrix             
  Xy = Xy/n;  
  error = Y1 - X1*beta;     // residuals
  double qhat = sum( square(error) )/n;   // average of squared residuals
  
  Iter = 0;
  
  while (Iter < maxIter){  // while loop until convergence beta^(i)
    arma::vec beta_old = beta;
    Iter = Iter + 1;
    for(int j = 0; j < p; j++){ // loop through beta^(i)_j
      // compute the Shoot and Update the variable
       double S0 = sum(XX.row(j) * beta) - XX(j,j) * beta[j] - Xy[j];
      // error = y - X*beta +  X(:,j)*beta(j);
       if (fabs(beta[j]) > 0){
         error = error + X.col(j)*beta[j];
         qhat = sum( square(error) )/n;
       }
      
      // Note that by C-S
      // S0^2 <= Qhat * XX(j,j), so that  Qhat >= S0^2/XX(j,j)  :)
      if ( pow(n,2) < pow(lambda1,2) / XX(j,j) ) {
        beta[j] = 0;
      }

      double tmp = (lambda1/n)  * sqrt(qhat);
      double qqhat = qhat - (pow(S0,2)/XX(j,j));
      if (qqhat<0){
        qqhat = 0; //max(qqhat,0)
      }
      if (S0 > tmp) {// Optimal beta(j) < 0
        beta[j] = (  ( lambda1 / sqrt( pow(n,2) - pow(lambda1,2) / XX(j,j)  ) ) * sqrt(qqhat)  - S0 )   /   XX(j,j);
        error = error - X.col(j)*beta[j];
      }
      if (S0 < - tmp){// Optimal beta(j) > 0
        beta[j] = ( - ( lambda1 / sqrt( pow(n,2) - pow(lambda1,2) / XX(j,j)  ) ) * sqrt(qqhat)  - S0 )   /   XX(j,j);
        error = error - X.col(j)*beta[j];
      }

      if (fabs(S0) <= tmp){ // Optimal beta(j) = 0
        beta[j] = 0;
      }
    } // end loop beta^(i)_j

    // Update primal and dual value
    double fobj = sqrt( sum(pow(X*beta-Y1,2))/n )  +  sum(lambda1*abs(beta)/n);
    
    double ErrorNorm = norm(error);
    if (  ErrorNorm > MaxErrorNorm ) { 
      arma::vec aaa  = (sqrt(n)*error/ErrorNorm);
      arma::vec bbb = abs(  lambda1/n - abs(X.t()*aaa/n) );
      double dual = sum(aaa%Y1)/n - sum(bbb%abs(beta));
      // check convergence
      if (fobj - dual  < OptTolObj) {
        if (norm(beta - beta_old,1) < OptTolNorm) {
          Iter = maxIter + 1;
        }
      }
    } else {
      double dual = sum(lambda1*abs(beta)/n);
      // check convergence
      if (fobj - dual  < OptTolObj) {
        if (norm(beta - beta_old,1) < OptTolNorm) {
          Iter = maxIter + 1;
        }
      }
    }
  } // end while loop
  return beta;
}

// [[Rcpp::export]]
arma::mat sqrtlassogrid(const arma::mat& X, const arma::vec& Y, const arma::vec& Lambda, double maxIter, double OptTolNorm, double OptTolObj, double standardize) {

  // preliminaries:
  int  p = X.n_cols, lam1 = Lambda.size();
  arma::mat BETA(p,lam1);
  arma::vec beta0(p,1);
  // compute first beta to initialize at it, then use warm starts for remaining beta's:
  arma::vec beta = sqrtlasso(X,Y,Lambda[0],maxIter,OptTolNorm,OptTolObj,standardize,beta0,1);
  for (int l1 = 0; l1 < lam1; l1++){ // loop through Lambda1
    // select lambda:
    double lambda1 = Lambda[l1];
    // compute solution for lambda1, initialize at previous estimate of beta:
    beta = sqrtlasso(X,Y,lambda1,maxIter,OptTolNorm,OptTolObj,standardize,beta,0);
    BETA.col(l1) = beta;//BETA  store beta 
  }  // end loop Lambda1   
  return BETA;
}
