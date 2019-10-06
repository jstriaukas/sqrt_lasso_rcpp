clear; clc;

load data
[n,p]=size(X);

lam = 1.1*sqrt(n)*norminv(1-(0.05/(2*p)), 0,1);
beta_mt = sqrtlasso(X,Y,lam,ones(200,1));

load beta
[beta, beta_mt]
norm(beta-beta_mt)


D = [[1 -1 0 0 0]; [0 1 -1 0 0]; [0 0 1 -1 0];[ 0 0 0 1 -1]];
x = randn(100,5);
y = randn(100,1);

Dtilde = D