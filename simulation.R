#This document shows the R code for our simulation studies (for a particular scenario)
#You need to download the R file functions.R in your working path.
source("functions.R")

#The following R code is for only one simulation. You can run the file repeatedly to estimate errors more accurately.
n <- 200
T <- 20
k <- 2
  
#generate Z
Z <- matrix(runif(n*k*(2^k),-1,1),n*(2^k),k)
Z <- Z[apply(Z,1,function(x) sum(x^2)) < 1,][1:n,]
Z <- Z - rep(1,n) %*% t(rep(1,n)) %*% Z / n
Z <- sqrt(n) * Z / (sum((Z%*%t(Z))^2))^0.25
  
#generate alpha
alpha <- matrix(runif(n*T,-2,0),n,T)

#generate A
P <- array(0,dim = c(n,n,T))
A <- array(0,dim = c(n,n,T))
for(t in 1:T){
  P[,,t] <- exp(outer(alpha[,t],alpha[,t],'+') + Z %*% t(Z))
  A[,,t] <- matrix(rpois(n*n,as.vector(P[,,t])),n,n)
  A[,,t][upper.tri(A[,,t])] <- t(A[,,t])[upper.tri(A[,,t])]
}

#estimation
est.c <- PGD.G(A,0.2*sqrt(n*T))
G_hat <- est.c$G
est.nc <- PGD.panel(A,ncol(Z))
Z_hat <- Newton(A,est.nc$Z,est.nc$alpha)
  
#error of G/Z for convex
err <- NULL
err <- c(err,sum((G_hat - Z %*% t(Z))^2)/n)
G.eigen <- eigen(G_hat)
Z_hatc <- G.eigen$vectors %*% rbind(diag(sqrt(G.eigen$values[1:ncol(Z)])),matrix(0,n-ncol(Z),ncol(Z)))
Oc <- (svd(t(Z_hatc) %*% Z)$v) %*% t(svd(t(Z_hatc) %*% Z)$u)
err <- c(err,sum((Z_hatc - Z %*% Oc)^2))

#error of G/Z for non-convex
err <- c(err,sum((Z_hat %*% t(Z_hat) - Z %*% t(Z))^2)/n)
O <- (svd(t(Z_hat) %*% Z)$v) %*% t(svd(t(Z_hat) %*% Z)$u)
err <- c(err,sum((Z_hat - Z %*% O)^2))
