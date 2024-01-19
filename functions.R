#R function for implementing one step newton update
#Input: longitudinal network data A (n*n*T dimensional array), initial estimates for Z (n*k dimensional matrix) and alpha (n*T dimensional matrix)
#Output: The proposed semiparametric one-step estimator of Z (n*k dimensional matrix)
Newton <- function(A,Z,alpha){
  if(length(dim(A)) == 2)
    A <- array(as.vector(A),c(nrow(A),nrow(A),1))
  T <- dim(A)[3]
  n <- nrow(Z)
  k <- ncol(Z)
  
  U <- svd(Z)$u
  U <- cbind(U,rep(1,n)/sqrt(n))
  U <- cbind(U,svd(U,nu=n)$u[,(k+2):n])
  V <- t(svd(Z)$v)
  D <- svd(Z)$d
  
  UZ <- NULL
  for(i in 1:k)
    UZ <- cbind(UZ,as.vector(V[i,] %*% t(U[,i])))
  for(j in (k+2):n)
    for(i in 1:k)
      UZ <- cbind(UZ,as.vector(V[i,] %*% t(U[,j])))
  for(i in 1:(k-1))
    for(j in (i+1):k)
      UZ <- cbind(UZ,(D[j]*as.vector(V[j,] %*% t(U[,i]))+D[i]*as.vector(V[i,] %*% t(U[,j])))/sqrt(D[i]^2+D[j]^2))
  
  Ieff <- matrix(0,n*k,n*k)
  Seff <- rep(0,n*k)
  for(t in 1:T){
    theta <- alpha[,t] %*% t(rep(1,n)) + rep(1,n) %*% t(alpha[,t]) + Z %*% t(Z)
    mu <- exp(theta)
    M <- A[,,t] - mu
      
    Lzz <- matrix(0,n*k,n*k)
    for(i in 1:n)
      for(j in 1:n){
        if(i == j){
          Lzz[((i-1)*k+1):(i*k),((j-1)*k+1):(j*k)] <- 3*mu[i,i] * Z[i,] %*% t(Z[i,])
          for(jj in 1:n)
            Lzz[((i-1)*k+1):(i*k),((j-1)*k+1):(j*k)] <- Lzz[((i-1)*k+1):(i*k),((j-1)*k+1):(j*k)] + mu[i,jj]*Z[jj,] %*% t(Z[jj,])
        }
        else
          Lzz[((i-1)*k+1):(i*k),((j-1)*k+1):(j*k)] <- mu[i,j] * Z[j,] %*% t(Z[i,])
      }
    
    Lza <- matrix(0,n*k,n)
    for(i in 1:n)
      for(j in 1:n){
        if(i == j){
          Lza[((i-1)*k+1):(i*k),j] <- 3*mu[i,i] * Z[i,]
          for(jj in 1:n)
            Lza[((i-1)*k+1):(i*k),j] <- Lza[((i-1)*k+1):(i*k),j] + mu[i,jj] * Z[jj,] 
        }
        else
          Lza[((i-1)*k+1):(i*k),j] <- mu[i,j] * Z[j,]
      }
    
    Laa <- matrix(0,n,n)
    for(i in 1:n)
      for(j in 1:n){
        if(i == j){
          Laa[i,j] <- 3*mu[i,i]
          for(jj in 1:n)
            Laa[i,j] <- Laa[i,j] + mu[i,jj]
        }
        else
          Laa[i,j] <- mu[i,j]
      }
    
    Lz <- rep(0,n*k)
    for(i in 1:n){
      Lz[((i-1)*k+1):(i*k)] <- Z[i,] * M[i,i]
      for(j in 1:n)
        Lz[((i-1)*k+1):(i*k)] <- Lz[((i-1)*k+1):(i*k)] + Z[j,] * M[i,j]
    }
    
    La <- rep(0,n)
    for(i in 1:n){
      La[i] <- M[i,i]
      for(j in 1:n)
        La[i] <- La[i] + M[i,j]
    }
    
    Ieff <- Ieff + Lzz - Lza %*% solve(Laa) %*% t(Lza)
    Seff <- Seff + Lz - Lza %*% solve(Laa) %*% La
  }
  vz <- as.vector(t(Z))
  vz <- vz + UZ %*% solve(t(UZ) %*% Ieff %*% UZ) %*% t(UZ) %*% Seff
  return(t(matrix(vz,k,n)))
}

#R function for implementing the initial algorithm
#Input: longitudinal network data A (n*n*T dimensional array), latent space dimension k (numeric)
#Output: initial estimates for Z (n*k dimensional matrix) and alpha (n*T dimensional matrix)
PGD.panel <- function(A,k){
  if(length(dim(A)) == 2)
    A <- array(as.vector(A),c(nrow(A),nrow(A),1))
  
  n <- dim(A)[1]
  T <- dim(A)[3]
  N <- apply(A,c(1,2),sum)
  
  #initial of initial
  p_hat <- sum(N)/(n^2)
  tau <- sqrt(n * p_hat)
  svdN <- svd(N)
  P_hat <- svdN$u %*% diag(svdN$d * I(svdN$d > tau)) %*% t(svdN$v)
  P_hat <- 0.01*T * I(P_hat < 0.01*T) + P_hat * I(0.01*T <= P_hat & P_hat <= T) + T * I(P_hat > T)
  Theta_hat <- log(P_hat/T)
  alpha0 <- solve(n*diag(rep(1,n))+rep(1,n)%*%t(rep(1,n))) %*% Theta_hat %*% rep(1,n)
  J <- diag(rep(1,n)) - rep(1,n) %*% t(rep(1,n)) / n
  R <- Theta_hat - (alpha0%*%t(rep(1,n)) + rep(1,n)%*%t(alpha0))
  R <- eigen(R)$vectors %*% diag(eigen(R)$values * I(eigen(R)$values > 0)) %*% t(eigen(R)$vectors)
  EV <- eigen(R)$values[1:10]
  if(k == 1)
    Z0 <- eigen(R)$vectors %*% rbind(diag(as.matrix(sqrt(eigen(R)$values[1:k]))),matrix(0,n-k,k))
  else
    Z0 <- eigen(R)$vectors %*% rbind(diag(sqrt(eigen(R)$values[1:k])),matrix(0,n-k,k))
  Z.init <- Z0
  alpha0 <- matrix(rep(alpha0,T),n,T)
  
  #gradient descent (BB step size)
  eta <- 0.1
  eta.Z <- eta / (T*svd(Z0)$d[1]^2)
  eta.a <- eta / (2*n)
  
  inner <- exp(Z0 %*% t(Z0))
  cum.in <- exp(alpha0) %*% t(exp(alpha0))
  M <- inner * cum.in
  Q <- apply(alpha0,2,function(x) apply(inner * exp(outer(x,x,'+')),2,sum))
  grad0.Z <- (N - M) %*% Z0
  grad0.a <- apply(A,c(1,3),sum) - Q
  
  Z1 <- Z0 + eta.Z * grad0.Z
  Z1 <- J %*% Z1
  alpha1 <- alpha0 + eta.a * grad0.a
  
  for(i in 1:999){
    inner <- exp(Z1 %*% t(Z1))
    cum.in <- exp(alpha1) %*% t(exp(alpha1))
    M <- inner * cum.in
    Q <- apply(alpha1,2,function(x) apply(inner * exp(outer(x,x,'+')),2,sum))
    grad1.Z <- (N - M) %*% Z1
    grad1.a <- apply(A,c(1,3),sum) - Q
    
    eta.Z <- -sum(diag(t(Z1-Z0) %*% (grad1.Z-grad0.Z)))/sum(diag(t(grad1.Z-grad0.Z) %*% (grad1.Z-grad0.Z)))
    eta.a <- -sum(diag(t(alpha1-alpha0) %*% (grad1.a-grad0.a)))/sum(diag(t(grad1.a-grad0.a) %*% (grad1.a-grad0.a)))
    Z0 <- Z1
    alpha0 <- alpha1
    grad0.Z <- grad1.Z
    grad0.a <- grad1.a
    Z1 <- Z0 + eta.Z * grad0.Z
    Z1 <- J %*% Z1
    alpha1 <- alpha0 + eta.a * grad0.a
  }
  return(list(Z = Z1,alpha = alpha1))
}

#R function for implementing the initial algorithm (version without diagonal data)
#Input: longitudinal network data A (n*n*T dimensional array), latent space dimension k (numeric)
#Output: initial estimates for Z (n*k dimensional matrix) and alpha (n*T dimensional matrix)
PGD.panel2 <- function(A,k){
  if(length(dim(A)) == 2)
    A <- array(as.vector(A),c(nrow(A),nrow(A),1))
  
  n <- dim(A)[1]
  T <- dim(A)[3]
  N <- apply(A,c(1,2),sum)
  
  #initial value
  p_hat <- sum(N)/(n^2)
  tau <- sqrt(n * p_hat)
  svdN <- svd(N)
  P_hat <- svdN$u %*% diag(svdN$d * I(svdN$d > tau)) %*% t(svdN$v)
  P_hat <- 0.01*T * I(P_hat < 0.01*T) + P_hat * I(0.01*T <= P_hat & P_hat <= T) + T * I(P_hat > T)
  Theta_hat <- log(P_hat/T)
  alpha0 <- solve(n*diag(rep(1,n))+rep(1,n)%*%t(rep(1,n))) %*% Theta_hat %*% rep(1,n)
  J <- diag(rep(1,n)) - rep(1,n) %*% t(rep(1,n)) / n
  R <- Theta_hat - (alpha0%*%t(rep(1,n)) + rep(1,n)%*%t(alpha0))
  R <- eigen(R)$vectors %*% diag(eigen(R)$values * I(eigen(R)$values > 0)) %*% t(eigen(R)$vectors)
  EV <- eigen(R)$values[1:10]
  if(k == 1)
    Z0 <- eigen(R)$vectors %*% rbind(diag(as.matrix(sqrt(eigen(R)$values[1:k]))),matrix(0,n-k,k))
  else
    Z0 <- eigen(R)$vectors %*% rbind(diag(sqrt(eigen(R)$values[1:k])),matrix(0,n-k,k))
  Z.init <- Z0
  alpha0 <- matrix(rep(alpha0,T),n,T)
  
  #gradient descent (BB)
  eta <- 0.1
  eta.Z <- eta / (T*svd(Z0)$d[1]^2)
  eta.a <- eta / (2*n)
  
  inner <- exp(Z0 %*% t(Z0)) - diag(diag(exp(Z0 %*% t(Z0))))
  cum.in <- exp(alpha0) %*% t(exp(alpha0))
  M <- inner * cum.in
  Q <- apply(alpha0,2,function(x) apply(inner * exp(outer(x,x,'+')),2,sum))
  grad0.Z <- (N - M) %*% Z0
  grad0.a <- apply(A,c(1,3),sum) - Q
  
  Z1 <- Z0 + eta.Z * grad0.Z
  Z1 <- J %*% Z1
  alpha1 <- alpha0 + eta.a * grad0.a
  
  for(i in 1:999){
    inner <- exp(Z1 %*% t(Z1)) - diag(diag(exp(Z1 %*% t(Z1))))
    cum.in <- exp(alpha1) %*% t(exp(alpha1))
    M <- inner * cum.in
    Q <- apply(alpha1,2,function(x) apply(inner * exp(outer(x,x,'+')),2,sum))
    grad1.Z <- (N - M) %*% Z1
    grad1.a <- apply(A,c(1,3),sum)- Q
    
    eta.Z <- -sum(diag(t(Z1-Z0) %*% (grad1.Z-grad0.Z)))/sum(diag(t(grad1.Z-grad0.Z) %*% (grad1.Z-grad0.Z)))
    eta.a <- -sum(diag(t(alpha1-alpha0) %*% (grad1.a-grad0.a)))/sum(diag(t(grad1.a-grad0.a) %*% (grad1.a-grad0.a)))
    Z0 <- Z1
    alpha0 <- alpha1
    grad0.Z <- grad1.Z
    grad0.a <- grad1.a
    Z1 <- Z0 + eta.Z * grad0.Z
    Z1 <- J %*% Z1
    alpha1 <- alpha0 + eta.a * grad0.a
  }
  return(list(Z = Z1,alpha = alpha1))
}
               
#R function for optimizing the penalized log-likelihood
#Input: longitudinal network data A (n*n*T dimensional array), tuning parameter lambda (numeric)
#Output: penalized mle for G (n*n dimensional matrix) and alpha (n*T dimensional matrix)
PGD.G <- function(A,lambda){
  if(length(dim(A)) == 2)
    A <- array(as.vector(A),c(nrow(A),nrow(A),1))
  
  n <- dim(A)[1]
  T <- dim(A)[3]
  N <- apply(A,c(1,2),sum)
  
  #initial value
  p_hat <- sum(N)/(n^2)
  tau <- sqrt(n * p_hat)
  svdN <- svd(N)
  P_hat <- svdN$u %*% diag(svdN$d * I(svdN$d > tau)) %*% t(svdN$v)
  P_hat <- 0.01*T * I(P_hat < 0.01*T) + P_hat * I(0.01*T <= P_hat & P_hat <= T) + T * I(P_hat > T)
  Theta_hat <- log(P_hat/T)
  alpha0 <- solve(n*diag(rep(1,n))+rep(1,n)%*%t(rep(1,n))) %*% Theta_hat %*% rep(1,n)
  J <- diag(rep(1,n)) - rep(1,n) %*% t(rep(1,n)) / n
  G0 <- Theta_hat - (alpha0%*%t(rep(1,n)) + rep(1,n)%*%t(alpha0))
  alpha0 <- matrix(rep(alpha0,T),n,T)
  
  #proximal gradient descent (BB)
  eta <- 0.1
  eta.G <- eta / T
  eta.a <- eta / (2*n)
  mu <- 1
  
  inner <- exp(G0)
  cum.in <- exp(alpha0) %*% t(exp(alpha0))
  M <- inner * cum.in
  Q <- apply(alpha0,2,function(x) apply(inner * exp(outer(x,x,'+')),2,sum))
  grad0.G <- (N - M)
  grad0.a <- apply(A,c(1,3),sum) - Q
  
  alpha1 <- alpha0 + eta.a * grad0.a
  G1.t <- G0 + eta.G * grad0.G
  #ADMM
  Y <- G1.t
  W <- G1.t
  Lambda <- diag(rep(1,n))
  for(j in 1:10){
    Y <- J %*% (W + Lambda/mu) %*% J
    W.t <- (G1.t/eta.G + mu*Y - Lambda)/(mu + 1/eta.G)
    W.eigen <- eigen(W.t)
    SW.values <- (W.eigen$values - eta.G*lambda) * I((W.eigen$values - eta.G*lambda)>0)
    W <- W.eigen$vectors %*% diag(SW.values) %*% t(W.eigen$vectors)
    Lambda <- Lambda - mu*(Y - W)
  }
  G1 <- Y
  
  for(i in 1:999){
    inner <- exp(G1)
    cum.in <- exp(alpha1) %*% t(exp(alpha1))
    M <- inner * cum.in
    Q <- apply(alpha1,2,function(x) apply(inner * exp(outer(x,x,'+')),2,sum))
    grad1.G <- (N - M)
    grad1.a <- apply(A,c(1,3),sum) - Q
    
    eta.G <- -sum(diag(t(G1-G0) %*% (grad1.G-grad0.G)))/sum(diag(t(grad1.G-grad0.G) %*% (grad1.G-grad0.G)))
    eta.a <- -sum(diag(t(alpha1-alpha0) %*% (grad1.a-grad0.a)))/sum(diag(t(grad1.a-grad0.a) %*% (grad1.a-grad0.a)))
    G0 <- G1
    alpha0 <- alpha1
    grad0.G <- grad1.G
    grad0.a <- grad1.a
    alpha1 <- alpha0 + eta.a * grad0.a
    G1.t <- G0 + eta.G * grad0.G
    Y <- G1.t
    W <- G1.t
    Lambda <- diag(rep(1,n))
    for(j in 1:10){
      Y <- J %*% (W + Lambda/mu) %*% J
      W.t <- (G1.t/eta.G + mu*Y - Lambda)/(mu + 1/eta.G)
      W.eigen <- eigen(W.t)
      SW.values <- (W.eigen$values - eta.G*lambda) * I((W.eigen$values - eta.G*lambda)>0)
      W <- W.eigen$vectors %*% diag(SW.values) %*% t(W.eigen$vectors)
      Lambda <- Lambda - mu*(Y - W)
    }
    G1 <- Y
  }
  return(list(G = G1,alpha = alpha1))
}
