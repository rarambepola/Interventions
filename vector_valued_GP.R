#try out both expressions for GP regression
library(MASS)
library(ggplot2)
set.seed(1)

plot.new()
par(mfrow=c(2, 2))

#polynomial regression
phi <- function(x){
  return(c(1, x, x^2, x^3, x^4, x^4))
}

N <- 50
x_max <- 3
x <- runif(N, 0, x_max)
y <- sin(-x) + rnorm(N, sd=0.1) + 1

beta_treat <- - 0.5
W <- rbinom(N, size=1, prob= 1 / (1 + exp(-1 * x)))

y[W == 1] <- y[W == 1] + beta_treat

plot(x[W==1], y[W==1], pch=19, col="red", xlim=c(min(x), max(x)),
     ylim=c(min(y), max(y)),
     main="fit separately (dashed is RFF approx)")
points(x[W==0], y[W==0], pch=19, col="blue")

#fit separate GPs to treated and untreated
k <- function(x, y, lambda, sigma){
  d <- sum((x-y)^2)
  return(sigma^2 * exp(-d / (2 * lambda^2)))
}
lambda_rbf <- 0.5
sigma_rbf <- 2

N_pred <- 200
x_pred <- sort(runif(N_pred, 0, x_max))

#treated data
N_treat <- sum(W)
x_treat <- x[W==1]
K <- matrix(NA, N_treat, N_treat)
for(i in 1:N_treat){
  for(j in 1:N_treat){
    K[i, j] <- k(x_treat[i], x_treat[j], lambda_rbf, sigma_rbf)
  }
}

K_star <- matrix(NA, N_pred, N_treat)
for(i in 1:N_pred){
  for(j in 1:N_treat){
    K_star[i, j] <- k(x_pred[i], x_treat[j], lambda_rbf, sigma_rbf)
  }
}

K_star_star <- matrix(NA, N_pred, N_pred)
for(i in 1:N_pred){
  for(j in 1:N_pred){
    K_star_star[i, j] <- k(x_pred[i], x_pred[j], lambda_rbf, sigma_rbf)
  }
}

f_pred_treat_mean <- K_star %*% solve(K + diag(N_treat)*sigma_n^2, y[W==1])

##untreated
#treated data
N_untreat <- sum(!W)
x_untreat <- x[W==0]
K <- matrix(NA, N_untreat, N_untreat)
for(i in 1:N_untreat){
  for(j in 1:N_untreat){
    K[i, j] <- k(x_untreat[i], x_untreat[j], lambda_rbf, sigma_rbf)
  }
}

K_star <- matrix(NA, N_pred, N_untreat)
for(i in 1:N_pred){
  for(j in 1:N_untreat){
    K_star[i, j] <- k(x_pred[i], x_untreat[j], lambda_rbf, sigma_rbf)
  }
}

f_pred_untreat_mean <- K_star %*% solve(K + diag(N_untreat)*sigma_n^2, y[W==0])
lines(x_pred, f_pred_treat_mean, col="red")
lines(x_pred, f_pred_untreat_mean, col="blue")


##do separate fits with RFF
N_features <- 300
w <- rnorm(N_features, mean=0, sd= 1 / lambda_rbf)

phi <- function(x, w, sigma){
  return(sigma * c(sin(x*w), cos(x*w)) / sqrt(length(w)))
}

#treated
Phi <- do.call("cbind", lapply(x[W==1], phi, w, sigma_rbf))

A <- (Phi %*% t(Phi)) / (sigma_n^2) + diag(dim(Phi)[1]) 
Phi_star <- do.call("cbind", lapply(x_pred, phi, w, sigma_rbf))
f_pred_treat_mean <- t(Phi_star) %*% solve(A, Phi %*% y[W==1]) / sigma_n^2

#untreated
#treated
Phi <- do.call("cbind", lapply(x[W==0], phi, w, sigma_rbf))

A <- (Phi %*% t(Phi)) / (sigma_n^2) + diag(dim(Phi)[1]) 
Phi_star <- do.call("cbind", lapply(x_pred, phi, w, sigma_rbf))
f_pred_untreat_mean <- t(Phi_star) %*% solve(A, Phi %*% y[W==0]) / sigma_n^2


lines(x_pred, f_pred_treat_mean, col="red", lty=2)
lines(x_pred, f_pred_untreat_mean, col="blue", lty=2)





#fit both with the same GP

beta00 <- 1
beta01 <- 1
rho <- 0.9 * beta00 * beta01

k_both <- function(x, y, w_x, w_y, lambda, sigma){
  d <- sum((x-y)^2)
  k_baseline <- k(x, y, lambda, sigma)
  
  if(w_x==w_y){
    if(w_x == 1) return(k_baseline * beta01^2) else return(k_baseline * beta00^2)
  }else{
    return(k_baseline * rho)
  }
}

K <- matrix(NA, N, N)
for(i in 1:N){
  for(j in 1:N){
    K[i, j] <- k_both(x[i], x[j], W[i], W[j], lambda_rbf, sigma_rbf)
  }
}

##treated
K_star <- matrix(NA, N_pred, N)
for(i in 1:N_pred){
  for(j in 1:N){
    K_star[i, j] <- k_both(x_pred[i], x[j], 1, W[j], lambda_rbf, sigma_rbf)
  }
}

K_star_star <- matrix(NA, N_pred, N_pred)
for(i in 1:N_pred){
  for(j in 1:N_pred){
    K_star_star[i, j] <- k_both(x_pred[i], x_pred[j], 1, 1, lambda_rbf, sigma_rbf)
  }
}

f_pred_treat_mean <- K_star %*% solve(K + diag(N)*sigma_n^2, y)
f_pred_treat_sd <- K_star_star - K_star %*% solve(K + diag(N)*sigma_n^2, t(K_star))

##untreated
K_star <- matrix(NA, N_pred, N)
for(i in 1:N_pred){
  for(j in 1:N){
    K_star[i, j] <- k_both(x_pred[i], x[j], 0, W[j], lambda_rbf, sigma_rbf)
  }
}

K_star_star <- matrix(NA, N_pred, N_pred)
for(i in 1:N_pred){
  for(j in 1:N_pred){
    K_star_star[i, j] <- k_both(x_pred[i], x_pred[j], 1, 1, lambda_rbf, sigma_rbf)
  }
}

f_pred_untreat_mean <- K_star %*% solve(K + diag(N)*sigma_n^2, y)
f_pred_untreat_sd <- K_star_star - K_star %*% solve(K + diag(N)*sigma_n^2, t(K_star))

plot(x[W==1], y[W==1], pch=19, col="red", xlim=c(min(x), max(x)),
     ylim=c(min(y), max(y)),
     main="fit jointly (dashed is RFF)")
points(x[W==0], y[W==0], pch=19, col="blue")

lines(x_pred, f_pred_treat_mean, col="red")
lines(x_pred, f_pred_untreat_mean, col="blue")




##do joint fit with RFF
N_features <- 300
w <- rnorm(N_features, mean=0, sd= 1 / lambda_rbf)

phi <- function(x, w, sigma, W){
  baseline_features <- sigma * c(sin(x*w), cos(x*w)) / sqrt(length(w))
  m <- length(baseline_features)
  if(W == 0) return(c(baseline_features * rho, 
                      baseline_features * sqrt(beta01^2 - rho),
                      rep(0, m)))
  
  return(c(baseline_features * rho, 
           rep(0, m),
           baseline_features * sqrt(beta00^2 - rho)))
}

#treated
Phi <- c()
for(i in 1:N){
  Phi <- cbind(Phi,
               phi(x[i], w, sigma_rbf, W[i]))
}

Phi_star <- c()
for(i in 1:N_pred){
  Phi_star <- cbind(Phi_star,
                    phi(x_pred[i], w, sigma_rbf, 1))
}
A <- (Phi %*% t(Phi)) / (sigma_n^2) + diag(dim(Phi)[1]) 
f_pred_treat_mean <- t(Phi_star) %*% solve(A, Phi %*% y) / sigma_n^2
f_pred_treat_sd_rff <- t(Phi_star) %*% solve(A, Phi_star)

#not treated
Phi_star <- c()
for(i in 1:N_pred){
  Phi_star <- cbind(Phi_star,
                    phi(x_pred[i], w, sigma_rbf, 0))
}
A <- (Phi %*% t(Phi)) / (sigma_n^2) + diag(dim(Phi)[1]) 
f_pred_untreat_mean <- t(Phi_star) %*% solve(A, Phi %*% y) / sigma_n^2
f_pred_untreat_sd_rff <- t(Phi_star) %*% solve(A, Phi_star)

lines(x_pred, f_pred_treat_mean, col="red", lty=2)
lines(x_pred, f_pred_untreat_mean, col="blue", lty=2)

plot(f_pred_treat_sd, f_pred_treat_sd_rff, cex=0.5)
abline(0, 1)
plot(f_pred_untreat_sd, f_pred_untreat_sd_rff, cex=0.5)
abline(0, 1)

par(mfrow=c(1, 1))

