setwd(paste0(Sys.getenv("HOME"), "/Interventions/"))
library(ggplot2)

set.seed(2)

##mock data for testing
N_obs <- 100
x_max <- 10
x <- sort(runif(N_obs, 0, x_max))
y <- rnorm(N_obs, mean=0, sd=0.2)
# plot(x, y)
w_prob <- rep(0.5, N_obs)
w_prob[x > 8] <- 0.01
w <- rbinom(N_obs, 1, w_prob)
# w <- rbinom(N_obs, 1, 0.5)

treat_effect <- 0.75
y_obs <- rnorm(N_obs, mean= w*treat_effect, sd=0.2)

# 
# ##do in R first to get code right
# #normal GP regression
# make_rff <- function(w, x, lambda=1, sigma=1){
#   if(!is.matrix(x)) x <- matrix(x, ncol=1)
#   w <- w / lambda
#   wx <- x %*% w
#   n_rff <- dim(w)[2]
#   return(sigma * cbind(sin(wx), cos(wx))/sqrt(n_rff))
# }
# 
# sigma_n <- 1
# lambda <- 1
# sigma <- 1
# n_rff <- 10
# ws <- matrix(rnorm(n_rff), nrow=1)
# rff <- make_rff(ws, x, lambda=lambda, sigma=sigma)
# A <- (t(rff) %*% rff + diag(n_rff * 2)) / (sigma_n^2)
# y_fit <- rff %*% solve(A, t(rff) %*% y) / (sigma_n^2)

# N_pred <- 1000
# x_pred <- sort(seq(0, x_max, length.out=N_pred))
# rff_pred <- make_rff(ws, x_pred, lambda, sigma)
# y_pred <- rff_pred %*% solve(A, t(rff) %*% y) / (sigma_n^2)
# 
# plot(x, y)
# lines(x_pred, y_pred, col="red")

#regression with treatment effect
make_rff2 <- function(w, b, x, treat_vec, lambda=1, sigma=1, rho=0){
  if(!is.matrix(x)) x <- matrix(x, ncol=1)
  w <- w / lambda
  n_rff <- dim(w)[2]
  wx1 <- (x %*% w) 
  wx1 <- wx1 + runif(length(wx1), 0, 2*pi)
  wx2 <- wx1
  wx3 <- wx1
  
  wx1 <- wx1 * sqrt(rho)
  coeff <- sqrt(sigma^2 - rho)
  wx2 <- wx2 * coeff
  wx2[treat_vec == 1, ] <- 0
  wx3 <- wx3 * coeff
  wx3[treat_vec == 0,] <- 0
  wx <- cbind(wx1, wx2, wx3)
  return(sigma*cbind(sin(wx), cos(wx)) / sqrt(3 * n_rff))
}

sigma_n <- 1
lambda <- 0.001
sigma <- 1
rho_scale <- 0
rho <- sigma^2 * rho_scale
n_rff <- 10
ws <- matrix(rnorm(n_rff), nrow=1)
bs <- runif(n_rff, 0, pi) 
rff <- make_rff2(ws, bs, x, w, lambda=lambda, sigma=sigma, rho=rho)
A <- (t(rff) %*% rff + diag(n_rff * 6)) / (sigma_n^2)
y_fit <- rff %*% solve(A, t(rff) %*% y_obs) / (sigma_n^2)


N_pred <- 2
x_pred <- rep(sort(seq(0, x_max, length.out=N_pred)), 2)
treat_pred <- rep(0:1, each=N_pred)
rff_pred <- make_rff2(ws, bs, x_pred, treat_pred, lambda, sigma, rho)
y_pred <- rff_pred %*% solve(A, t(rff) %*% y_obs) / (sigma_n^2)

plot(x, y_obs, cex=0.001)
points(x[w==0], y_obs[w==0], pch=19, col="red", cex=0.5)
points(x[w==1], y_obs[w==1], pch=19, col="blue", cex=0.5)
lines(x_pred[treat_pred == 0], y_pred[treat_pred == 0], col="red")
lines(x_pred[treat_pred == 1], y_pred[treat_pred == 1], col="blue")
