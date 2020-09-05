#try out both expressions for GP regression
library(MASS)
library(ggplot2)
# set.seed(1)
#polynomial regression
phi <- function(x){
  return(c(1, x, x^2, x^3, x^4, x^4))
}

N <- 100
x <- runif(N, 0, 2)
y <- sin(x) + rnorm(N, sd=0.5)

plot(x, y)

sigma_p <- 0.4
sigma_n <- 0.6
Phi <- do.call("cbind", lapply(x, phi))

N_pred <- 100
x_pred <- sort(runif(N_pred, 0, 2))

A <- (Phi %*% t(Phi)) / (sigma_n^2) + diag(dim(Phi)[1]) / sigma_p^2
Phi_star <- do.call("cbind", lapply(x_pred, phi))
f_pred_mean <- t(Phi_star) %*% solve(A, Phi %*% y) / sigma_n^2

lines(x_pred, f_pred_mean)

f_pred_sd <- t(Phi_star) %*% solve(A, Phi_star)
N_sample <- 100
y_sample <- mvrnorm(N_sample, f_pred_mean, f_pred_sd)

pred_df <- data.frame(x=x_pred,
                      y=f_pred_mean,
                      upper=apply(y_sample, 2, quantile, 0.975),
                      lower=apply(y_sample, 2, quantile, 0.025)
)

obs_df <- data.frame(x=x, 
                     y=y)

p <- ggplot() + geom_line(data=pred_df, aes(x=x, y=y)) + 
  geom_ribbon(data=pred_df, aes(x=x, ymin=lower, ymax=upper), alpha=0.3) + 
  geom_point(data=obs_df, aes(x=x, y=y)) + 
  theme_minimal()
# print(p)

#now other way round
k <- function(x, y){
  return(sigma_p^2 * sum(phi(x) * phi(y)))
}

K <- matrix(NA, N, N)
for(i in 1:N){
  for(j in 1:N){
    K[i, j] <- k(x[i], x[j])
  }
}

K_star <- matrix(NA, N_pred, N)
for(i in 1:N_pred){
  for(j in 1:N){
    K_star[i, j] <- k(x_pred[i], x[j])
  }
}

K_star_star <- matrix(NA, N_pred, N_pred)
for(i in 1:N_pred){
  for(j in 1:N_pred){
    K_star_star[i, j] <- k(x_pred[i], x_pred[j])
  }
}

f_pred_mean2 <- K_star %*% solve(K + diag(N)*sigma_n^2, y)
# plot(f_pred_mean, f_pred_mean2)
plot(x, y, main="polynomial")
lines(x_pred, f_pred_mean2)
lines(x_pred, f_pred_mean, col="red")

f_pred_sd2 <- K_star_star - K_star %*% solve(K + diag(N)*sigma_n^2, t(K_star))

# plot(f_pred_sd, f_pred_sd2)
# abline(0, 1)

###step 2, do same thing with rbf kernel
#define features
N_features <- 100
lambda_rbf <- 0.5
sigma_rbf <- 2

#1D so w can be a vector
# w <- rnorm(N_features, mean=0, sd= 1 / lambda_rbf)
w <- rnorm(N_features) / lambda_rbf

phi <- function(x, w, sigma){
  return(sigma * c(sin(x*w), cos(x*w)) / sqrt(length(w)))
}

Phi <- do.call("cbind", lapply(x, phi, w, sigma_rbf))

A <- (Phi %*% t(Phi)) / (sigma_n^2) + diag(dim(Phi)[1]) / sigma_p^2
Phi_star <- do.call("cbind", lapply(x_pred, phi, w, sigma_rbf))
f_pred_mean <- t(Phi_star) %*% solve(A, Phi %*% y) / sigma_n^2
f_pred_sd <- t(Phi_star) %*% solve(A, Phi_star)

# plot(x, y, main="rbf")
# lines(x_pred, f_pred_mean)


#with exact kernel
k <- function(x, y, lambda, sigma){
  d <- sum((x-y)^2)
  return(sigma_p^2 * sigma^2 * exp(-d / (2 * lambda^2)))
}

K <- matrix(NA, N, N)
for(i in 1:N){
  for(j in 1:N){
    K[i, j] <- k(x[i], x[j], lambda_rbf, sigma_rbf)
  }
}

K_star <- matrix(NA, N_pred, N)
for(i in 1:N_pred){
  for(j in 1:N){
    K_star[i, j] <- k(x_pred[i], x[j], lambda_rbf, sigma_rbf)
  }
}

K_star_star <- matrix(NA, N_pred, N_pred)
for(i in 1:N_pred){
  for(j in 1:N_pred){
    K_star_star[i, j] <- k(x_pred[i], x_pred[j], lambda_rbf, sigma_rbf)
  }
}

f_pred_mean2 <- K_star %*% solve(K + diag(N)*sigma_n^2, y)
# plot(f_pred_mean, f_pred_mean2)
plot(x, y, main="rbf (approximation in red)")
lines(x_pred, f_pred_mean2)
lines(x_pred, f_pred_mean, col="red")




