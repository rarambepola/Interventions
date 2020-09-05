setwd(paste0(Sys.getenv("HOME"), "/Interventions/"))
library(ggplot2)

set.seed(1)
N_obs <- 100
x_max <- 10
x <- sort(runif(N_obs, 0, x_max))
y <- rnorm(N_obs, mean=3, sd=0.2)
# plot(x, y)
w_prob <- rep(0.5, N_obs)
w_prob[x > 3 & x < 8] <- 0.01
w <- rbinom(N_obs, 1, w_prob)
# w <- rbinom(N_obs, 1, 0.5)

treat_effect <- 0.75
y_obs <- rnorm(N_obs, mean=3 + w*treat_effect, sd=0.2)

obs_df <- data.frame(x=x, y=y_obs, w=as.factor(w))

p <- ggplot(obs_df) + geom_point(aes(x=x, y=y, col=w))
# print(p)


#define a GP
lambda <- 20
sigma <- 2
rho <- sigma^2 * 0.9

k <- function(x1, x2, w1, w2){
  d <- sum((x1 - x2)^2)
  corr <- rho + (sigma^2 - rho) * (w1 == w2)
  return(corr * exp(- d / lambda^2))
}

K <- matrix(NA, N_obs, N_obs)
for(i in 1:N_obs){
  for(j in 1:N_obs){
    K[i, j] <- k(x[i], x[j], w[i], w[j])
  }
}

N_pred <- 500
x_pred <- seq(0, x_max, length.out = N_pred)
x_pred <- rep(x_pred, 2)
w_pred <- rep(0:1, each=N_pred)
N_pred <- N_pred * 2

K_pred <- matrix(NA, N_pred, N_obs)
for(i in 1:N_pred){
  for(j in 1:N_obs){
    K_pred[i, j] <- k(x_pred[i], x[j], w_pred[i], w[j])
  }
}

lambda_GP <- 0.001
y_pred <- K_pred %*% solve(K + lambda_GP * diag(N_obs), y_obs)

pred_df <- data.frame(x=x_pred, y=y_pred, w=as.factor(w_pred))

p <- p + geom_line(data=pred_df, aes(x=x, y=y, col=w)) + theme_minimal()

# print(p)


#######
w_prob <- rep(0.5, N_obs)
w_prob[x > 7] <- 0.01
w <- rbinom(N_obs, 1, w_prob)
# w <- rbinom(N_obs, 1, 0.5)

treat_effect <- 0.75
y_obs <- rnorm(N_obs, mean= -2  + w*treat_effect, sd=0.2)

obs_df <- data.frame(x=x, y=y_obs, w=as.factor(w))

p <- ggplot(obs_df) + geom_point(aes(x=x, y=y, col=w))
# print(p)


#define a GP
lambda <- 10
sigma <- 1
rho <- sigma^2 * 0.5
k <- function(x1, x2, w1, w2){
  d <- sum((x1 - x2)^2)
  corr <- rho + (sigma^2 - rho) * (w1 == w2)
  return(corr * exp(- d / lambda^2))
}

K <- matrix(NA, N_obs, N_obs)
for(i in 1:N_obs){
  for(j in 1:N_obs){
    K[i, j] <- k(x[i], x[j], w[i], w[j])
  }
}

N_pred <- 200
x_pred <- seq(0, x_max, length.out = N_pred)
x_pred <- rep(x_pred, 2)
w_pred <- rep(0:1, each=N_pred)
N_pred <- N_pred * 2

K_pred <- matrix(NA, N_pred, N_obs)
for(i in 1:N_pred){
  for(j in 1:N_obs){
    K_pred[i, j] <- k(x_pred[i], x[j], w_pred[i], w[j])
  }
}

K_pred_pred <- matrix(NA, N_pred, N_pred)
for(i in 1:N_pred){
  for(j in 1:N_pred){
    K_pred_pred[i, j] <- k(x_pred[i], x_pred[j], w_pred[i], w_pred[j])
  }
}


lambda_GP <- 0.1
y_pred <- K_pred %*% solve(K + lambda_GP * diag(N_obs), y_obs)
y_var <- sqrt(diag(K_pred_pred - K_pred %*% solve(K + lambda_GP * diag(N_obs), t(K_pred))))


pred_df <- data.frame(x=x_pred, y=y_pred, w=as.factor(w_pred), ymin=y_pred - y_var,
                      ymax=y_pred+y_var)

p <- p + geom_ribbon(data=pred_df, aes(x=x, ymin=ymin, ymax=ymax, group=w), alpha=0.1) + 
  geom_line(data=pred_df, aes(x=x, y=y, col=w)) + theme_minimal()
  
print(p)


