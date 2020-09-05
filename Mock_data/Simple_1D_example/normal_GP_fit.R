setwd(paste0(Sys.getenv("HOME"), "/Interventions/"))
library(ggplot2)

set.seed(2)

##mock data 
N_obs <- 100
x_max <- 10
x_min <- 0
x <- sort(runif(N_obs, x_min, x_max))
y <- rnorm(N_obs, mean=0, sd=0.2)

##compute in R
#rbf kernel
lambda <- 2
sigma <- 1
sigma_obs <- 0.1
dist_mat <- matrix(NA, N_obs, N_obs)
for(i in 1:N_obs){
  for(j in 1:N_obs){
    dist_mat[i, j] <- sum((x[i] - x[j])^2)
  }
}

k <- function(d2, sigma=sigma, lambda=lambda){
  return(sigma^2 * exp(-d2 / (lambda^2)))
}

K <- k(dist_mat, sigma, lambda)

N_pred <- 100
x_pred <- seq(x_min, x_max, length.out = N_pred)

dist_mat_pred <- matrix(NA, N_pred, N_obs)
for(i in 1:N_pred){
  for(j in 1:N_obs){
    dist_mat_pred[i, j] <- sum((x_pred[i] - x[j])^2)
  }
}

K_pred <- k(dist_mat_pred, sigma, lambda)

y_fit <- K %*% solve(K + sigma_obs^2 * diag(N_obs), y)
y_pred <- K_pred %*% solve(K + sigma_obs^2 * diag(N_obs), y)

plot(x, y)
lines(x_pred, y_pred, col="red")


##fit in torch
reticulate::use_python("/usr/bin/python", required=T)
library(reticulate)
print(py_config())
np <- import("numpy")
torch <- import("torch")
Variable <-  torch$autograd$Variable

to_torch <- function(m){
  return(torch$from_numpy(m)$float())
}

to_torch_vector <- function(m){
  return(torch$from_numpy(matrix(m))$float()$squeeze())
}

from_torch <- function(x){
  x$detach()$numpy()
}

'%+_t%' <- function(x, y) torch$add(x, y)
'%-_t%' <- function(x, y) torch$sub(x, y)
t_torch <- function(X) X$transpose(as.integer(1), as.integer(0))

dist_mat_torch <- to_torch(dist_mat)
y_torch <- to_torch_vector(y)

#hyperparameters
log_sigma <- Variable(torch$ones(as.integer(1)), requires_grad=TRUE)
log_lambda <- Variable(torch$ones(as.integer(1)), requires_grad=TRUE)

define_priors <- function(log_sigma_mean, 
                          log_sigma_sd,
                          log_lambda_mean,
                          log_lambda_sd){
  log_prior <- function(log_sigma, 
                        log_lambda){
    log_prob_out <- torch$distributions$Normal(log_sigma_mean, log_sigma_sd)$log_prob(log_sigma)
    log_prob_out <- log_prob_out %+_t% 
      torch$distributions$Normal(log_lambda_mean, log_lambda_sd)$log_prob(log_lambda)
    return(log_prob_out)
  }
  return(log_prior)
}


log_prior <- define_priors(0, 0.5, 2, 0.5)

optimizer <- torch$optim$Adam(list(log_sigma,
                                   log_lambda),
                              lr=0.1)
n_iter <- 100
param_matrix <- matrix(NA, nrow=2, ncol=n_iter)
prior_vec <- c()
likelihood_vec <- c()
loss_vec <- c()
for(i in 1:n_iter){
  if(i %% 10 == 0){
    cat(i)
    cat(" ")
  }

  ptm <- proc.time()
  param_matrix[1, i] <- from_torch(log_sigma)
  param_matrix[2, i] <- from_torch(log_lambda)
  
  #prior
  nll <- 0
  prior_val <- log_prior(log_sigma, log_lambda)
  # priot_val <- torch$zeros(as.integer(1))

  
  # nll <- nll %-_t% log_prior(log_sigma, log_lambda)
  
  #likelihood
  sigma2 <- torch$pow(torch$exp(log_sigma), 2)
  lambda2 <- torch$pow(torch$exp(log_lambda), 2)
  
  K <- torch$mul(sigma2, torch$exp(torch$mul(-1, torch$div(dist_mat_torch, lambda2))))
  cov_mat <- K %+_t% torch$diag(torch$mul(sigma_obs^2, torch$ones(as.integer(N_obs))))
  
  likelihood_val <- torch$distributions$MultivariateNormal(torch$zeros(as.integer(N_obs)),
                                                              cov_mat)$log_prob(y_torch)
  
  # nll <- nll %-_t% torch$distributions$MultivariateNormal(torch$zeros(as.integer(N_obs)),
  #                                                         cov_mat)$log_prob(y_torch)
  
  # y_fit <- torch$matmul(K, torch$solve(y_torch$unsqueeze(as.integer(1)), cov_mat)[[0]])
  
  loss_val <- torch$mul(-1, prior_val %+_t% likelihood_val)
  nll <- loss_val
  
  loss_vec[i] <- from_torch(loss_val)
  prior_vec[i] <- from_torch(prior_val)
  likelihood_vec[i] <- from_torch(likelihood_val)
  optimizer$zero_grad()
  nll$backward(retain_graph=TRUE)
  optimizer$step()
  # print(proc.time() - ptm)
}


#plot fit
lambda <- as.vector(exp(from_torch(log_lambda)))
sigma <- as.vector(exp(from_torch(log_sigma)))

K_pred <- k(dist_mat_pred, sigma, lambda)
K <- k(dist_mat, sigma, lambda)
y_fit <- K %*% solve(K + sigma_obs^2 * diag(N_obs), y)
y_pred <- K_pred %*% solve(K + sigma_obs^2 * diag(N_obs), y)

plot(x, y)
lines(x_pred, y_pred, col="red")





