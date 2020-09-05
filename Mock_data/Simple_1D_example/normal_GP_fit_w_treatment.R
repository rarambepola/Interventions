setwd(paste0(Sys.getenv("HOME"), "/Interventions/"))
library(ggplot2)

set.seed(2)

##mock data 
N_obs <- 100
x_max <- 13
x_min <- 0
x <- sort(runif(N_obs, x_min, x_max))
y_baseline <- rnorm(N_obs, mean=0, sd=0.5)
treat_effect <- -2
treat_vec <- rbinom(N_obs, 1, 0.5)
x_cutoff <- 7
treat_vec[x > x_cutoff] <- rbinom(sum(x > x_cutoff), 1, 0.05)

f <- function(x){
  return(0.25*x + sin(x/1.2) + 1)
}

y <- y_baseline + treat_effect * treat_vec + f(x)


##compute in R
#rbf kernel
lambda <- 4
sigma <- 1
sigma_obs <- 0.1
rho_scale <- 0.99
rho <- sigma^2 * rho_scale
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
x_pred <- rep(seq(x_min, x_max, length.out = N_pred), 2)
treat_pred <- rep(0:1, each=N_pred)
N_pred <- length(x_pred)

dist_mat_pred <- matrix(NA, N_pred, N_obs)
for(i in 1:N_pred){
  for(j in 1:N_obs){
    dist_mat_pred[i, j] <- sum((x_pred[i] - x[j])^2)
  }
}

treatment_same_mat <- matrix(0, N_obs, N_obs)
treatment_diff_mat <- matrix(0, N_obs, N_obs)

for(i in 1:N_obs){
  treatment_same_mat[i, treat_vec[i] == treat_vec] <- 1
  treatment_diff_mat[i, treat_vec[i] != treat_vec] <- 1
}

treatment_same_mat_pred <- matrix(0, N_pred, N_obs)
treatment_diff_mat_pred <- matrix(0, N_pred, N_obs)

for(i in 1:N_pred){
  treatment_same_mat_pred[i, treat_pred[i] == treat_vec] <- 1
  treatment_diff_mat_pred[i, treat_pred[i] != treat_vec] <- 1
}

K_pred <- k(dist_mat_pred, sigma, lambda)
K_use <- K * treatment_same_mat + (K * rho * treatment_diff_mat/ sigma^2)
K_pred_use <- K_pred * treatment_same_mat_pred + (K_pred * rho * treatment_diff_mat_pred/ sigma^2)

y_fit <- K_use %*% solve(K_use + sigma_obs^2 * diag(N_obs), y)
y_pred <- K_pred_use %*% solve(K_use + sigma_obs^2 * diag(N_obs), y)

plot(x, y, cex=0.001)
points(x[treat_vec == 0], y[treat_vec==0], col="red", pch=19, cex=0.5)
points(x[treat_vec == 1], y[treat_vec==1], col="blue", pch=19, cex=0.5)
lines(x_pred[treat_pred == 0], y_pred[treat_pred == 0], col="red")
lines(x_pred[treat_pred == 1], y_pred[treat_pred == 1], col="blue")
lines(x_pred[treat_pred == 0], f(x_pred[treat_pred == 0]), col="red", lty=3)
lines(x_pred[treat_pred == 1], f(x_pred[treat_pred == 1]) + treat_effect, col="blue", lty=3)


make_predictions <- function(dist_mat_pred, dist_mat, treat_pred, treat_vec, sigma, lambda, rho){
  K <- k(dist_mat, sigma, lambda)
  K_pred <- k(dist_mat_pred, sigma, lambda)
  # print(dim(K))
  K_use <- K * treatment_same_mat + (K * rho * treatment_diff_mat/ sigma^2)
  K_pred_use <- K_pred * treatment_same_mat_pred + (K_pred * rho * treatment_diff_mat_pred/ sigma^2)
  y_pred <- K_pred_use %*% solve(K_use + sigma_obs^2 * diag(N_obs), y)
  return(y_pred)
}

# 
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
treatment_same_mat_torch <- to_torch(treatment_same_mat)
treatment_diff_mat_torch <- to_torch(treatment_diff_mat)
y_torch <- to_torch_vector(y)
# 
#hyperparameters
log_sigma <- Variable(torch$ones(as.integer(1)), requires_grad=TRUE)
log_lambda <- Variable(torch$ones(as.integer(1))$add(-2), requires_grad=TRUE)
logit_rho_scale <- Variable(torch$zeros(as.integer(1)), requires_grad=TRUE)

define_priors <- function(log_sigma_mean,
                          log_sigma_sd,
                          log_lambda_mean,
                          log_lambda_sd,
                          logit_rho_scale_mean,
                          logit_rho_scale_sd){
  log_prior <- function(log_sigma,
                        log_lambda,
                        logit_rho_scale){
    log_prob_out <- torch$distributions$Normal(log_sigma_mean, log_sigma_sd)$log_prob(log_sigma)
    log_prob_out <- log_prob_out %+_t%
      torch$distributions$Normal(log_lambda_mean, log_lambda_sd)$log_prob(log_lambda)
    log_prob_out <- log_prob_out %+_t%
      torch$distributions$Normal(logit_rho_scale_mean, logit_rho_scale_sd)$log_prob(logit_rho_scale)
    return(log_prob_out)
  }
  return(log_prior)
}


log_prior <- define_priors(-3, 0.1, 4, 0.1, 10, 0.1)

optimizer <- torch$optim$Adam(list(log_sigma,
                                   log_lambda,
                                   logit_rho_scale),
                              lr=0.1)
n_iter <- 100
param_matrix <- matrix(NA, nrow=3, ncol=n_iter)
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
  prior_val <- log_prior(log_sigma, log_lambda, logit_rho_scale)
  # prior_val <- torch$zeros(as.integer(1))
  
  # nll <- nll %-_t% log_prior(log_sigma, log_lambda)

  #likelihood
  sigma2 <- torch$pow(torch$exp(log_sigma), 2)
  lambda2 <- torch$pow(torch$exp(log_lambda), 2)
  rho <- torch$mul(sigma2, torch$div(1, torch$add(1, torch$exp(torch$mul(-1, logit_rho_scale)))))
  # rho <- torch$zeros(as.integer(1))

  K_base <- torch$exp(torch$mul(-1, torch$div(dist_mat_torch, lambda2)))
  K_full <- torch$mul(sigma2, torch$mul(treatment_same_mat_torch, K_base)) %+_t%
    torch$mul(rho, torch$mul(treatment_diff_mat_torch, K_base))
  
  # print(K_full)
  
  cov_mat <- K_full %+_t% torch$diag(torch$mul(sigma_obs^2, torch$ones(as.integer(N_obs))))

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
  
  if(i < n_iter){
    optimizer$zero_grad()
    nll$backward(retain_graph=TRUE)
    optimizer$step()
  }
  # print(proc.time() - ptm)
  
  if(i %% 10 == 0){
    #plot fit
    pred <- make_predictions(dist_mat_pred, dist_mat, treat_pred, treat_vec, 
                             as.numeric(sqrt(sigma2$detach()$numpy())), 
                             as.numeric(sqrt(lambda2$detach()$numpy())), 
                             as.numeric(rho$detach()$numpy()))
    
    plot(x, y, cex=0.001, main=paste0("i: ", i, " likelihood: ", round(as.numeric(from_torch(likelihood_val)), 2)),
         sub=paste0("sigma: ", round(as.numeric(sqrt(sigma2$detach()$numpy())), 2),
                    " lambda: ", round(as.numeric(sqrt(lambda2$detach()$numpy())), 2),
                    " rho: ", round(as.numeric(rho$detach()$numpy()), 2)))
    points(x[treat_vec == 0], y[treat_vec==0], col="red", pch=19, cex=0.5)
    points(x[treat_vec == 1], y[treat_vec==1], col="blue", pch=19, cex=0.5)
    lines(x_pred[treat_pred == 0], pred[treat_pred == 0], col="red")
    lines(x_pred[treat_pred == 1], pred[treat_pred == 1], col="blue")
    lines(x_pred[treat_pred == 0], f(x_pred[treat_pred == 0]), col="red", lty=3)
    lines(x_pred[treat_pred == 1], f(x_pred[treat_pred == 1]) + treat_effect, col="blue", lty=3)
    
    
  }
}

lambda <- as.vector(exp(from_torch(log_lambda)))
sigma <- as.vector(exp(from_torch(log_sigma)))
rho <- as.vector(from_torch(rho))
# rho <- 0
# rho <- 0
# lambda <- 1


K <- k(dist_mat, sigma, lambda)
K_pred <- k(dist_mat_pred, sigma, lambda)
K_pred_use <- K_pred * treatment_same_mat_pred + (K_pred * rho * treatment_diff_mat_pred/ sigma^2)
K_use <- K * treatment_same_mat + (K * rho * treatment_diff_mat/ sigma^2)
y_pred <- K_pred_use %*% solve(K_use + sigma_obs^2 * diag(N_obs), y)

plot(x, y, cex=0.001)
points(x[treat_vec == 0], y[treat_vec==0], col="red", pch=19, cex=0.5)
points(x[treat_vec == 1], y[treat_vec==1], col="blue", pch=19, cex=0.5)
lines(x_pred[treat_pred == 0], y_pred[treat_pred == 0], col="red")
lines(x_pred[treat_pred == 1], y_pred[treat_pred == 1], col="blue")
lines(x_pred[treat_pred == 0], f(x_pred[treat_pred == 0]), col="red", lty=3)
lines(x_pred[treat_pred == 1], f(x_pred[treat_pred == 1]) + treat_effect, col="blue", lty=3)


cov_mat <- K_use + diag(N_obs) * sigma_obs^2
ll <- - 0.5 * N_obs * log(2*pi) - 0.5 * log(det(cov_mat)) - 0.5 * (t(y) %*% solve(cov_mat, y))
cat("\n")
print(as.vector(ll))
print(likelihood_vec[n_iter])
print(- 0.5 * log(det(cov_mat)) - 0.5 * (t(y) %*% solve(cov_mat, y)))
print(log(det(cov_mat)))

cov_mat_old <- cov_mat
K_use_old <- K_use

 
# #plot fit
# lambda <- as.vector(exp(from_torch(log_lambda)))
# sigma <- as.vector(exp(from_torch(log_sigma)))
# 
# K_pred <- k(dist_mat_pred, sigma, lambda)
# K <- k(dist_mat, sigma, lambda)
# y_fit <- K %*% solve(K + sigma_obs^2 * diag(N_obs), y)
# y_pred <- K_pred %*% solve(K + sigma_obs^2 * diag(N_obs), y)
# 
# plot(x, y)
# lines(x_pred, y_pred, col="red")
# 
# 
# 
# 
# 
