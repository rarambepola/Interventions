setwd(paste0(Sys.getenv("HOME"), "/Interventions"))
library(ggplot2)
library(raster)
library(INLA)
library(TMB)

load("Mock_data/Mdg/elevation_selection_bias_nonlinear.RData")
# load("Mock_data/Mdg/no_selection_bias.RData")

X <- covs[, 1:3]
# X <- covs[, 1, drop=FALSE]
W <- covs[, 4]
N <- cluster_size
N_pos <- cluster_pos
p_obs <- N_pos / N

Y <- log((N_pos + 0.5) / (N - N_pos + 0.5))
Y_var <- (1 / (N_pos + 0.5)) + (1 / (N - N_pos + 0.5))
# Y_var <- 1

#create distance matrix
n_obs <- length(N)
dist_mat <- matrix(NA, n_obs, n_obs)

for(i in 1:n_obs){
  for(j in 1:i){
    d <- sum((X[i, ] - X[j, ])^2)
    dist_mat[i, j] <- d
    dist_mat[j, i] <- d
  }
}
dist_mat2 <- dist_mat
# #hyperparameters are variables
t_t_ind <- matrix(0, n_obs, n_obs)
t_t_ind[W==1, W==1] <- 1
nt_nt_ind <- matrix(0, n_obs, n_obs)
nt_nt_ind[W==0, W==0] <- 1
t_nt_ind <- matrix(0, n_obs, n_obs)
t_nt_ind[W==1, W==0] <- 1
nt_t_ind <- matrix(0, n_obs, n_obs)
nt_t_ind[W==0, W==1] <- 1

reticulate::use_python("/usr/bin/python", required=T)
library(reticulate)
print(py_config())

np <- import("numpy")
torch <- import("torch")
Variable <-  torch$autograd$Variable

log_beta00 <- Variable(torch$ones(as.integer(1))$add(-1), requires_grad=TRUE)
log_beta01 <- Variable(torch$ones(as.integer(1))$add(-1), requires_grad=TRUE)
logit_rho_scale <- Variable(torch$ones(as.integer(1))$add(1), requires_grad=TRUE)
log_lambda <- Variable(torch$ones(as.integer(1)), requires_grad=TRUE)
log_sigma_0 <- Variable(torch$ones(as.integer(1)), requires_grad=TRUE)
log_sigma_1 <- Variable(torch$ones(as.integer(1)), requires_grad=TRUE)
# sigmas <- c(sigma_0, sigma_1)


W_inv <- 1 - W
W <- torch$from_numpy(matrix(W))$float()$squeeze()
W_inv <- torch$from_numpy(matrix(W_inv))$float()$squeeze()
Y_var <- torch$from_numpy(matrix(Y_var))$float()$squeeze()
Y <- torch$from_numpy(matrix(Y))$float()$squeeze()

t_t_ind <- torch$from_numpy(t_t_ind)$float()
nt_nt_ind <- torch$from_numpy(nt_nt_ind)$float()
t_nt_ind <- torch$from_numpy(t_nt_ind)$float()
nt_t_ind <- torch$from_numpy(nt_t_ind)$float()


optimizer <- torch$optim$Adam(list(log_beta00,
                                   log_beta01,
                                   logit_rho_scale,
                                   log_lambda,
                                   log_sigma_0,
                                   log_sigma_1),
                              lr=0.1)

params <- c(log_beta00, log_beta01, logit_rho_scale, log_lambda, log_sigma_0, log_sigma_1)
dist_mat <- torch$from_numpy(dist_mat)$float()
dist_mat3 <- dist_mat$detach()$numpy()
loss_vec <- c()
param_matrix <- c()
n_iter <- 20


'%+_t%' <- function(x, y) torch$add(x, y)


for (i in 1:n_iter) {
  print(i)
  ptm <- proc.time()
  K <- torch$exp(torch$mul(-1, torch$div(dist_mat, torch$exp(log_lambda))))
  param_matrix <- cbind(param_matrix, sapply(params, function(p) p$detach()$numpy()))
  
  beta00 <- torch$exp(log_beta00)
  beta01 <- torch$exp(log_beta01)
  lambda <- torch$exp(log_lambda)
  sigma_0 <- torch$exp(log_sigma_0)
  sigma_1 <- torch$exp(log_sigma_1)
  
  rho_scale <- torch$div(1, torch$add(1, torch$exp(torch$mul(-1, logit_rho_scale))))
  rho <- torch$mul(beta00, torch$mul(rho_scale, beta01))
  
  beta012 <- torch$pow(beta01, 2)
  beta002 <- torch$pow(beta00, 2)
  
  K_obs_obs <- torch$mul(K,
                         (torch$mul(t_t_ind, beta002) %+_t% torch$mul(nt_nt_ind, beta012)) %+_t% 
                           (torch$mul(t_nt_ind, rho) %+_t% torch$mul(nt_t_ind, rho))
  )
  
  K_cf_cf <- torch$mul(K,
                       (torch$mul(t_t_ind, beta012) %+_t% torch$mul(nt_nt_ind, beta002)) %+_t% 
                         (torch$mul(t_nt_ind, rho) %+_t% torch$mul(nt_t_ind, rho))
  )
  
  K_obs_cf <- torch$mul(K,
                        (torch$mul(t_t_ind, rho) %+_t% torch$mul(nt_nt_ind, rho)) %+_t% 
                          (torch$mul(t_nt_ind, beta012) %+_t% torch$mul(nt_t_ind, beta002))
  )
  
  # K_cf_cf <- torch$mul(K, torch$pow(beta01, 2))
  # K_obs_obs <- torch$mul(K, torch$pow(beta00, 2))
  # K_obs_cf <- torch$mul(K, rho)
  Sigma_cf <- torch$diag(torch$mul(torch$add(torch$mul(W_inv, sigma_0), torch$mul(W, sigma_1)), Y_var))
  Sigma_obs <- torch$diag(torch$mul(torch$add(torch$mul(W_inv, sigma_1), torch$mul(W, sigma_0)), Y_var))
  Var_cf <- torch$trace(
    torch$add(
      K_cf_cf,
        torch$mul(
          -1,
          torch$matmul(K_obs_cf,
                       torch$solve(
                         K_obs_cf$transpose(as.integer(1), as.integer(0)),
                         torch$add(K_obs_obs, Sigma_obs)
                       )[0]
                       )
        )
    )
  )
  
  dist_sum_torch <- torch$zeros(as.integer(1))
  
  for(j in 1:n_obs){
  # for(j in 1:1){
    cat(j)
    cat(" ")
    minus_i_index <- setdiff((0:(n_obs-1)), j-1)
    K_12 <- K_obs_obs[minus_i_index]$transpose(as.integer(1), as.integer(0))[j-1]
    K_22 <- K_obs_obs[minus_i_index]$transpose(as.integer(1), as.integer(0))[minus_i_index]
    Sigma_22 <- torch$add(K_22, Sigma_obs[minus_i_index]$transpose(as.integer(1), as.integer(0))[minus_i_index])
    Y_fit_i <- torch$matmul(K_12$unsqueeze(as.integer(0)), torch$solve(Y[minus_i_index]$unsqueeze(as.integer(1)), Sigma_22)[0])
    dist_sum_torch <- dist_sum_torch$add(torch$pow(torch$add(torch$mul(-1, Y[j-1]), Y_fit_i), 2))
    # Y_fit[i] <- Y_fit_i$detach()$numpy()
  }
  cat("\n")
  
  loss <- torch$add(Var_cf, dist_sum_torch)
  
  loss_vec[i] <- loss$detach()$numpy()
  if(i < n_iter){
    optimizer$zero_grad()
    loss$backward(retain_graph=TRUE)
    optimizer$step()
  }
  print(proc.time() - ptm)
}

dist_mat_old <- dist_mat
K_old <- K$detach()$numpy()
params_r <- sapply(params, function(p) p$detach()$numpy())




#now make predictions
Sigma_22 <- torch$add(K_obs_obs, Sigma_obs)
proj_vec <- as.vector(torch$solve(Y$unsqueeze(as.integer(1)), Sigma_22)[0]$detach()$numpy())
K <- K$detach()$numpy()
beta00 <- as.vector(beta00$detach()$numpy())
beta01 <- as.vector(beta01$detach()$numpy())
rho <- as.vector(rho$detach()$numpy())
W <- as.vector(W$detach()$numpy())
sigma_0 <- as.vector(sigma_0$detach()$numpy())
sigma_1 <- as.vector(sigma_1$detach()$numpy())
lambda <- as.vector(lambda$detach()$numpy())

K_obs_treated <- K
K_obs_treated[, W==1] <- K_obs_treated[, W==1] * beta00^2
K_obs_treated[, W==0] <- K_obs_treated[, W==0] * rho

K_obs_not_treated <- K
K_obs_not_treated[, W==1] <- K_obs_not_treated[, W==1] * rho
K_obs_not_treated[, W==0] <- K_obs_not_treated[, W==0] * beta01^2

Y_pred_treated <- as.vector(K_obs_treated %*% proj_vec)
Y_pred_not_treated <- as.vector(K_obs_not_treated %*% proj_vec)
Y_fit <- as.vector(K_obs_obs$detach()$numpy() %*% proj_vec)

plot(X[, 3], Y_pred_not_treated, cex=0.5)
points(X[, 3], Y_pred_treated, cex=0.5, col="red")

print(cor(environmental_risk, Y_pred_not_treated))
print(cor(true_risk, Y_pred_not_treated))
print(cor(true_risk, Y_fit))

plot(Y_pred_treated, Y_pred_not_treated)
abline(0, 1)

# save(list=c("beta00", "beta01", "lambda", "sigma_0", "sigma_1", "rho", "n_iter"),
     # file="Mock_data/Mdg/no_selection_bias_fit.RData")

save(list=c("beta00", "beta01", "lambda", "sigma_0", "sigma_1", "rho", "n_iter"),
     file="Mock_data/Mdg/elevation_selection_bias_nonlinear_fit.RData")

# #change into normal parameters
# params_untransformed <- rep(NA, 6)
# params_untransformed[c(1:2, 4:6)] <- exp(param_matrix[c(1:2, 4:6), n_iter])
# # rho_scale <- 1 / (1 + exp(-param_matrix[3, n_iter]))
# # rho <- params_untransformed[1] * params_untransformed[2] * rho_scale
# params_untransformed[3] <- rho
# params_r <- params_untransformed
# print(loss)
# 


# Y_W <- torch$matmul(K_obs_obs, proj_vec)
# Y_n_W <- torch$matmul(K_obs_cf, proj_vec)
# Y_W <- Y_W$detach()$numpy()
# Y_n_W <- Y_n_W$detach()$numpy()
# 
# W <- W$detach()$numpy()
# Y_test <- torch$matmul(K, proj_vec)
# 
# 
# plot(X[W==1], Y_W[W==1], cex=0.5)
# points(X[W==0], Y_W[W==0], cex=0.5, col="red")
# 
# 
# plot(X, Y_W, cex=0.5)
# points(X, Y_n_W, cex=0.5, col="green")
# # points(X, Y_test$detach()$numpy(), cex=0.5, col="blue")
# plot(X, Y_test$detach()$numpy(), cex=0.5, col="blue")
# 
# Y <- Y$detach()$numpy()
#
# #can do this in R
# K_obs_obs <- K_obs_obs$detach()$numpy()
# K_cf
#
# 
# 
# Y_fit <- c()
# for(j in 1:n_obs){
#   cat(j)
#   cat(" ")
#   minus_i_index <- setdiff((0:(n_obs-1)), j-1)
#   K_12 <- K[minus_i_index]$transpose(as.integer(1), as.integer(0))[j-1]
#   K_22 <- K[minus_i_index]$transpose(as.integer(1), as.integer(0))[minus_i_index]
#   Sigma_22 <- torch$add(K_22, Sigma_obs[minus_i_index]$transpose(as.integer(1), as.integer(0))[minus_i_index])
#   Y_fit_i <- torch$matmul(K_12$unsqueeze(as.integer(0)), torch$solve(Y[minus_i_index]$unsqueeze(as.integer(1)), Sigma_22)[0])
#   dist_sum_torch <- dist_sum_torch$add(torch$pow(torch$add(torch$mul(-1, Y[j-1]), Y_fit_i), 2))
#   Y_fit[j] <- Y_fit_i$detach()$numpy()
# }

# 
# Y <- Y$detach()$numpy()
# W <- W$detach()$numpy()
# plot(Y, Y_fit)
# cor(Y, Y_fit)
# 
# # plot(X, Y)
# # plot(X, Y_fit)
# i <- 1
# plot(X[W == 1, i], Y_fit[W==1], cex=0.5)
# points(X[W == 0, i], Y_fit[W==0], cex=0.5, col="red")
# 
# 
# # g_inv <- function(x) (1.5*exp(x) - 0.5) / (1 + exp(x))
# g_inv <- function(y, n) return((exp(y) - exp(y)/(2*n) - 1/(2*n)) / (1 + exp(y)))
# Y_inv <- 1 / (1 + exp(-Y_fit))
# Y_inv2 <- g_inv(Y_fit, N)
# print(mean(Y_inv[W==1]))
# print(mean(Y_inv[W==0]))
# print(mean(Y_inv2[W==1]))
# print(mean(Y_inv2[W==0]))
# 


# #compute GP surfaces
# Sigma_22 <- torch$add(K_obs_obs, Sigma_obs)
# proj_vec <- torch$solve(Y$unsqueeze(as.integer(1)), Sigma_22)[0]


# # ##compute LOO-CV
# Y_fit <- c()
# # i <- 3
# dist_sum_torch <- torch$zeros(as.integer(1))
# for(i in 1:n_obs){
#   print(i)
#   minus_i_index <- setdiff((0:(n_obs-1)), i-1)
#   K_12 <- K_obs_obs[minus_i_index]$transpose(as.integer(1), as.integer(0))[i-1]
#   K_22 <- K_obs_obs[minus_i_index]$transpose(as.integer(1), as.integer(0))[minus_i_index]
#   Sigma_22 <- torch$add(K_22, Sigma_obs[minus_i_index]$transpose(as.integer(1), as.integer(0))[minus_i_index])
#   Y_fit_i <- torch$matmul(K_12$unsqueeze(as.integer(0)), torch$solve(Y[minus_i_index]$unsqueeze(as.integer(1)), Sigma_22)[0])
#   dist_sum_torch <- dist_sum_torch$add(torch$pow(torch$add(torch$mul(-1, Y[i-1]), Y_fit_i), 2))
#   Y_fit[i] <- Y_fit_i$detach()$numpy()
# }
# # # K_12_R <- K[-i, i]
# # # Y_fit_R[i] <- K_12_R %*% solve(K[-i, -i] + sigma_obs[-i, -i], Y[-i])
# # 
# # 
# # 
