# args <- commandArgs(trailingOnly=TRUE)
args <- c("no_selection_bias", 1)

#check if using google cloud
#check whether on google cloud or not
google_cloud <- (Sys.getenv("HOME") == "/root")

#set up directories
if(google_cloud){
  print("run_model.R")
  path_output <- Sys.getenv("path_output")
  path_input <- Sys.getenv("path_input")
  setwd(path_input)
  print(getwd())
  output_folder <- paste0(path_output, "/")
}else{
  setwd(paste0(Sys.getenv("HOME"), "/Interventions"))
  output_folder <- NULL
}

if(google_cloud){
  install.packages("randtoolbox")
  install.packages("reticulate")
  reticulate::use_python("/usr/bin/python3", required=T)
}else{
  reticulate::use_python("/usr/bin/python", required=T)
}

library(reticulate)
print(py_config())

library(ggplot2)
library(raster)
library(randtoolbox)

mock_data_scenario <- args[1]

# load("Mock_data/Mdg/elevation_selection_bias.RData")
load(paste0("Mock_data/Mdg/", mock_data_scenario, ".RData"))
# load("Mock_data/Mdg/elevation_selection_bias_nonlinear.RData")

#change inputs into covariance matrix and weights vector
X <- covs[, 1:3]
W <- covs[, 4]
N <- cluster_size
N_pos <- cluster_pos
p_obs <- N_pos / N
N_obs <- dim(X)[1]


emp_logit_mean <- function(p, n){
  return(log((p + 1/(2*n)) / (1 - p + 1/(2*n))))
}
#convert observations using empircal logit
Y <- emp_logit_mean(N_pos/N, N)
Y_var <- (1 / (N_pos + 0.5)) + (1 / (N - N_pos + 0.5))
#convert true baseline and treatment prevalence too
Y_baseline <- emp_logit_mean(environmental_prevalence, N)
Y_treated <- emp_logit_mean(with_treatment_prevalence, N)
Y_treatment_effect <- Y_treated - Y_baseline

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
nt_nt_ind <- matrix(0, n_obs, n_obs)
t_nt_ind <- matrix(0, n_obs, n_obs)
nt_t_ind <- matrix(0, n_obs, n_obs)

for(i in 1:n_obs){
  for(j in 1:n_obs){
    if(W[i] == 1){
      if(W[j] == 1){
        t_t_ind[i, j] <- 1
      }else{
        t_nt_ind[i, j] <- 1
      }
    }else{
      if(W[j] == 1){
        nt_t_ind[i, j] <- 1
      }else{
        nt_nt_ind[i, j] <- 1
      }
    }
  }
}


#create weights for RFFs
N_covs <- dim(X)[2]
N_rff <- 100
use_RFF <- TRUE

#weights have sd 1 to start with
w <- matrix(qnorm(halton(N_covs*N_rff)), nrow=N_covs)

#treatment index matrices
any_index_mat <- matrix(0, nrow=2*N_rff*3, ncol=N_obs)
treat_index_mat <- any_index_mat
not_treat_index_mat <- any_index_mat
treat_index_mat_cf <- any_index_mat
not_treat_index_mat_cf <- any_index_mat
any_index_mat[1:(2*N_rff), ] <- 1
treat_index_mat[2*N_rff + (1:(2*N_rff)), W==1] <- 1
treat_index_mat_cf[2*N_rff + (1:(2*N_rff)), W==0] <- 1
not_treat_index_mat[4*N_rff + (1:(2*N_rff)), W==0] <- 1
not_treat_index_mat_cf[4*N_rff + (1:(2*N_rff)), W==1] <- 1


np <- import("numpy")
torch <- import("torch")
Variable <-  torch$autograd$Variable

log_beta00 <- Variable(torch$ones(as.integer(1))$add(1), requires_grad=TRUE)
log_beta01 <- Variable(torch$ones(as.integer(1))$add(0.5), requires_grad=TRUE)
logit_rho_scale <- Variable(torch$ones(as.integer(1))$add(-3), requires_grad=TRUE)
log_lambda <- Variable(torch$ones(as.integer(1)), requires_grad=TRUE)
log_sigma_0 <- Variable(torch$ones(as.integer(1)), requires_grad=TRUE)
log_sigma_1 <- Variable(torch$ones(as.integer(1)), requires_grad=TRUE)
# sigmas <- c(sigma_0, sigma_1)


W_inv <- 1 - W
W <- torch$from_numpy(matrix(W))$float()$squeeze()
W_inv <- torch$from_numpy(matrix(W_inv))$float()$squeeze()
Y_var <- torch$from_numpy(matrix(Y_var))$float()$squeeze()
Y <- torch$from_numpy(matrix(Y))$float()$squeeze()
w <- torch$from_numpy(w)$float()
X <- torch$from_numpy(X)$float()

any_index_mat <- torch$from_numpy(any_index_mat)$float()
treat_index_mat <- torch$from_numpy(treat_index_mat)$float()
not_treat_index_mat <- torch$from_numpy(not_treat_index_mat)$float()
treat_index_mat_cf <- torch$from_numpy(treat_index_mat_cf)$float()
not_treat_index_mat_cf <- torch$from_numpy(not_treat_index_mat_cf)$float()

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
# n_iter <- 2
n_iter <- args[2]

n_obs_torch <- torch$from_numpy(matrix(n_obs))$float()$squeeze()

'%+_t%' <- function(x, y) torch$add(x, y)
'%-_t%' <- function(x, y) torch$sub(x, y)
t_torch <- function(X) X$transpose(as.integer(1), as.integer(0))

#function to compute variance of counterfactuals
compute_var <- function(K_cf_cf, K_obs_cf, K_obs_obs, Sigma_obs){
  Var_cf <- torch$trace(K_cf_cf %-_t% torch$matmul(K_obs_cf, torch$solve(t_torch(K_obs_cf),
                                                                         K_obs_obs %+_t% Sigma_obs)[0]))
  return(Var_cf)
}

compute_var_Rff <- function(X, y, Sigma_inv, X_star, N_features){
  A <- torch$matmul(torch$matmul(X, Sigma_inv), t_torch(X)) %+_t% torch$diag(torch$ones(as.integer(N_features)))
  return(torch$trace(torch$matmul(t_torch(X_star), torch$solve(X_star, A)[0])))
}

#function to compute leave one out prediction
compute_LOO_pred <- function(j, K_obs_obs, Sigma_obs, Y){
  minus_j_index <- setdiff((0:(n_obs-1)), j-1)
  K_12 <- t_torch(K_obs_obs[minus_j_index])[j-1]
  K_22 <- t_torch(K_obs_obs[minus_j_index])[minus_j_index]
  Sigma_22 <- K_22 %+_t% t_torch(Sigma_obs[minus_j_index])[minus_j_index]
  Y_fit_j <- torch$matmul(K_12$unsqueeze(as.integer(0)), torch$solve(Y[minus_j_index]$unsqueeze(as.integer(1)), Sigma_22)[0])
  return(Y_fit_j)
}

# A <- X_rff %*% sigma_obs_inv %*% t(X_rff) + diag(N_rff * 2 * N_covs)
# w_bar_1 <- solve(A, X_rff) %*% sigma_obs_inv %*% Y
# mean_cf_rff <- t(X_rff_cf) %*% w_bar_1

compute_LOO_pred_rff <- function(j, X, Sigma_inv, Y){
  minus_j_index <- setdiff((0:(n_obs-1)), j-1)
  X_use <- t_torch(t_torch(X)[minus_j_index])
  Sigma_inv <- t_torch(Sigma_inv[minus_j_index])[minus_j_index]
  X_star <- t_torch(X)[j-1]
  Y <- Y[minus_j_index]
  A <- torch$matmul(torch$matmul(X_use, Sigma_inv), t_torch(X_use)) %+_t% torch$diag(torch$ones(as.integer(N_rff * 2 * N_covs)))
  w_bar <- torch$matmul(torch$matmul(torch$solve(X_use, A)[0], Sigma_inv), Y)
  mean_pred <- torch$mul(w_bar, X_star)$sum()
  return(mean_pred)
}


for (i in 1:n_iter) {
  print(i)
  ptm <- proc.time()
  K <- torch$exp(torch$mul(-1, torch$div(dist_mat, torch$mul(0.5, torch$exp(log_lambda)))))
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
  
  if(use_RFF){
    #make RFFs
    lambda_sq <- torch$pow(lambda, 0.5)
    treat_coeff <- torch$pow(torch$sub(beta012, rho), 0.5)
    not_treat_coeff <- torch$pow(torch$sub(beta002, rho), 0.5)
    w_use <- torch$div(w, lambda_sq)
    wX <- t_torch(torch$matmul(X, w_use))
    X_rff_base <- torch$div(torch$cat(c(torch$sin(wX), torch$cos(wX)), as.integer(0)), sqrt(N_rff))
    X_rff_base <- torch$cat(c(X_rff_base, X_rff_base, X_rff_base), as.integer(0))
    X_rff <- torch$mul(torch$mul(X_rff_base, any_index_mat), rho) %+_t% 
      torch$mul(torch$mul(X_rff_base, treat_index_mat), treat_coeff) %+_t%
      torch$mul(torch$mul(X_rff_base, not_treat_index_mat), not_treat_coeff)
    
    X_rff_cf <- torch$mul(torch$mul(X_rff_base, any_index_mat), rho) %+_t% 
      torch$mul(torch$mul(X_rff_base, treat_index_mat_cf), treat_coeff) %+_t%
      torch$mul(torch$mul(X_rff_base, not_treat_index_mat_cf), not_treat_coeff)
    
    Sigma_obs_inv <- torch$diag(torch$div(1, torch$mul(torch$add(torch$mul(W_inv, sigma_1), torch$mul(W, sigma_0)), Y_var)))
    
    Var_cf <- compute_var_Rff(X_rff, y, Sigma_obs_inv, X_rff_cf, N_rff*2*N_covs)
    # K_cf_cf <- torch$mul(K, torch$pow(beta01, 2))
  }else{
    K_obs_obs <- torch$mul(K,
                           (torch$mul(t_t_ind, beta012) %+_t% torch$mul(nt_nt_ind, beta002)) %+_t% 
                             (torch$mul(t_nt_ind, rho) %+_t% torch$mul(nt_t_ind, rho))
    )
    
    K_cf_cf <- torch$mul(K,
                         (torch$mul(t_t_ind, beta002) %+_t% torch$mul(nt_nt_ind, beta012)) %+_t% 
                           (torch$mul(t_nt_ind, rho) %+_t% torch$mul(nt_t_ind, rho))
    )
    
    K_obs_cf <- torch$mul(K,
                          (torch$mul(t_t_ind, rho) %+_t% torch$mul(nt_nt_ind, rho)) %+_t% 
                            (torch$mul(t_nt_ind, beta012) %+_t% torch$mul(nt_t_ind, beta002))
    )
    
    Sigma_cf <- torch$diag(torch$mul(torch$add(torch$mul(W_inv, sigma_0), torch$mul(W, sigma_1)), Y_var))
    Sigma_obs <- torch$diag(torch$mul(torch$add(torch$mul(W_inv, sigma_1), torch$mul(W, sigma_0)), Y_var))
    
    Var_cf <- compute_var(K_cf_cf, t_torch(K_obs_cf), K_obs_obs, Sigma_obs)
  }

  

  # K_obs_obs <- torch$mul(K, torch$pow(beta00, 2))
  # K_obs_cf <- torch$mul(K, rho)

 

  # print(Var_cf)
  # print(Var_cf_rff)
  
  dist_sum_torch <- torch$zeros(as.integer(1))
  # dist_sum_torch_rff <- torch$zeros(as.integer(1))
  # for(j in 1:n_obs){
  for(j in 1:10){
    if(j %% 10 == 0){
      cat(j)
      cat(" ")
    }
    if(use_RFF){
      Y_fit_j <- compute_LOO_pred_rff(j, X_rff, Sigma_obs_inv, Y)
    }else{
      Y_fit_j <- compute_LOO_pred(j, K_obs_obs, Sigma_obs, Y)
    }
    dist_sum_torch <- dist_sum_torch$add(torch$pow(torch$add(torch$mul(-1, Y[j-1]), Y_fit_j), 2))
    # dist_sum_torch_rff <- dist_sum_torch_rff$add(torch$pow(torch$add(torch$mul(-1, Y[j-1]), Y_fit_j_rff), 2))
  }
  cat("\n")
  
  # print("dist sum")
  # print(dist_sum_torch)
  # print(dist_sum_torch_rff)
  
  param_size <- beta002 %+_t% beta012 %+_t% sigma_0$pow(2) %+_t% 
    sigma_1$pow(2) %+_t% rho$pow(2) %+_t% torch$mul(lambda$pow(2), 3)
  
  #10 times number of hyperparameters = 60... I think
  loss <- torch$mul(torch$log(Var_cf %+_t% dist_sum_torch), n_obs_torch) %+_t%
    torch$mul(60, torch$log(param_size))
  
  loss_vec[i] <- loss$detach()$numpy()
  #dont take a step the last time
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
Y_pred_treatment_effect <- Y_pred_treated - Y_pred_not_treated

# plot(X[, 3], Y_pred_not_treated, cex=0.5)
# points(X[, 3], Y_pred_treated, cex=0.5, col="red")
# 
# print(cor(environmental_risk, Y_pred_not_treated))
# print(cor(true_risk, Y_pred_not_treated))
# print(cor(true_risk, Y_fit))
# 
# plot(Y_pred_treated, Y_pred_not_treated)
# abline(0, 1)
# 

# 
# plot(Y_treatment_effect, Y_pred_treatment_effect)

save(list=c("beta00", "beta01", "lambda", "sigma_0", "sigma_1", "rho", "n_iter",
            "Y_pred_treated", "Y_pred_not_treated", "Y_treatment_effect", "Y_pred_treatment_effect",
            "Y_fit", "param_matrix", "loss_vec"),
file=paste0(output_folder, "Mock_data/Mdg/", mock_data_scenario, "fit_", n_iter, ".RData"))


