setwd(paste0(Sys.getenv("HOME"), "/Interventions"))
load("Mock_data/Mdg/elevation_selection_bias.RData")

#fix some hyperparameter vals and compute risk
beta00 <- 3
beta01 <- 0.5
# rho <- 0.3
rho_scale <- 0.2
lambda <- 2
sigma_0 <- 3
sigma_1 <- 1
sigmas <- c(sigma_0, sigma_1)

# beta00 <- 1
# beta01 <- 1
# rho <- 1
# lambda <- 1
# sigma_0 <- 1
# sigma_1 <- 1

# beta00 <- params_r[1]
# beta01 <- params_r[2]
# rho <- params_r[3]
# lambda <- params_r[4]
# sigma_0 <- params_r[5]
# sigma_1 <- params_r[6]

# rho_scale <- 1 / (1 + exp(-params_r[3]))
rho <- sqrt(beta00 * beta01) * rho_scale

sigmas <- c(sigma_0, sigma_1)



X <- covs[, 1:3]
W <- covs[, 4]
# W <- c(1, 1, 0)
N <- cluster_size
N_pos <- cluster_pos
N_covs <- dim(X)[2]
N_obs <- length(W)

Y <- log((N_pos + 0.5) / (N - N_pos + 0.5))
Y_var <- (1 / (N_pos + 0.5)) + (1 / (N - N_pos + 0.5))
# Y_var <- 1
# sigmas <- c(1,1)

#create distance matrix
n_obs <- length(W)
dist_mat <- matrix(NA, n_obs, n_obs)

for(i in 1:n_obs){
  for(j in 1:i){
    d <- sum((X[i, ] - X[j, ])^2)
    dist_mat[i, j] <- d
    dist_mat[j, i] <- d
  }
}

#compute variance
k <- function(d2, lambda){
  return(exp(- d2/(2*lambda)))
}

K <- k(dist_mat, lambda)
t_t_ind <- matrix(0, n_obs, n_obs)
# t_t_ind[W==1, W==1] <- 1
nt_nt_ind <- matrix(0, n_obs, n_obs)
# nt_nt_ind[W==0, W==0] <- 1
t_nt_ind <- matrix(0, n_obs, n_obs)
# t_nt_ind[W==1, W==0] <- 1
nt_t_ind <- matrix(0, n_obs, n_obs)
# nt_t_ind[W==0, W==1] <- 1

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

K_obs_obs <- K * (t_t_ind * beta01^2 + nt_nt_ind * beta00^2 + t_nt_ind * rho + nt_t_ind * rho)
K_cf_cf <- K * (t_t_ind * beta00^2 + nt_nt_ind * beta01^2 + t_nt_ind * rho + nt_t_ind * rho)
K_obs_cf <- K * (t_t_ind * rho + nt_nt_ind * rho + t_nt_ind * beta01^2 + nt_t_ind *beta00^2)

# K_cf_cf <- k(dist_mat, lambda) * (beta01^2)
sigma_cf <- diag(Y_var * sigmas[2 - W])
# K_obs_cf <- k(dist_mat, lambda) * rho
# K_obs_obs <- k(dist_mat, lambda) * (beta00^2)
sigma_obs <- diag(Y_var * sigmas[W + 1])
# sigma_obs <- diag(rep(1, n_obs))

cov_cf <- K_cf_cf - t(K_obs_cf) %*% solve(K_obs_obs + sigma_obs, K_obs_cf)
Var_cf <- diag(cov_cf)
cov_obs <- K_obs_obs - t(K_obs_obs) %*% solve(K_obs_obs + sigma_obs, K_obs_obs)

print(sum(Var_cf))
##make predictions for mean to check if they're the same
mean_cf <- t(K_obs_cf) %*% solve(K_obs_obs + sigma_obs, Y)
mean_obs <- t(K_obs_obs) %*% solve(K_obs_obs + sigma_obs, Y)

####now compute with RFF####
N_rff <- 200
#weights have sd 1 to start with
# w <- matrix(qnorm(halton(N_covs*N_rff)), nrow=N_covs)
w <- matrix(rnorm(N_covs*N_rff), nrow=N_covs)
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
sigma_obs_inv <- diag(1 / (Y_var * sigmas[W + 1]))
# sigma_obs_inv <- diag(rep(1, n_obs))

##for now assume rho is smaller than beta01^2 and beta00^2
lambda_sq <- sqrt(lambda)
treat_coeff <- sqrt(beta01^2 - rho)
not_treat_coeff <- sqrt(beta00^2 - rho)
w_use <- w / lambda_sq
wX <- t(X %*% w_use)
X_rff_base1 <- rbind(cos(wX), sin(wX)) / sqrt(N_rff)
X_rff_base <- rbind(X_rff_base1, X_rff_base1, X_rff_base1)

# X_rff <- X_rff_base * any_index_mat * rho + X_rff_base * treat_index_mat * treat_coeff + 
#   X_rff_base * not_treat_index_mat * not_treat_coeff


# w_use <- torch$div(w, lambda_sq)
# wX <- t_torch(torch$matmul(X, w_use))
# X_rff_base <- torch$div(torch$cat(c(torch$sin(wX), torch$cos(wX)), as.integer(0)), sqrt(N_rff))
# X_rff_base <- torch$cat(c(X_rff_base, X_rff_base, X_rff_base), as.integer(0))

X_rff <- X_rff_base * (any_index_mat * sqrt(rho) + treat_index_mat * treat_coeff + not_treat_index_mat * not_treat_coeff)
X_rff_cf <- X_rff_base * (any_index_mat * sqrt(rho) + treat_index_mat_cf * treat_coeff + not_treat_index_mat_cf * not_treat_coeff)

K_obs_approx <- t(X_rff) %*% X_rff
K_cf_approx <- t(X_rff_cf) %*% X_rff_cf
K_obs_cf_approx <- t(X_rff) %*% X_rff_cf
plot(K_obs_obs[1:10, ], K_obs_approx[1:10, ], main="K_obs")
abline(0, 1)
plot(K_cf_cf[1:10, ], K_cf_approx[1:10, ], main="K_cf")
abline(0, 1)
plot(K_obs_cf[1:10, ], K_obs_cf_approx[1:10, ], main="K_obs_cf")
abline(0, 1)
# print(cor(as.vector(K_obs_obs[1:10, ]), as.vector(K_obs_approx[1:10, ])))
K_approx <- t(X_rff_base1) %*% X_rff_base1
# plot(K_approx[1:10, ], K[1:10, ])
# abline(0, 1)
# print(median(K_approx[1:10, ] /  K[1:10, ]))

# X_rff <- torch$mul(torch$mul(X_rff_base, any_index_mat), rho) %+_t% 
#   torch$mul(torch$mul(X_rff_base, treat_index_mat), treat_coeff) %+_t%
#   torch$mul(torch$mul(X_rff_base, not_treat_index_mat), not_treat_coeff)
# 
# X_rff_cf <- torch$mul(torch$mul(X_rff_base, any_index_mat), rho) %+_t% 
#   torch$mul(torch$mul(X_rff_base, treat_index_mat), not_treat_coeff) %+_t%
#   torch$mul(torch$mul(X_rff_base, not_treat_index_mat), treat_coeff)

A <- X_rff %*% sigma_obs_inv %*% t(X_rff) + diag(N_rff * 2 * N_covs)
cov_cf_rff <- t(X_rff_cf) %*% solve(A, X_rff_cf)
cov_obs_rff <- t(X_rff) %*% solve(A, X_rff)
# cov_cf_rff <- t(X_rff_cf) %*% A_inv %*% X_rff_cf
Var_cf_rff <- diag(cov_cf_rff)
print(sum(Var_cf_rff))


w_bar_1 <- solve(A, X_rff) %*% sigma_obs_inv %*% Y
# w_bar_2 <- solve(A, t(Y %*% sigma_obs_inv %*% t(X_rff)))
par(mfrow=c(1,2))
plot(cov_obs[1:10, ], cov_obs_rff[1:10, ])
abline(0, 1)
plot(cov_cf[1:10, ], cov_cf_rff[1:10, ])
abline(0, 1)

mean_cf_rff <- t(X_rff_cf) %*% w_bar_1
mean_obs_rff <- t(X_rff) %*% w_bar_1

plot(mean_cf, mean_cf_rff)
abline(0, 1)
plot(mean_obs, mean_obs_rff)
abline(0, 1)
par(mfrow=c(1,1))
# print(sum(Var_cf_rff < 0))
# print(sum(Var_cf < 0))

par(mfrow=c(1, 2))
# diag(cov_cf) <- 0
# diag(cov_cf_rff) <- 0
plot(cov_cf[1:10, ], cov_cf_rff[1:10, ])
abline(0, 1)
plot(Var_cf, Var_cf_rff)
abline(0, 1)
par(mfrow=c(1, 1))
