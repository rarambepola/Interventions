setwd(paste0(Sys.getenv("HOME"), "/Interventions"))
load("Mock_data/Mdg/elevation_selection_bias.RData")

#fix some hyperparameter vals and compute risk
# beta00 <- 0.9
# beta01 <- 1.1
# rho <- 0.3
# lambda <- 1
# sigma_0 <- 1
# sigma_1 <- 0.5
# sigmas <- c(sigma_0, sigma_1)

# beta00 <- 1
# beta01 <- 1
# rho <- 1
# lambda <- 1
# sigma_0 <- 1
# sigma_1 <- 1

beta00 <- params_r[1]
beta01 <- params_r[2]
rho <- params_r[3]
lambda <- params_r[4]
sigma_0 <- params_r[5]
sigma_1 <- params_r[6]

# rho_scale <- 1 / (1 + exp(-params_r[3]))
# rho <- sqrt(beta00 * beta01) * rho_scale

sigmas <- c(sigma_0, sigma_1)



X <- covs[, 1:3]
W <- covs[, 4]
N <- cluster_size
N_pos <- cluster_pos

Y <- log((N_pos + 0.5) / (N - N_pos + 0.5))
Y_var <- (1 / (N_pos + 0.5)) + (1 / (N - N_pos + 0.5))
# Y_var <- 1
# sigmas <- c(1,1)

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

#compute variance
k <- function(d2, lambda){
  return(exp(- d2/(1*lambda)))
}

K <- k(dist_mat, lambda)
K_cf_cf <- k(dist_mat, lambda) * (beta01^2)
sigma_cf <- diag(Y_var * sigmas[2 - W])
K_obs_cf <- k(dist_mat, lambda) * rho
K_obs_obs <- k(dist_mat, lambda) * (beta00^2)
sigma_obs <- diag(Y_var * sigmas[W + 1])

Var_cf <- sum(diag(K_cf_cf - K_obs_cf %*% solve(K_obs_obs + sigma_obs, t(K_obs_cf))))
print(Var_cf)

# #compute LOO-CV diff
# risk_sum <- 0
# for(i in 1:n_obs){
#   K_obs_all_minus_i <- k(dist_mat[-i, i], lambda) * beta00 
#   K_cf_minus_i_minus_i <- k(dist_mat[-i, -i], lambda) * beta01
#   exp_f <- K_obs_all_minus_i %*% solve(K_cf_minus_i_minus_i + sigma_cf[-i, -i], Y[-i])
#   risk_sum <- risk_sum + (exp_f - Y[i]) ^ 2
# }


##LOO-CV
Y_fit_R <- c()
dist_sum <- 0
for(i in 1:n_obs){
  K_12_R <- K[-i, i]
  Y_fit_R[i] <- K_12_R %*% solve(K[-i, -i] + sigma_obs[-i, -i], Y[-i])
  dist_sum <- dist_sum + (Y_fit_R[i] - Y[i])^2
}
dist_sum <- dist_sum / n_obs


