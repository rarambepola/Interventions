setwd(paste0(Sys.getenv("HOME"), "/Interventions"))
library(ggplot2)
library(raster)
library(INLA)
library(TMB)
library(Rcpp)


# load("Mock_data/Mdg/Fits/no_selection_bias_fit_n_iter_7.RData")

# load("Mock_data/Mdg/elevation_selection_bias_nonlinear.RData")

load("Mock_data/Mdg/Fits/no_selection_bias_fit_n_iter_100.RData")
# load("Mock_data/Mdg/Fits/elevation_selection_bias_fit_n_iter_100.RData")
load("Mock_data/Mdg/no_selection_bias.RData")
# load("Mock_data/Mdg/elevation_selection_bias.RData")

# load("Mock_data/Mdg/Fits/elevation_selection_bias_nonlinear_fit_n_iter_30.RData")
# load("Mock_data/Mdg/elevation_selection_bias_nonlinear.RData")


# 
# load("Mock_data/Mdg/Fits/elevation_selection_bias_fit_n_iter_50.RData")
# load("Mock_data/Mdg/elevation_selection_bias.RData")

N <- cluster_size
N_pos <- cluster_pos
p_obs <- N_pos / N


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


plot(Y_treatment_effect, Y_pred_treatment_effect)
print(cor(Y_treatment_effect, Y_pred_treatment_effect))

plot(Y_baseline, Y_pred_not_treated)
abline(0, 1)

plot(Y_treated, Y_pred_treated)
abline(0, 1)

plot(Y_treated - Y_baseline, Y_pred_treated - Y_pred_not_treated)
# plot(Y_pred_not_treated, Y_pred_treated)
# abline(0, 1)
# 
# plot(Y_pred_not_treated, environmental_prevalence)


#do inverse with normal logit
p_pred_treated <- 1 / (1 + exp(-Y_pred_treated))
p_pred_not_treated <- 1 / (1 + exp(-Y_pred_not_treated))

plot(p_pred_treated - p_pred_not_treated,
     with_treatment_prevalence - environmental_prevalence)
print(cor(p_pred_treated - p_pred_not_treated,
    with_treatment_prevalence - environmental_prevalence))


#do inverse with empirical logit
inv_emp_logit <- function(y, n){
  q <- 1 / (2*n)
  # return((exp(y) + q*exp(y) - q) / (1 + exp(y)))
  return(1 / (1+exp(-y)))
}

p_pred_treated <- inv_emp_logit(Y_pred_treated, N)
p_pred_not_treated <- inv_emp_logit(Y_pred_not_treated, N)

plot(p_pred_treated - p_pred_not_treated,
     with_treatment_prevalence - environmental_prevalence)
print(cor(p_pred_treated - p_pred_not_treated,
          with_treatment_prevalence - environmental_prevalence))


##make map
rain <- stack("Mock_data/Mdg/rain_stack.tif")[[1]]
temp <- stack("Mock_data/Mdg/temp_stack.tif")[[1]]
elevation <- raster("Mock_data/Mdg/elevation.tif")
mainland <- raster("../Map_madagascar/mainland_raster.tif")

rain <- rain * mainland
temp <- temp * mainland
elevation <- elevation * mainland
W <- covs[, 4]

X <- cbind(values(rain),
            values(temp),
            values(elevation))
which_na <- is.na(rowSums(X))
X <- X[!which_na, ]

X_obs <- covs[, 1:3]
n_obs <- dim(X_obs)[1]
n_pixel <- dim(X)[1]
sourceCpp('Mock_data/Mdg/make_dist_mat.cpp')
#do calculations n_each at a time
n_each <- 10000

N_use <- 1

Y_treated <- rep(NA, n_pixel)
Y_not_treated <- rep(NA, n_pixel)

n_iter <- floor(n_pixel / n_each)
for(i in 1:n_iter){
  print(paste0(i, " out of ", n_iter))
  if(i == n_iter){
    index_i <- (1 + (i-1) * n_each):n_pixel
  }else{
    index_i <- 1:n_each + (i-1) * n_each
  }
  dist_mat <- make_dist_mat(X[index_i, ], X_obs)
  K <- exp(-(dist_mat / lambda))
  
  K_not_treated <- K
  K_not_treated[, W==0] <- K_not_treated[, W==0] * beta00^2
  K_not_treated[, W==1] <- K_not_treated[, W==1] * rho
  
  K_treated <- K
  K_treated[, W==0] <- K_treated[, W==0] * rho
  K_treated[, W==1] <- K_treated[, W==1] * beta01^2
  
  Y_treated[index_i] <- inv_emp_logit(as.vector(K_treated %*% proj_vec), N_use)
  Y_not_treated[index_i] <- inv_emp_logit(as.vector(K_not_treated %*% proj_vec), N_use)
}

treated_r <- rain
not_treated_r <- rain
values(treated_r) <- NA
values(treated_r)[!which_na] <- Y_treated
values(not_treated_r) <- NA
values(not_treated_r)[!which_na] <- Y_not_treated


plot(treated_r, main="treated")
plot(not_treated_r, main="not treated")
par(mfrow=c(1, 2))
# treat_effect_true <- raster("Mock_data/Mdg/treat_effect_nonlinear.tif")
treat_effect_true <- raster("Mock_data/Mdg/treat_effect_linear.tif")
plot(not_treated_r - treated_r, main="treatment_effect pred")
plot(treat_effect_true, main="treatment effect true")
par(mfrow=c(1, 1))
# writeRaster(not_treated_r - treated_r, filename = "Mock_data/Mdg/Fits/treat_pred_effect_nonlinear.tif")
# writeRaster(not_treated_r - treated_r, filename = "Mock_data/Mdg/Fits/treat_pred_effect_linear_elevation_bias.tif")
# writeRaster(not_treated_r - treated_r, filename = "Mock_data/Mdg/Fits/treat_pred_effect_linear_no_bias.tif",
#             overwrite=TRUE)
writeRaster(not_treated_r - treated_r, filename = "Mock_data/Mdg/Fits/treat_pred_effect_linear_elevation_bias.tif",
            overwrite=TRUE)
