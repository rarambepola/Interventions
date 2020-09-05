setwd(paste0(Sys.getenv("HOME"), "/Interventions/"))
rm(list = ls())
library(ggplot2)


# set.seed(2)

# par(mfrow=c(2,2))

##mock data 
N_obs <- 50
x_max <- 13
x_min <- 0
x <- sort(runif(N_obs, x_min, x_max))
f <- function(x){
  return(0.25*x + sin(x/1.2) + 1)
}

N_plot <- 300
x_plot <- seq(x_min, x_max, length.out=N_plot)

y_baseline <- rnorm(N_obs, mean=0, sd=0.2)
# treat_effect <- -2
treat_ratio <- 0.6
treat_vec <- rbinom(N_obs, 1, 0.5)
x_cutoff <- 7
treat_vec[x > x_cutoff] <- rbinom(sum(x > x_cutoff), 1, 0.1)

# y <- y_baseline + treat_effect * treat_vec + f(x)
y <- y_baseline + f(x) * (treat_vec  * treat_ratio + !treat_vec)


plot_df <- data.frame(x=x, y=y, type=as.factor(treat_vec))
true_df <- data.frame(x=rep(x_plot, 2),
                      y = c(f(x_plot), f(x_plot) * treat_ratio),
                      # y=c(f(x_plot), f(x_plot) + treat_effect),
                      type=as.factor(rep(0:1, each=N_plot)))

plot(plot_df[, 1:2], cex=0.001)
points(plot_df[plot_df$type == 0, 1:2], col="red", pch=19, cex=0.5)
points(plot_df[plot_df$type == 1, 1:2], col="blue", pch=19, cex=0.5)
lines(true_df[true_df$type == 0, 1:2], col="red", lty=3)
lines(true_df[true_df$type == 1, 1:2], col="blue", lty=3)

##choose hyperparameters and fit GP
sigma <- 0.2
lambda <- 0.1
rho_scale <- 0
rho <- rho_scale * sigma^2 
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


pred_df <- data.frame(x=x_pred, y=y_pred, type=treat_pred)

# plot(plot_df[, 1:2], cex=0.001)
# points(plot_df[treat_vec == 0, 1:2], col="red", pch=19, cex=0.5)
# points(plot_df[treat_vec == 1, 1:2], col="blue", pch=19, cex=0.5)
# lines(true_df[true_df$type == 0, 1:2], col="red", lty=3)
# lines(true_df[true_df$type == 1, 1:2], col="blue", lty=3)
# lines(pred_df[pred_df$type == 0, 1:2], col="red", lty=1)
# lines(pred_df[pred_df$type == 1, 1:2], col="blue", lty=1)


##add CIs
treatment_same_mat_pred_pred <- matrix(0, N_pred, N_pred)
treatment_diff_mat_pred_pred <- matrix(0, N_pred, N_pred)

dist_mat_pred_pred <- matrix(NA, N_pred, N_pred)
for(i in 1:N_pred){
  for(j in 1:N_pred){
    dist_mat_pred_pred[i, j] <- sum((x_pred[i] - x_pred[j])^2)
  }
}

for(i in 1:N_pred){
  treatment_same_mat_pred_pred[i, treat_pred[i] == treat_pred] <- 1
  treatment_diff_mat_pred_pred[i, treat_pred[i] != treat_pred] <- 1
}


K_pred_pred <- k(dist_mat_pred_pred, sigma, lambda)


K_pred_pred_use <- K_pred_pred * treatment_same_mat_pred_pred + (K_pred_pred * rho * treatment_diff_mat_pred_pred/ sigma^2)

y_sd_pred <- sqrt(diag(K_pred_pred - K_pred %*% solve(K + sigma_obs * diag(N_obs), t(K_pred))))

mean_df <- data.frame(x=x_pred,
                      y=y_pred,
                      type=as.factor(treat_pred))
CI_df <- data.frame(x=x_pred,
                    ymin=y_pred - 2 * y_sd_pred,
                    ymax=y_pred + 2 * y_sd_pred,
                    type=as.factor(treat_pred))
p <- ggplot() + geom_ribbon(data=CI_df, aes(x=x, ymin=ymin, ymax=ymax, group=type), alpha=0.1) + 
  geom_line(data=mean_df, aes(x=x, y=y, col=type)) + theme_minimal() + 
  geom_point(data=plot_df, aes(x=x, y=y, col=type)) + 
  geom_line(data=true_df, aes(x=x, y=y, group=type), linetype="dashed", alpha=0.9)
# print(p)

effect_df <- data.frame(x=c(true_df$x[true_df$type==1], mean_df$x[mean_df$type==1]),
                        y=c(true_df$y[true_df$type==1] - true_df$y[true_df$type==0],
                            mean_df$y[mean_df$type==1] - mean_df$y[mean_df$type==0]),
                        type=rep(c("true", "pred"), c(N_plot, N_pred/2)))

effect_p <- ggplot(data=effect_df) + geom_line(aes(x=x, y=y, col=type)) + theme_minimal()
# print(effect_p)

# print(ggpubr::ggarrange(p, effect_p))

bias_df <- data.frame(x=plot_df$x,
                      y=abs(y_fit - plot_df$y),
                      type=plot_df$type)

bias_p <- ggplot(bias_df) + geom_point(aes(x=x, y=y)) + 
  geom_smooth(aes(x=x, y=y), se=FALSE, alpha=0.1) + theme_minimal()
print(ggpubr::ggarrange(p, effect_p, bias_p, nrow=1))

