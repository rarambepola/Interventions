setwd(paste0(Sys.getenv("HOME"), "/Interventions"))
library(ggplot2)
library(raster)

set.seed(1)

# rain_stack <- stack("../Map_madagascar/Model_2020/rain_stack.tif")
# temp_stack <- stack("../Map_madagascar/Model_2020/lst_stack.tif")
# 
# 
# #normalise and save 12 months
# rain_mean <- mean(values(rain_stack$rain_stack.3), na.rm=T)
# rain_sd <- sd(values(rain_stack$rain_stack.3), na.rm=T)
# temp_mean <- mean(values(temp_stack$lst_stack.3), na.rm=T)
# temp_sd <- sd(values(temp_stack$lst_stack.3), na.rm=T)
# 
# rain_raster_list <- list()
# temp_raster_list <- list()
# 
# for(i in 1:12){
#   rain_raster_list[[i]] <- (rain_stack[[i]] - rain_mean) / rain_sd
#   temp_raster_list[[i]] <- (temp_stack[[i]] - temp_mean) / temp_sd
# }
# static_stack <- stack("../Map_madagascar/Model_2020/raster_covs_static.tif")
# elevation_r <- (static_stack$raster_covs_static.3 - mean(values(static_stack$raster_covs_static.3), na.rm=T)) / sd(values(static_stack$raster_covs_static.3), na.rm=T)
# writeRaster(elevation_r, filename = "Mock_data/Mdg/elevation.tif", overwrite=TRUE)
# 
# 
# writeRaster(stack(rain_raster_list), filename = "Mock_data/Mdg/rain_stack.tif", overwrite=TRUE)
# writeRaster(stack(temp_raster_list), filename = "Mock_data/Mdg/temp_stack.tif", overwrite=TRUE)

rain <- stack("Mock_data/Mdg/rain_stack.tif")[[1]]
temp <- stack("Mock_data/Mdg/temp_stack.tif")[[1]]
elevation <- raster("Mock_data/Mdg/elevation.tif")
mainland <- raster("../Map_madagascar/mainland_raster.tif")

# rain_na <- which(is.na(values(rain)) & !is.na(values(mainland)))
# temp_na <- which(is.na(values(temp)) & !is.na(values(mainland)))
# elevation_na <- which(is.na(values(elevation)) & !is.na(values(mainland)))
rain_na <- which(is.na(values(rain)))
temp_na <- which(is.na(values(temp)))
elevation_na <- which(is.na(values(elevation)))

n_blur <- 5
w <- matrix(1/(n_blur * n_blur), n_blur, n_blur)
values(temp)[temp_na] <- values(focal(temp, w, sum, na.rm=T))[temp_na]
values(rain)[rain_na] <- values(focal(rain, w, sum, na.rm=T))[rain_na]
values(elevation)[elevation_na] <- values(focal(elevation, w, sum, na.rm=T))[elevation_na]


#choose parameters
beta0 <- runif(1, -2, 0)
beta <- rnorm(2, sd=0.5)
beta_elev <- -0.5

i <- 1
logit_p <- beta0 + rain * beta[1] + temp * beta[2] + elevation * beta_elev
p_raster <- 1 / (1 + exp(-logit_p))
plot(logit_p)

n_obs <- dim(cluster_coords)[1]



#use same cluster locations as 2016 DHS
load("Mock_data/Mdg/cluster_coords.RData")


#make some data with no selection bias for treatment
#binary treatment with fixed effect for now
beta_treat <- -0.8
environmental_risk <- extract(logit_p, cluster_coords)
#remove areas with NAs
which_na <- which(is.na(environmental_risk))
environmental_risk <- environmental_risk[-which_na]
cluster_coords <- cluster_coords[-which_na, ]
cluster_size <- cluster_size[-which_na]
n_locs <- dim(cluster_coords)[1]
treatment_vec <- sample(c(0, 1), n_locs, replace = TRUE)
true_risk <- environmental_risk + beta_treat * treatment_vec
with_treatment_risk <- environmental_risk + beta_treat

logit_p_w_treatment <- logit_p + beta_treat
p_treat_raster <- 1 / (1 + exp(-logit_p_w_treatment))
  
treat_effect <- p_raster - p_treat_raster
writeRaster(treat_effect, file="Mock_data/Mdg/treat_effect_linear.tif",
            overwrite=TRUE)

plot_df <- data.frame(x=cluster_coords[, 1],
                      y=cluster_coords[, 2],
                      logit_p=true_risk,
                      p=1 / (1 + exp(-true_risk)),
                      treatment=treatment_vec)
p <- ggplot(plot_df) + geom_point(aes(x=x, y=y, col=p)) + theme_minimal() + 
  coord_fixed() + scale_color_viridis_c() + facet_wrap(.~treatment) + ggtitle("True risk")
print(p)

cluster_pos <- rbinom(n_locs, size = cluster_size, prob = 1 / (1 + exp(-true_risk)))

plot_df$pfpr <- cluster_pos / cluster_size
p <- ggplot(plot_df) + geom_point(aes(x=x, y=y, col=pfpr)) + theme_minimal() + 
  coord_fixed() + scale_color_viridis_c() + facet_wrap(.~treatment) + ggtitle("pfpr")
print(p)

#make design matrix
covs <- cbind(extract(rain, cluster_coords),
              extract(temp, cluster_coords),
              extract(elevation, cluster_coords),
              treatment_vec)

environmental_prevalence <- 1 / (1 + exp(-environmental_risk))
true_prevalence <- 1 / (1 + exp(-true_risk))
with_treatment_prevalence <- 1 / (1 + exp(-with_treatment_risk))

save(list=c("environmental_risk", "true_risk", "beta_treat", "cluster_coords", "cluster_size",
            "cluster_pos", "covs", "beta0", "beta", "environmental_prevalence",
            "true_prevalence", "with_treatment_prevalence"),
     file="Mock_data/Mdg/no_selection_bias.RData")

#biased to be lower in areas of high elevation
elevation <- extract(elevation, cluster_coords)
treatment_vec <- rbinom(n = n_locs, size = 1, prob = 1 / (1+exp(elevation)))
true_risk <- environmental_risk + beta_treat * treatment_vec
cluster_pos <- rbinom(n_locs, size = cluster_size, prob = 1 / (1 + exp(-true_risk)))
true_prevalence <- 1 / (1 + exp(-true_risk))
covs <- cbind(covs[, 1:3],
              treatment_vec)

save(list=c("environmental_risk", "true_risk", "beta_treat", "cluster_coords", "cluster_size",
            "cluster_pos", "covs", "beta0", "beta", "environmental_prevalence",
            "true_prevalence", "with_treatment_prevalence"),
     file="Mock_data/Mdg/elevation_selection_bias.RData")


#non linear treatment effect
treatment_frac <- 0.5
treatment_effect <- rep(1, n_obs)
treatment_effect[treatment_vec == 1] <- treatment_frac
# true_risk <- environmental_risk * treatment_effect
environmental_prevalence <- 1 / (1 + exp(-environmental_risk))
true_prevalence <- environmental_prevalence  * treatment_effect
cluster_pos <- rbinom(n_locs, size = cluster_size, prob = true_prevalence)
with_treatment_prevalence <- environmental_prevalence * treatment_frac

treat_effect <- 0.5 * p_raster
writeRaster(treat_effect, file="Mock_data/Mdg/treat_effect_nonlinear.tif",
            overwrite=TRUE)

save(list=c("environmental_risk", "true_risk", "treatment_frac", "cluster_coords", "cluster_size",
            "cluster_pos", "covs", "beta0", "beta", "environmental_prevalence",
            "true_prevalence", "with_treatment_prevalence"),
     file="Mock_data/Mdg/elevation_selection_bias_nonlinear.RData")
