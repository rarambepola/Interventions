setwd(paste0(Sys.getenv("HOME"), "/Interventions"))
library(ggplot2)
library(raster)

save_plot <- TRUE


#palette
# bin_upper_vals <- c(0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3)
bin_upper_vals <- c(0.05, 0.075, 0.1, 0.15)
bin_lower_vals <- c(0, bin_upper_vals[-length(bin_upper_vals)])
bin_names <- paste0(bin_lower_vals, " - ", bin_upper_vals)
bin_names[length(bin_names) + 1] <- paste0("> ", bin_upper_vals[length(bin_upper_vals)])
n_bins <- length(bin_names)
pal <- c(rev(RColorBrewer::brewer.pal(n_bins, "RdYlBu")))
names(pal) <- bin_names


##MDG
true_r <- raster("Mock_data/Mdg/treat_effect_linear.tif")
# true_r <- raster("Mock_data/Mdg/treat_effect_nonlinear.tif")
pred_r <- raster("Mock_data/Mdg/Fits/treat_pred_effect_linear_no_bias.tif")
# pred_r <- raster("Mock_data/Mdg/Fits/treat_pred_effect_linear_elevation_bias.tif")
# pred_r <- raster("Mock_data/Mdg/Fits/treat_pred_effect_nonlinear.tif")

which_na <- which(is.na(values(true_r)))
coords_r <- coordinates(true_r)[-which_na, ]
n_pixel <- length(values(true_r)) - length(which_na)

#make data frame
plot_df <- data.frame(x=coords_r[, 1],
                      y=coords_r[, 2],
                      z=c(values(true_r)[-which_na], values(pred_r)[-which_na]),
                      type=rep(c("True", "Predicted"), each=n_pixel)
)

plot_df$"vals" <- sapply(plot_df$z, function(val) bin_names[sum(val > bin_upper_vals) + 1])

# p <- ggplot(plot_df, aes(x=x, y=y, fill=vals)) + geom_tile() + 
p <- ggplot(plot_df, aes(x=x, y=y, fill=z)) + geom_tile() + 
  theme_minimal() + coord_fixed() + theme(panel.grid.major = NULL, panel.grid.minor = NULL) + 
  # scale_fill_manual("Treatment effect", values = pal, na.translate=FALSE, position="bottom") + 
  scale_fill_viridis_c("Treatment effect") + 
  facet_wrap(~type, nrow=1) + 
  xlab(NULL) + ylab(NULL)

if(save_plot){
  # ggsave(p, filename = "Mock_data/Mdg/Fits/elevation_bias_nonlinear.png", width=8, units="in")
  ggsave(p, filename = "Mock_data/Mdg/Fits/no_bias_linear_cont.png", width=8, units="in")
  # ggsave(p, filename = "Mock_data/Mdg/Fits/elevation_bias_linear_cont.png", width=8, units="in")
  
}

