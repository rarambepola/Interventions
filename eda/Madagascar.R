setwd(paste0(Sys.getenv("HOME"), "/Interventions"))

library(ggplot2)

#
itn_dhs <- read.csv("DHS/Madagascar/Svy_505_ITN_HH_Res.csv")
itn_dhs_u5 <- itn_dhs[itn_dhs$n_pop_u5 > 0, ]


#get cluster coordinates
n_cluster <- max(itn_dhs$clusterid)
cluster_size <- sapply(1:n_cluster, function(i) sum(itn_dhs_u5$n_pop_u5[itn_dhs_u5$clusterid == i]))
cluster_coords <- cbind(itn_dhs$longitude,
                        itn_dhs$latitude)[match(1:n_cluster, itn_dhs$clusterid), ]
save(list=c("cluster_coords", "cluster_size"),
     file = "Mock_data/Mdg/cluster_coords.RData")
plot(cluster_coords, asp=1)

#calculate cluster itn coverage
itn_coverage_u5 <- sapply(1:n_cluster, function(i) sum(itn_dhs_u5$n_u5_under_itn[itn_dhs_u5$clusterid == i]) / sum(itn_dhs_u5$n_pop_u5[itn_dhs_u5$clusterid == i]))
itn_u5_df <- data.frame(x=cluster_coords[, 1],
                        y=cluster_coords[, 2],
                        coverage=itn_coverage_u5,
                        age="u5")

itn_coverage_all <- sapply(1:n_cluster, function(i) sum(itn_dhs$n_slept_under_itn[itn_dhs$clusterid == i]) / sum(itn_dhs$n_defacto_pop[itn_dhs$clusterid == i]))
itn_all_df <- data.frame(x=cluster_coords[, 1],
                        y=cluster_coords[, 2],
                        coverage=itn_coverage_all,
                        age="all")

p <- ggplot(rbind(itn_u5_df, itn_all_df)) + geom_point(aes(x=x, y=y, col=coverage)) + scale_color_viridis_c() + 
  coord_fixed() + theme_minimal() + facet_wrap(.~age)
print(p)

##look at cube surfaces


#look at useful categories
bin_upper_vals <- c(0.1, 0.2, 0.6, 0.7, 0.8, 0.9)
bin_lower_vals <- c(0, bin_upper_vals[-length(bin_upper_vals)])
bin_names <- paste0(bin_lower_vals, " - ", bin_upper_vals)
bin_names[length(bin_names) + 1] <- paste0("> ", bin_upper_vals[length(bin_upper_vals)])
n_bins <- length(bin_names)
pal <- c(rev(RColorBrewer::brewer.pal(n_bins, "RdYlBu")))
names(pal) <- bin_names


itn_all_df$"coverage_bin" <- sapply(itn_all_df$coverage, function(val) bin_names[sum(val > bin_upper_vals) + 1])

p <- ggplot(itn_all_df) + geom_point(aes(x=x, y=y, col=coverage_bin)) + scale_fill_manual("Coverage", values = pal) + 
  coord_fixed() + theme_minimal() 
print(p)
