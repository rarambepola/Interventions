setwd(paste0(Sys.getenv("HOME"), "/Interventions"))

itn_dhs <- read.csv("DHS/Madagascar/Svy_505_ITN_HH_Res.csv")
# pr_dhs <- read.csv("Real_data/Mdg/Svy_505_PR_CLUST_RDT_Res.csv")
pr_dhs <- read.csv("Real_data/Mdg/Madagascar_2016.csv")
pr_dhs <- pr_dhs[, -c(18:36)]

#get cluster prevalence
pr_examined <- c()
pr_pos <- c()

cluster_nos <- sort(unique(pr_dhs$cluster_number))
n_cl <- length(cluster_nos)

for(i in 1:n_cl){
  cluster_no <- cluster_nos[i]
  pr_cluster <- pr_dhs[pr_dhs$cluster_number == i, ]
  
  pr_examined[i] <- sum(pr_cluster$malaria_rapid_test_result_description %in% c("negative", "positive"),
                        na.rm=T)
  pr_pos[i] <- sum(pr_cluster$malaria_rapid_test_result_description == "positive",
                   na.rm=T)
}

pr_pos <- pr_pos[pr_examined > 0]
cluster_nos <- cluster_nos[pr_examined > 0]
pr_examined <- pr_examined[pr_examined > 0]
pfpr <- pr_pos / pr_examined
