setwd(paste0(Sys.getenv("HOME"), "/Interventions"))
library(ggplot2)
library(raster)

rain_stack <- stack("../Map_madagascar/Model_2020/rain_stack.tif")
temp_stack <- stack("../Map_madagascar/Model_2020/lst_stack.tif")


#normalise and save 12 months
rain_mean <- mean(values(rain_stack$rain_stack.3), na.rm=T)
rain_raster_list <- list()
for(i in 1:12){
  
}

#choose parameters
beta0 <- runif(1, -6, -4)
beta <- rnorm()
