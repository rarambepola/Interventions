setwd(paste0(Sys.getenv("HOME"), "/Interventions"))
library(ggplot2)
library(raster)
library(INLA)
library(TMB)

# load("Mock_data/Mdg/no_selection_bias.RData")
# load("Mock_data/Mdg/elevation_selection_bias.RData")
load("Mock_data/Mdg/elevation_selection_bias_nonlinear.RData")
mesh <- inla.mesh.2d(loc = cluster_coords, cutoff = 0.25, max.edge = c(1, 3))
plot(mesh, asp=1)
mesh_coords <- mesh$loc[, 1:2]
n_mesh <- mesh$n

alpha <- 2
nu <- alpha - 1
spde <- (inla.spde2.matern(mesh=mesh, alpha=alpha)$param.inla)[c("M0","M1","M2")]
A <- inla.spde.make.A(mesh=mesh, loc=as.matrix(cluster_coords))


#compile model
model_folder <- "Mock_data/Mdg/"
model_name <- "simple"
model_path <- paste0(model_folder, model_name)
tryCatch(dyn.unload(dynlib(model_path)),
         error = function(e) print(e))
compile(paste0(model_path, ".cpp"))
dyn.load(dynlib(model_path))

N_covs <- dim(covs)[2]
model_silent <- TRUE
m <- MakeADFun(
  data = list(X=covs,
              Y=cluster_pos,
              pops=cluster_size,
              A=A,
              spde=spde,
              prior_rho_min = 1,
              prior_rho_prob = 0.1,
              prior_sigma_max = 2,
              prior_sigma_prob = 0.1,
              nu=nu
  ),
  parameters = list(beta_0=runif(1, -1, 1), 
                    beta=rep(0, N_covs),
                    S=rnorm(n_mesh),
                    log_rho=1.0,
                    log_sigma=0.0
  ),
  random="S",
  DLL = model_name,
  silent=model_silent
)


ptm <- proc.time()
fit <- nlminb(m$par, m$fn, m$gr, control=list(iter.max=300,eval.max=300))
ptm2 <- proc.time()
print("Time to fit")
print(ptm2 - ptm)

print(fit$par[5])
