# Three areas with yearly movement
# Data
sim_x <- create_release_matrix()
sim_z <- mmmIndex(3)
sim_f <- as.matrix(0.04)
sim_m <- 0.1
sim_h <- 0.02
sim_u <- 0.1
# Movement rates
sim_r <- matrix(c(0.9, 0.1, 0.0,
                  0.1, 0.8, 0.1,
                  0.0, 0.3, 0.7), byrow = TRUE, nrow = 3)
dim(sim_r) <- c(3, 3, 1, 1)
# Movement parameters
sim_p <- create_movement_parameters(sim_r, sim_z)
# Arguments
data <- list(x = sim_x, z = sim_z, f = sim_f, m = sim_m, h = sim_h, u = sim_u)
parameters <- list(p = sim_p)
# Simulate
sim <- mmmSim(data, parameters)
# Augment data
data$y <- sim$simulation$y
# Simulation data
sim_data <- list(
  x = data$x,
  y = data$y,
  z = data$z,
  f = data$f,
  m = data$m,
  h = data$h,
  u = data$u
)
# Write to data/sim_data.rda
# usethis::use_data(sim_data, overwrite = TRUE)
