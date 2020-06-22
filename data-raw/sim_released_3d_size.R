# Define
nt <- dim(sim_released_3d)[1]
na <- dim(sim_released_3d)[2]
ng <- 3
# Assign
sim_released_3d_size <- array(0, dim = c(nt, na, ng))
sim_released_3d_size[, , 1] <- sim_released_3d
sim_released_3d_size[, , 2] <- sim_released_3d
sim_released_3d_size[, , 3] <- sim_released_3d
# Write to data/
usethis::use_data(sim_released_3d_size, overwrite = TRUE)
