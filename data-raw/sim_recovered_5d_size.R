# Define
nt <- dim(sim_recovered_5d)[1]
na <- dim(sim_recovered_5d)[3]
ng <- 3
# Assign
sim_recovered_5d_size <- array(0, dim = c(nt, nt, na, na, ng))
sim_recovered_5d_size[, , , , 1] <- sim_recovered_5d
sim_recovered_5d_size[, , , , 2] <- sim_recovered_5d
sim_recovered_5d_size[, , , , 3] <- sim_recovered_5d
# Write to data/
usethis::use_data(sim_recovered_5d_size, overwrite = TRUE)
