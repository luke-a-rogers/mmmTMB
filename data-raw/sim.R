# Simulate data

# Times: 30
# Areas: 03
# Groups: 01

# Data
sim_mT <- create_release_matrix()
sim_mI <- mmmIndex(3)
sim_mF <- as.matrix(0.04)
sim_sM <- 0.1
# Movement rates
sim_aK <- matrix(c(0.9, 0.1, 0.0,
                   0.1, 0.8, 0.1,
                   0.0, 0.3, 0.7), byrow = TRUE, nrow = 3)
dim(sim_aK) <- c(3, 3, 1, 1)
# Arguments
d <- list(mT = sim_mT, mI = sim_mI, aK = sim_aK, mF = sim_mF, sM = sim_sM)
# Simulate
sim_01 <- mmmSim(d)
sim <- sim_01$simulation

# Write to data/
usethis::use_data(sim, overwrite = TRUE)






# Write to data/
usethis::use_data(sim_mT, overwrite = TRUE)








# Simulate mR

# Times: 30
# Areas: 03
# Groups: 01

sim_mF <- as.matrix(0.04)
sim_sM <- 0.1
d <- list(mT = sim_mT, mI = sim_mI, aK = sim_aK, mF = sim_mF, sM = sim_sM)
s1 <- mmmSim(d)

sim_mR <- s1$data$mR

# Write to data/
usethis::use_data(sim_mR, overwrite = TRUE)
