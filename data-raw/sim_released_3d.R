# Load
load("~/github/mmmTMBarchive2/data/sim_released_3d_yrs_03a_all.rda")
# Assign
sim_released_3d <- sim_released_3d_yrs_03a_all
# Write to data/
usethis::use_data(sim_released_3d, overwrite = TRUE)
