# Load
load("~/github/mmmTMBarchive2/data/sim_recovered_5d_yrs_03a_all.rda")
# Assign
sim_recovered_5d <- sim_recovered_5d_yrs_03a_all
# Write to data/
usethis::use_data(sim_recovered_5d, overwrite = TRUE)
