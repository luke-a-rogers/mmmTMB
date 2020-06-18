# Load
load("~/github/mmmTMBarchive2/data/sim_capture_rate_2d_yrs_03a.rda")
# Assign
sim_capture_rate_2d <- sim_capture_rate_2d_yrs_03a
# Write to data/
usethis::use_data(sim_capture_rate_2d, overwrite = TRUE)
