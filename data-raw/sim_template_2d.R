# Load
load("~/github/mmmTMBarchive2/data/sim_template_2d_yrs_03a.rda")
# Assign
sim_template_2d <- sim_template_2d_yrs_03a
# Write to data/
usethis::use_data(sim_template_2d, overwrite = TRUE)
