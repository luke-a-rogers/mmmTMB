# Load
load("~/github/mmmTMBarchive2/data/sim_report_ratio_2d_yrs_03a.rda")
# Assign
sim_report_ratio_2d <- sim_report_ratio_2d_yrs_03a
# Write to data/
usethis::use_data(sim_report_ratio_2d, overwrite = TRUE)
