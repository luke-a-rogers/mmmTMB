# Fit
fit_list <- mmmTMB(
  released_3d = sim_released_3d,
  recovered_5d = sim_recovered_5d,
  capture_rate_2d = sim_capture_rate_2d,
  report_ratio_2d = sim_report_ratio_2d,
  tag_loss_rate = 0.02,
  imm_loss_ratio = 0.1,
  template_2d = sim_template_2d,
  openmp_cores = floor(parallel::detectCores() / 2))
# Extract
smp <- subset_by_name(fit_list$model$par, "movement_parameters_3d")
smp_3d <- fit_list$adfun$report()$movement_parameters_3d
smp_4d <- fit_list$adfun$report()$movement_probability_4d
# Create
sysdata <- list(sim_movement_parameters = smp,
                sim_movement_parameters_3d = smp_3d,
                sim_movement_probability_4d = smp_4d)
# Write to R/sysdata.rda
usethis::use_data(sysdata, internal = TRUE, overwrite = TRUE)
