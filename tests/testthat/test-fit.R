context("test-fit.R")

test_that("mmmTMB() returns a list of class 'mmmTMB'", {
  # Arguments
  data_list <- list(
    released_3d = sim_released_3d,
    recovered_5d = sim_recovered_5d,
    capture_rate_2d =  sim_capture_rate_2d,
    report_ratio_2d = sim_report_ratio_2d,
    template_2d = sim_template_2d,
    tag_loss_rate = 0.02,
    imm_loss_ratio = 0.1
  )
  structure_list <- list(
    recapture_delay = 1,
    error_family = 1,
    result_units = 1
  )
  optimizer_list <- list(
    newton_steps = 1,
    nlminb_loops = 5,
    nlminb_control = mmmTMBcontrol(),
    openmp_cores = floor(parallel::detectCores() / 2)
  )

	# Fit
	fit_list <- mmmTMB(data_list = data_list,
	                   structure_list = structure_list,
	                   optimizer_list)
	testthat::expect_type(fit_list, "list")
	testthat::expect_s3_class(fit_list, "mmmTMB")
})


