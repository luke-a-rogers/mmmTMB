context("test-fit.R")

test_that("mmmTMB() returns a list of class 'mmmTMB'", {
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
	# Tests
	testthat::expect_type(fit_list, "list")
	testthat::expect_s3_class(fit_list, "mmmTMB")
})


