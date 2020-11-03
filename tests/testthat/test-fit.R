context("test-fit.R")

# test_that("mmmTMB() returns a list of class 'mmmTMB'", {
# 	# Fit
# 	fit_list <- mmmTMB(
# 	  released_3d = sim_released_3d,
# 	  recovered_5d = sim_recovered_5d,
# 	  capture_rate_2d = sim_capture_rate_2d,
# 	  report_ratio_2d = sim_report_ratio_2d,
# 	  tag_loss_rate = 0.02,
# 	  imm_loss_ratio = 0.1,
# 	  template_2d = sim_template_2d,
# 	  openmp_cores = floor(parallel::detectCores() / 2))
# 	# Tests
# 	testthat::expect_type(fit_list, "list")
# 	testthat::expect_s3_class(fit_list, "mmmTMB")
# })
#
# test_that("mmmTMB() fits a time varying model", {
#   # Fit
#   fit_list <- mmmTMB(
#     released_3d = sim_released_3d,
#     recovered_5d = sim_recovered_5d,
#     capture_rate_2d = sim_capture_rate_2d,
#     report_ratio_2d = sim_report_ratio_2d,
#     tag_loss_rate = 0.02,
#     imm_loss_ratio = 0.1,
#     template_2d = sim_template_2d,
#     time_process = 1,
#     time_pattern = 1,
#     pattern_size = 4,
#     newton_steps = 4, # Optimizer
#     nlminb_loops = 4,
#     openmp_cores = floor(parallel::detectCores() / 2))
#   # Tests
#   testthat::expect_type(fit_list, "list")
#   testthat::expect_s3_class(fit_list, "mmmTMB")
# })
#
# test_that("mmmTMB() computes SEs for size class release and recovery data", {
#   # Fit
#   fit_list <- mmmTMB(
#     released_3d = sim_released_3d_size,
#     recovered_5d = sim_recovered_5d_size,
#     capture_rate_2d = sim_capture_rate_2d,
#     report_ratio_2d = sim_report_ratio_2d,
#     tag_loss_rate = 0.02,
#     imm_loss_ratio = 0.1,
#     template_2d = sim_template_2d,
#     nlminb_loops = 10,
#     newton_steps = 1,
#     openmp_cores = floor(parallel::detectCores() / 2))
#   # Tests
#   testthat::expect_true(sum(
#     which(is.na(fit_list$results$movement_probability_results))) == 0)
# })
