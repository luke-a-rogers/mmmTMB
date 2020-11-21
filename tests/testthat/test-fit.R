context("test-fit.R")

test_that("mmmFit() returns a list of class 'mmmFit'", {
  # Fit
  fit <- mmmFit(
    data = sysdata$data,
    parameters = NULL,
    settings = mmmSet(
      error_family = 0,
      nlminb_loops = 5,
      newton_iters = 5,
      openmp_cores = 6)
  )
  # Test
  testthat::expect_type(fit, "list")
  testthat::expect_s3_class(fit, "mmmFit")
})
