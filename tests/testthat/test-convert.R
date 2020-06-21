context("test-convert.R")

test_that("create_movement_parameters_3d() matches movement_parameters_3d", {
  # Assign
  smp <- sysdata$sim_movement_parameters
  smp_3d <- sysdata$sim_movement_parameters_3d
  # Convert
  nv <- dim(smp_3d)[1]
  np <- dim(smp_3d)[2]
  ng <- dim(smp_3d)[3]
  mp_3d <- create_movement_parameters_3d(pars = smp, nv, np, ng)
  # Tests
  testthat::expect_equal(dim(mp_3d), dim(smp_3d))
  testthat::expect_equal(mp_3d, smp_3d)
})

test_that("create_movement_probability_4d() matches movement_probability_4d", {
  # Assign
  smp_3d <- sysdata$sim_movement_parameters_3d
  smp_4d <- sysdata$sim_movement_probability_4d
  # Convert
  na <- dim(smp_4d)[1]
  nv <- dim(smp_4d)[3]
  ng <- dim(smp_4d)[4]
  mp_4d <- create_movement_probability_4d(na, nv, ng, smp_3d, sim_template_2d)
  # Tests
  testthat::expect_equal(dim(mp_4d), dim(smp_4d))
  testthat::expect_equal(mp_4d, smp_4d)
})
