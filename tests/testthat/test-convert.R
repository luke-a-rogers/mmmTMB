context("test-convert.R")

test_that("logit() and invlogit() are inverses of eachother", {
  # Inputs
  x1 <- 0
  x2 <- 0.5
  x3 <- 1
  y1 <- -Inf
  y2 <- 0
  y3 <- Inf
  # Outputs
  p1 <- invlogit(logit(x1))
  p2 <- invlogit(logit(x2))
  p3 <- invlogit(logit(x3))
  r1 <- logit(invlogit(y1))
  r2 <- logit(invlogit(y2))
  r3 <- logit(invlogit(y3))
  # Test
  testthat::expect_equal(p1, x1)
  testthat::expect_equal(p2, x2)
  testthat::expect_equal(p3, x3)
  testthat::expect_equal(r1, y1)
  testthat::expect_equal(r2, y2)
  testthat::expect_equal(r3, y3)
})
