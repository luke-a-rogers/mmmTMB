#' Logit (Log Odds Ratio)
#'
#' @param x Numeric between zero and one(inclusive)
#'
#' @return The log of the odds ratio
#' @export
#'
#' @examples
#'
logit <- function (x) {
  # Check arguments
  checkmate::assert_numeric(x, lower = 0, upper = 1)
  # Return logit(x)
  return(log(x / (1 - x)))
}

#' Inverse Logit
#'
#' @param y Numeric
#'
#' @return The inverse of the log odds ratio
#' @export
#'
#' @examples
#'
invlogit <- function (y) {
  # Check arguments
  checkmate::assert_numeric(y)
  # Compute invlogit(y)
  ans <- numeric(length = length(y))
  ans[which(y == Inf)] <- 1
  ans[which(y != Inf)] <- exp(y) / (1 + exp(y))
  dim(ans) <- dim(y)
  # Return invlogit(y)
  return(ans)
}

#' Create Movement Rate Array
#'
#' @param a Numeric array. Movement parameters. (Not transposed).
#' @param m Integer matrix. Movement indexes. (Not transposed).
#'
#' @return [array()] Movement rates
#' @export
#'
#' @examples
#'
create_movement_rates <- function (a, m) {

  # Check arguments
  checkmate::assert_array(a, mode = "numeric", any.missing = FALSE)
  checkmate::assert_array(a, d = 3)
  checkmate::assert_matrix(m, mode = "integerish", any.missing = FALSE)
  checkmate::assert_numeric(m, lower = 0, upper = 1)
  # Convert to transposes
  ta <- aperm(a, c(2,1,3))
  tm <- t(m)
  # Initialize index limits
  na <- nrow(tm) # Matrix tm has rows = na, and cols = na
  npt <- dim(ta)[2] # Array ta has dim = c(np, npt, ng)
  ng <- dim(ta)[3]
  # Initialize array
  k <- array(0, dim = c(na, na, npt, ng))
  # Compute movement rates from movement parameters
  for (mg in seq_len(ng)) {
    for (cpt in seq_len(npt)) {
      # Set the current parameter indexes to one
      cp_num <- 1L
      cp_den <- 1L
      for (pa in seq_len(na)) {
        # Set the denominator sum to zero for this row
        d_sum <- 0
        for (ca in seq_len(na)) {
          if (tm[ca, pa] != 0) {
            d_sum <- d_sum + exp(ta[cp_den, cpt, mg])
            cp_den <- cp_den + 1
          }
        }
        # Set the probability sum to zero
        k_sum <- 0L
        # Compute the probability for ca != pa
        for (ca in seq_len(na)) {
          if (tm[ca, pa] != 0) {
            k[pa, ca, cpt, mg] <- exp(ta[cp_num, cpt, mg]) / (1L + d_sum)
            k_sum <- k_sum + k[pa, ca, cpt, mg]
            cp_num <- cp_num + 1
          }
        }
        # Compute the probability for pa == ca
        k[pa, pa, cpt, mg] <- 1L - k_sum
      }
    }
  }
  return(k)
}

#' Create Movement Parameter Array
#'
#' @param x Numeric Array. aK.
#' @param m Integer matrix. Movement indexes. dim = c(pa, ca)
#'
#' @return [array()]
#' @export
#'
#' @examples
#'
#' # 2x2
#' x <- array(c(0.9, 0.2, 0.1, 0.8, 0.7, 0.4, 0.3, 0.6), dim = c(2, 2, 2, 1))
#' m <- array(c(0, 1, 1, 0), dim = c(2, 2))
#' a <- create_movement_parameters(x, m)
#' k <- create_movement_rates(a, m)
#' all(k - x < 1e-12)
#'
#' # 3x3 with missing adjacent rates
#' v <- c(0.9,0.1,0.3,0,0.8,0,0.1,0.1,0.7,0.6,0.2,0.6,0,0.5,0,0.4,0.3,0.4)
#' x <- array(v, dim = c(3,3,2,1))
#' m <- array(c(0,1,1,0,0,0,1,1,0), dim = c(3,3))
#' a <- create_movement_parameters(x, m)
#' k <- create_movement_rates(a, m)
#' all(k - x < 1e-12)
#'
create_movement_parameters <- function (x,
                                        m = NULL) {

  # Check arguments
  checkmate::assert_array(x, d = 4, any.missing = FALSE)
  checkmate::assert_matrix(m, mode = "integerish", any.missing = FALSE)
  # Initialize index limits
  np <- sum(m)
  na <- dim(x)[1]
  npt <- dim(x)[3]
  ng <- dim(x)[4]
  # Check x against m
  for (mg in seq_len(ng)) {
    for (cpt in seq_len(npt)) {
      for (ca in seq_len(na)) {
        for (pa in seq_len(na)) {
          if (pa != ca) {
            if (m[pa, ca] == 0) {
              if (x[pa, ca, npt, ng] != 0) {
                stop("movement rates must correspond to parameter index")
              }
            }
          }
        }
      }
    }
  }
  # Initialize array
  ta <- array(0, dim = c(np, npt, ng))
  # Compute movement parameters from movement rates
  for (mg in seq_len(ng)) {
    for (cpt in seq_len(npt)) {
      # Set the current parameter index to one
      cp_ind <- 1L
      # Iterate over rows
      for (pa in seq_len(na)) {
        # Set the movement rate sum to zero
        k_sum <- 0L
        # Iterate over columns
        for (ca in seq_len(na)) {
          if (m[pa, ca] != 0) {
            k_sum <- k_sum + x[pa, ca, cpt, mg]
          }
        }
        # Compute the movement parameters
        for (ca in seq_len(na)) {
          if (m[pa, ca] != 0) {
            ta[cp_ind, cpt, mg] <- log(x[pa, ca, cpt, mg]) - log(1 - k_sum)
            cp_ind <- cp_ind + 1L
          }
        }
      }
    }
  }
  # Untranspose
  a <- aperm(ta, c(2, 1, 3))
  # Return
  return(a)
}

#' Subset Numeric Vector or Matrix by Rowname, Colname, or Name
#'
#' @param x [numeric()] or [matrix()]
#' @param x_names [character()]
#'
#' @return
#'
#' @examples
#'
subset_by_name <- function (x, x_names) {
  if (is.numeric(x) & !is.matrix(x)) {
    x[which(is.element(names(x), x_names))]
  } else if (is.matrix(x)) {
    rows <- which(is.element(rownames(x), x_names))
    cols <- which(is.element(colnames(x), x_names))
    x[rows, cols]
  }
}

#' Matrix Power
#'
#' @param x [matrix()] Numeric matrix
#' @param n [integer()] Matrix power
#'
#' @return [matrix()]
#' @export
#'
#' @examples
#'
#' x <- matrix(c(1,2,3,4), nrow = 2)
#' n <- 2
#' matrix_power(x, n)
#'
matrix_power <- function (x, n) {
  # Check arguments
  checkmate::assert_matrix(x, mode = "numeric", any.missing = FALSE)
  checkmate::assert_integerish(n, lower = 1, len = 1, any.missing = FALSE)
  # Compute matrix power
  x_n <- x
  for (i in seq_len(n - 1)) x_n <- x_n %*% x
  return(x_n)
}

#' Create Movement Results
#'
#' @param v [numeric()] Vector of movement parameters in transpose order.
#' @param m [matrix()] Matrix of movement parameter covariances in
#'   transpose order.
#' @param np [integer()] Number of movement parameters per step and group.
#' @param npt [integer()] Number of movement parameter time steps.
#' @param ng [integer()] Number of groups.
#' @param mI [matrix()] Matrix index. See [mmmIndex()].
#' @param pow [integer()] Results step as a multiple of tag time step.
#'   A matrix power.
#' @param draws [integer()] Number of draws from which to bootstrap
#' movement result standard errors.
#'
#' @return A [list()]
#' @export
#'
#' @examples
#'
create_movement_results <- function (v, m, np, npt, ng, mI, pow, draws) {

  # Check arguments ------------------------------------------------------------

  # Vector v: vtaP
  checkmate::assert_numeric(v, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(v, len = np * npt * ng)
  # Matrix m: mtaP_cov
  checkmate::assert_matrix(m, any.missing = FALSE)
  checkmate::assert_true(nrow(m) == ncol(m))
  # Dimensions
  checkmate::assert_integerish(np, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(npt, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(ng, lower = 1, len = 1, any.missing = FALSE)
  # Matrix index
  checkmate::assert_matrix(mI, mode = "integerish", any.missing = FALSE)
  checkmate::assert_numeric(mI, lower = 0, upper = 1)
  checkmate::assert_true(nrow(mI) == ncol(mI))
  # Matrix power: result step as a multiple of tag step
  checkmate::assert_integerish(pow, lower = 1, len = 1, any.missing = FALSE)
  # Number of bootstrapping draws
  checkmate::assert_integerish(draws, lower = 1, len = 1, any.missing = FALSE)

  # Create movement parameters -------------------------------------------------

  taP_fit <- array(v, dim = c(np, npt, ng))
  aP_fit <- aperm(taP_fit, c(2, 1, 3))

  # Create movement rates ------------------------------------------------------

  aK_fit <- create_movement_rates(a = aP_fit, m = mI)
  aK_results <- array(0, dim = dim(aK_fit))
  for (cpt in seq_len(npt)) {
    for (mg in seq_len(ng)) {
      aK_results[, , cpt, mg] <- matrix_power(aK_fit[, , cpt, mg], pow)
    }
  }

  # Create movement rate standard errors ---------------------------------------

  # Initialize
  aK_draws <- array(0, dim = c(dim(aK_results), draws))
  aK_se <- array(0, dim = dim(aK_results))
  # Perform draws: each row is one draw in vector transpose order
  mtaP_draws <- MASS::mvrnorm(
    n = draws,
    mu = v,
    Sigma = m,
    empirical = FALSE
  )
  # Populate aK_draws
  for (i in seq_len(draws)) {
    # Create parameter array
    tap <- array(mtaP_draws[i, ], dim = c(np, npt, ng))
    ap <- aperm(tap, c(2, 1, 3))
    # Create movement rates
    aK_draws[, , , , i] <- create_movement_rates(a = ap, m = mI)
  }
  # Summarize SE
  na <- nrow(mI)
  for (mg in seq_len(ng)) {
    for (cpt in seq_len(npt)) {
      for (ca in seq_len(na)) {
        for (pa in seq_len(na)) {
          aK_se[pa, ca, cpt, mg] <- sd(aK_draws[pa, ca, cpt, mg, ], na.rm = TRUE)
        }
      }
    }
  }

  # Create movement results data frame -----------------------------------------

  # Initialize
  movement_results <- matrix(0, nrow = prod(dim(aK_results)), ncol = 6L)
  row_ind <- 1L
  # Populate
  for (mg in seq_len(ng)) {
    for (cpt in seq_len(npt)) {
      for (pa in seq_len(na)) {
        for (ca in seq_len(na)) {
          movement_results[row_ind, 1] <- pa
          movement_results[row_ind, 2] <- ca
          movement_results[row_ind, 3] <- cpt
          movement_results[row_ind, 4] <- mg
          movement_results[row_ind, 5] <- aK_results[pa, ca, cpt, mg]
          movement_results[row_ind, 6] <- aK_se[pa, ca, cpt, mg]
          row_ind <- row_ind + 1L
        }
      }
    }
  }
  # Set column names
  dim_names <- c("Area_From", "Area_To", "Result_Step", "Group")
  res_names <- c("Estimate", "SE")
  colnames(movement_results) <- c(dim_names, res_names)
  # Convert to data frame
  movement_results <- as.data.frame(movement_results)

  # Return list of results -----------------------------------------------------

  return(list(
    movement_results = movement_results,
    aP_fit = aP_fit,
    aK_fit = aK_fit,
    aK_results = aK_results,
    aK_se = aK_se
  ))
}
