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
#' @param a Numeric array. Movement parameters.
#' @param m Integer matrix. Movement indexes.
#' @param t Are \code{a} and \code{m} transposes?
#'
#' @return
#' @export
#'
#' @examples
#'
create_movement_rates <- function (a, m, t = FALSE) {

  # Check arguments
  checkmate::assert_array(a, mode = "numeric", any.missing = FALSE)
  checkmate::assert_array(a, d = 3)
  checkmate::assert_matrix(m, mode = "integerish", any.missing = FALSE)
  checkmate::assert_numeric(m, lower = 0, upper = 1)
  # Convert to transposes
  if (t) {
    ta <- a
    tm <- m
  } else {
    ta <- aperm(a, c(2,1,3))
    tm <- t(m)
  }
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

#' Generic Create Movement Parameter Array
#'
#' @param x A vector of parameter values or array of movement rates
#' @param ... Other arguments for class specific methods
#'
#' @return An array
#' @export
#'
#' @examples
#'
create_movement_parameters <- function (x, ...) {
  UseMethod("create_movement_parameters")
}

#' Default Create Movement Parameter Array
#'
#' @param x A vector of parameter values or array of movement rates
#' @param ... Other arguments for class specific methods
#'
#' @return
#'
#' @examples
#'
create_movement_parameters.default <- function (x, ...) {
  stop("methods only implemented for numeric and array classes")
}

#' Create Movement Parameter Array from a Vector of Parameters
#'
#' @param x Numeric vector. Movement parameters.
#' @param d Integer vector. Dimensions based on source of \code{x}. From
#'   \code{R} use \code{c(npt, np, ng)} or from \code{TMB} use
#'   \code{c(np, npt, ng)}.
#' @param t Logical. Transpose the first two dimensions of the parameter
#'   array?
#'
#' @return An array.
#'
#' @examples
#'
#' # aP
#' x <- rnorm(24, 0, 1)
#' d <- c(6, 4, 2)
#' aP <- create_movement_parameters(x, d)
#'
#' # taP
#' x <- rnorm(24, 0, 1)
#' d <- c(6, 4, 2)
#' taP <- create_movement_parameters(x, d, TRUE)
#'
create_movement_parameters.numeric <- function (x,
                                                d = NULL,
                                                t = FALSE) {
  # Check arguments
  checkmate::assert_numeric(x, finite = TRUE, any.missing = FALSE)
  checkmate::assert_integerish(d, lower = 0, len = 3, any.missing = FALSE)
  checkmate::assert_logical(t, len = 1L, any.missing = FALSE)
  # Assemble array
  a <- array(x, dim = d)
  # Transpose?
  if (t) a <- aperm(a, c(2, 1, 3))
  # Return
  return(a)
}

#' Create Movement Parameter Array from a Movement Rate Array
#'
#' @param x Numeric Array. aK.
#' @param m Integer matrix. Movement indexes. dim = c(pa, ca)
#' @param t Logical. Transpose the first two dimensions of the parameter
#'   array?
#'
#' @return An array
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
create_movement_parameters.array <- function (x,
                                              m = NULL,
                                              t = FALSE) {

  # Check arguments
  checkmate::assert_array(x, d = 4, any.missing = FALSE)
  checkmate::assert_null(d)
  checkmate::assert_logical(t, len = 1L, any.missing = FALSE)
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
  # Untranspose?
  if (t) {
    a <- ta
  } else {
    a <- aperm(ta, c(2, 1, 3))
  }
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

matrix_power <- function (x, n) {
  x_n <- x
  for (i in seq_len(n - 1)) {x_n <- x_n %*% x}
  x_n
}

create_movement_parameters_3d <- function (pars, nv, np, ng) {
  # Create movement parameter array
  array(pars, dim = c(nv, np, ng))
}

create_movement_probability_4d <- function(na, nv, ng, mp_3d, tp_2d) {

  #---------------- Check arguments -------------------------------------------#

  #---------------- Populate the movement probability array -------------------#

  # Initialize array
  mp_4d <- array(0, dim = c(na, na, nv, ng))
  # Compute probabilities
  for (mg in seq_len(ng)) {
    for (cv in seq_len(nv)) {
      par_num <- 1
      par_denom <- 1
      for (pa in seq_len(na)) {
        # Set the denom sum to zero for this row
        denom_sum <- 0
        # Compute the denom sum
        for (ca in seq_len(na)) {
          if (tp_2d[pa, ca] > 0) {
            denom_sum <- denom_sum + exp(mp_3d[cv, par_denom, mg])
            par_denom <- par_denom + 1
          }
        }
        # Set the prob sum to zero for this row
        prob_sum <- 0
        # Compute the probability for pa != ca
        for (ca in seq_len(na)) {
          if (tp_2d[pa, ca] > 0) {
            mp_4d[pa, ca, cv, mg] <- exp(mp_3d[cv, par_num, mg]) / (1 + denom_sum)
            prob_sum <- prob_sum + mp_4d[pa, ca, cv, mg]
            par_num <- par_num + 1
          }
        }
        # Compute the probability for pa == ca
        mp_4d[pa, pa, cv, mg] <- 1 - prob_sum
      }
    }
  }
  # Return the movement probability array
  mp_4d
}


create_movement_probability_results_4d <- function (mp_4d, result_units) {

  #---------------- Check arguments -------------------------------------------#

  #---------------- Extract dimensions ----------------------------------------#

  nv <- dim(mp_4d)[3]
  ng <- dim(mp_4d)[4]

  #---------------- Populate the movement probability results array -----------#

  mpr_4d <- array(0, dim = dim(mp_4d))
  for (mg in seq_len(ng)) {
    for (cv in seq_len(nv)) {
      mpr_4d[, , cv, mg] <- matrix_power(mp_4d[, , cv, mg], result_units)
    }
  }
  # Return
  mpr_4d
}



create_movement_probability_results_std_err_4d <- function (pars,
                                                            covs,
                                                            dims,
                                                            tp_2d,
                                                            result_units,
                                                            n_draws = 1000) {

  #---------------- Check arguments -------------------------------------------#



  #---------------- Extract dimensions ----------------------------------------#

  nv <- dims[1]
  np <- dims[2]
  nt <- dims[3]
  na <- dims[4]
  ng <- dims[5]

  #---------------- Perform draws ---------------------------------------------#

  mpr_se_5d <- array(0, dim = c(na, na, nv, ng, n_draws))
  for (i in seq_len(n_draws)) {

    #---------------- Draw parameters -----------------------------------------#

    par_devs <- MASS::mvrnorm(n = 1,
                              mu = rep(0, nv * np * ng),
                              Sigma = covs,
                              empirical = FALSE)
    par_draw <- pars + par_devs

    #---------------- Create movement parameter array -------------------------#

    mp_3d <- create_movement_parameters_3d(pars = par_draw,
                                           nv = nv,
                                           np = np,
                                           ng = ng)

    #---------------- Create movement probability array -----------------------#

    mp_4d <- create_movement_probability_4d(na = na,
                                            nv = nv,
                                            ng = ng,
                                            mp_3d = mp_3d,
                                            tp_2d = tp_2d)

    #---------------- Create movement probability results array ---------------#

    mpr_se_5d[, , , , i] <- create_movement_probability_results_4d(
      mp_4d = mp_4d,
      result_units = result_units
    )
  }

  #---------------- Summarize to SE -------------------------------------------#

  mpr_se_4d <- array(0, dim = c(na, na, nv, ng))
  for (pa in seq_len(na)) {
    for (ca in seq_len(na)) {
      for (cv in seq_len(nv)) {
        for (mg in seq_len(ng)) {
          mpr_se_4d[pa, ca, cv, mg] <- sd(mpr_se_5d[pa, ca, cv, mg, ], na.rm = TRUE)
        }
      }
    }
  }
  # Return
  mpr_se_4d
}



create_movement_probability_results <- function (pars,
                                                 covs,
                                                 dims,
                                                 tp_2d,
                                                 result_units,
                                                 n_draws = 1000) {

  #---------------- Check arguments -------------------------------------------#



  #---------------- Extract dimensions ----------------------------------------#

  nv <- dims[1]
  np <- dims[2]
  nt <- dims[3]
  na <- dims[4]
  ng <- dims[5]

  #---------------- Create movement parameter array ---------------------------#

  mp_3d <- create_movement_parameters_3d(pars = pars,
                                         nv = nv,
                                         np = np,
                                         ng = ng)

  #---------------- Create movement probability array -------------------------#

  mp_4d <- create_movement_probability_4d(na = na,
                                          nv = nv,
                                          ng = ng,
                                          mp_3d = mp_3d,
                                          tp_2d = tp_2d)

  #---------------- Create movement probability results array -----------------#

  mpr_4d <- create_movement_probability_results_4d(
    mp_4d = mp_4d,
    result_units = result_units
  )

  #---------------- Compute standard errors array -----------------------------#

  mpr_se_4d <- create_movement_probability_results_std_err_4d(
    pars = pars,
    covs = covs,
    dims = dims,
    tp_2d = tp_2d,
    result_units = result_units,
    n_draws = n_draws
  )

  #---------------- Create movement probability results data frame ------------#

  mp_results_mat <- matrix(0, nrow = na * na * nv * ng, ncol = 6)
  row_count <- 1
  for (mg in seq_len(ng)) {
    for (pa in seq_len(na)) {
      for (ca in seq_len(na)) {
        for (cv in seq_len(nv)) {
          mp_results_mat[row_count, 1] <- cv
          mp_results_mat[row_count, 2] <- mg
          mp_results_mat[row_count, 3] <- pa
          mp_results_mat[row_count, 4] <- ca
          mp_results_mat[row_count, 5] <- mpr_4d[pa, ca, cv, mg]
          mp_results_mat[row_count, 6] <- mpr_se_4d[pa, ca, cv, mg]
          row_count <- row_count + 1
        }
      }
    }
  }
  # Set column names
  colnames(mp_results_mat) <- c("Pattern_Time", "Class", "Area_From",
                                "Area_To", "Estimate", "SE")
  # Convert to data frame
  mp_results_df <- as.data.frame(mp_results_mat)

  # Return movement probability results data frame
  mp_results_df
}



