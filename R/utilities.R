#' Logit (Log Odds Ratio)
#'
#' @param x [numeric()] Value between zero and one (inclusive).
#'
#' @return [numeric()] The log of the odds ratio.
#' @export
#'
#' @examples
#' logit(0)
#' logit(0.5)
#' logit(1)
#' logit(matrix(c(0, 0.5, 0.5, 1), nrow = 2))
#'
logit <- function (x) {
  # Check arguments
  checkmate::assert_numeric(x, lower = 0, upper = 1)
  # Return logit(x)
  return(log(x / (1 - x)))
}

#' Inverse Logit
#'
#' @param y [numeric()] Finite or infinite scalar.
#'
#' @return The inverse of the log odds ratio
#' @export
#'
#' @examples
#' invlogit(Inf)
#' invlogit(0)
#' invlogit(-Inf)
#' invlogit(matrix(c(-Inf, 0, 0, Inf), nrow = 2))
#'
invlogit <- function (y) {
  # Check arguments
  checkmate::assert_numeric(y)
  # Compute invlogit(y)
  ans <- numeric(length = length(y))
  inf_ind <- which(y == Inf)
  fin_ind <- which(y != Inf)
  ans[inf_ind] <- 1
  ans[fin_ind] <- exp(y[fin_ind]) / (1 + exp(y[fin_ind]))
  dim(ans) <- dim(y)
  # Return invlogit(y)
  return(ans)
}

#' Create Index Vector
#'
#' @param a [integer()] Length of object to index
#' @param b [integer()] Length of index vector
#' @param from [integer()] Value from which to index
#'
#' @return An [integer()] vector
#' @export
#'
#' @examples
#' # a == 1
#' create_index_vector(1, 3) # Repeated
#'
#' # a == b
#' create_index_vector(3, 3) # Sequence
#'
#' # a < b
#' create_index_vector(2, 3) # Truncated sequence of repeated values
#' create_index_vector(2, 4) # Sequence of repeated values
#' create_index_vector(5, 9) # Truncated sequence of repeated values
#'
#' # a > b
#' create_index_vector(5, 3) # Truncated sequence
#'
create_index_vector <- function (a, b, from = 0) {
  # Check arguments
  checkmate::assert_integerish(a, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(b, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(from, lower = 0, len = 1, any.missing = FALSE)
  # Compute index vector
  vec <- rep(c(seq_len(a) - 1L), each = ceiling(b / a))[seq_len(b)]
  vec <- as.integer(vec + from)
  # Return
  return(vec)
}


#' Create Movement Rate Array
#'
#' @param p [array()] Movement parameters.
#' @param z [matrix()] Movement index..
#'
#' @return [array()] Movement rates.
#' @export
#'
#' @examples
#' # 2x2
#' x <- array(c(0.9, 0.2, 0.1, 0.8, 0.7, 0.4, 0.3, 0.6), dim = c(2, 2, 2, 1))
#' z <- array(c(0, 1, 1, 0), dim = c(2, 2))
#' p <- create_movement_parameters(x, z)
#' r <- create_movement_rates(p, z)
#' all(x - r < 1e-12)
#'
#' # 3x3 with missing adjacent rates
#' v <- c(0.9,0.1,0.3,0,0.8,0,0.1,0.1,0.7,0.6,0.2,0.6,0,0.5,0,0.4,0.3,0.4)
#' x <- array(v, dim = c(3,3,2,1))
#' z <- array(c(0,1,1,0,0,0,1,1,0), dim = c(3,3))
#' p <- create_movement_parameters(x, z)
#' r <- create_movement_rates(p, z)
#' all(x - r < 1e-12)
#'
create_movement_rates <- function (p, z) {

  # Check arguments
  checkmate::assert_array(p, mode = "numeric", any.missing = FALSE)
  checkmate::assert_array(p, d = 3)
  checkmate::assert_matrix(z, mode = "integerish", any.missing = FALSE)
  checkmate::assert_numeric(z, lower = 0, upper = 1)
  # Convert to transposes
  tp <- aperm(p, c(2,1,3))
  tz <- t(z)
  # Initialize index limits
  na <- nrow(tz) # Matrix tz has rows = na, and cols = na
  npt <- dim(tp)[2] # Array tp has dim = c(np, npt, ng)
  ng <- dim(tp)[3]
  # Initialize array
  r <- array(0, dim = c(na, na, npt, ng))
  # Compute movement rates from movement parameters
  for (mg in seq_len(ng)) {
    for (cpt in seq_len(npt)) {
      # Set the current parameter indexes to one
      cp_num <- 1L
      cp_den <- 1L
      for (pa in seq_len(na)) {
        # Set the denominator sum to zero for this row
        den_sum <- 0
        for (ca in seq_len(na)) {
          if (tz[ca, pa] != 0) {
            den_sum <- den_sum + exp(tp[cp_den, cpt, mg])
            cp_den <- cp_den + 1
          }
        }
        # Set the probability sum to zero
        rates_sum <- 0L
        # Compute the probability for ca != pa
        for (ca in seq_len(na)) {
          if (tz[ca, pa] != 0) {
            r[pa, ca, cpt, mg] <- exp(tp[cp_num, cpt, mg]) / (1L + den_sum)
            rates_sum <- rates_sum + r[pa, ca, cpt, mg]
            cp_num <- cp_num + 1
          }
        }
        # Compute the probability for pa == ca
        r[pa, pa, cpt, mg] <- 1L - rates_sum
      }
    }
  }
  return(r)
}

#' Create Movement Parameter Array
#'
#' @param r [array()] Movement rates.
#' @param z [matrix()] Movement index.
#'
#' @return [array()] Movement parameters.
#' @export
#'
#' @examples
#'
#' # 2x2
#' x <- array(c(0.9, 0.2, 0.1, 0.8, 0.7, 0.4, 0.3, 0.6), dim = c(2, 2, 2, 1))
#' z <- array(c(0, 1, 1, 0), dim = c(2, 2))
#' p <- create_movement_parameters(x, z)
#' r <- create_movement_rates(p, z)
#' all(x - r < 1e-12)
#'
#' # 3x3 with missing adjacent rates
#' v <- c(0.9,0.1,0.3,0,0.8,0,0.1,0.1,0.7,0.6,0.2,0.6,0,0.5,0,0.4,0.3,0.4)
#' x <- array(v, dim = c(3,3,2,1))
#' z <- array(c(0,1,1,0,0,0,1,1,0), dim = c(3,3))
#' p <- create_movement_parameters(x, z)
#' r <- create_movement_rates(p, z)
#' all(x - r < 1e-12)
#'
create_movement_parameters <- function (r, z = NULL) {

  # Check arguments
  checkmate::assert_array(r, d = 4, any.missing = FALSE)
  checkmate::assert_numeric(r, lower = 0, upper = 1)
  checkmate::assert_matrix(z, mode = "integerish", any.missing = FALSE)
  # Initialize index limits
  np <- sum(z)
  na <- dim(r)[1]
  npt <- dim(r)[3]
  ng <- dim(r)[4]
  # Check r against z
  for (mg in seq_len(ng)) {
    for (cpt in seq_len(npt)) {
      row_sums <- rowSums(r[,, cpt, mg])
      if (!all(row_sums == 1)) stop("row sums of r must all equal one\n")
      for (ca in seq_len(na)) {
        for (pa in seq_len(na)) {
          if (pa != ca) {
            if (z[pa, ca] == 0) {
              if (r[pa, ca, cpt, mg] != 0) {
                stop("movement rates must correspond to parameter index")
              }
            }
          }
        }
      }
    }
  }
  # Initialize array
  tp <- array(0, dim = c(np, npt, ng))
  # Compute movement parameters from movement rates
  for (mg in seq_len(ng)) {
    for (cpt in seq_len(npt)) {
      # Set the current parameter index to one
      cp_ind <- 1L
      # Iterate over rows
      for (pa in seq_len(na)) {
        # Set the movement rate sum to zero
        rates_sum <- 0L
        # Iterate over columns
        for (ca in seq_len(na)) {
          if (z[pa, ca] != 0) {
            rates_sum <- rates_sum + r[pa, ca, cpt, mg]
          }
        }
        # Compute the movement parameters
        for (ca in seq_len(na)) {
          if (z[pa, ca] != 0) {
            tp[cp_ind, cpt, mg] <- log(r[pa, ca, cpt, mg]) - log(1 - rates_sum)
            cp_ind <- cp_ind + 1L
          }
        }
      }
    }
  }
  # Untranspose
  p <- aperm(tp, c(2, 1, 3))
  # Return
  return(p)
}

#' Subset Numeric Vector or Matrix by Rowname, Colname, or Name
#'
#' @param x [numeric()] or [matrix()]
#' @param x_names [character()]
#'
#' @return
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
#' @param x [matrix()] Square numeric matrix.
#' @param pow [integer()] Matrix power.
#'
#' @return [matrix()]
#' @export
#'
#' @examples
#'
#' x <- matrix(c(1,2,3,4), nrow = 2)
#' pow <- 2
#' matrix_power(x, pow)
#'
matrix_power <- function (x, pow) {
  # Check arguments
  checkmate::assert_matrix(x, mode = "numeric", any.missing = FALSE)
  checkmate::assert_true(nrow(x) == ncol(x))
  checkmate::assert_integerish(pow, lower = 1, len = 1, any.missing = FALSE)
  # Compute matrix power
  x_pow <- x
  for (i in seq_len(pow - 1)) x_pow <- x_pow %*% x
  return(x_pow)
}

#' Create Movement Results
#'
#' @param vtp [numeric()] Vector of movement parameters in transpose order.
#' @param mtp [matrix()] Matrix of movement parameter covariances in
#'   transpose order.
#' @param np [integer()] Number of movement parameters per step and group.
#' @param npt [integer()] Number of movement parameter time steps.
#' @param ng [integer()] Number of groups.
#' @param z [matrix()] Matrix index. See [mmmIndex()].
#' @param pow [integer()] Results step as a multiple of tag time step.
#'   A matrix power.
#' @param draws [integer()] Number of draws from which to bootstrap
#' movement result standard errors.
#'
#' @return A [list()]
#'
create_movement_results <- function (vtp, mtp, np, npt, ng, z, pow, draws) {

  # Check arguments ------------------------------------------------------------

  # Vector vtp
  checkmate::assert_numeric(vtp, finite = TRUE, any.missing = FALSE)
  checkmate::assert_numeric(vtp, len = np * npt * ng)
  # Matrix mtp
  checkmate::assert_matrix(mtp, any.missing = FALSE)
  checkmate::assert_true(nrow(mtp) == ncol(mtp))
  # Dimensions
  checkmate::assert_integerish(np, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(npt, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(ng, lower = 1, len = 1, any.missing = FALSE)
  # Matrix index
  checkmate::assert_matrix(z, mode = "integerish", any.missing = FALSE)
  checkmate::assert_numeric(z, lower = 0, upper = 1)
  checkmate::assert_true(nrow(z) == ncol(z))
  # Matrix power: result step as a multiple of tag step
  checkmate::assert_integerish(pow, lower = 1, len = 1, any.missing = FALSE)
  # Number of bootstrapping draws
  checkmate::assert_integerish(draws, lower = 1, len = 1, any.missing = FALSE)

  # Create movement parameters -------------------------------------------------

  tp_fit <- array(vtp, dim = c(np, npt, ng))
  p_fit <- aperm(tp_fit, c(2, 1, 3))

  # Create movement rates ------------------------------------------------------

  r_fit <- create_movement_rates(p = p_fit, z = z)
  r_results <- array(0, dim = dim(r_fit))
  for (cpt in seq_len(npt)) {
    for (mg in seq_len(ng)) {
      r_results[, , cpt, mg] <- matrix_power(r_fit[, , cpt, mg], pow)
    }
  }

  # Create movement rate standard errors ---------------------------------------

  # Initialize
  r_draws <- array(0, dim = c(dim(r_results), draws))
  # Perform draws: each row is one draw in vector transpose order
  mtp_draws <- MASS::mvrnorm(
    n = draws,
    mu = vtp,
    Sigma = mtp,
    empirical = FALSE
  )
  # Populate r_draws
  for (i in seq_len(draws)) {
    # Create parameter array
    tp <- array(mtp_draws[i, ], dim = c(np, npt, ng))
    p <- aperm(tp, c(2, 1, 3))
    # Create movement rates
    r_draws[, , , , i] <- create_movement_rates(p = p, z = z)
  }
  # Convert to r_result_draws
  r_results_draws <- array(0, dim = c(dim(r_results), draws))
  r_results_se <- array(0, dim = dim(r_results))
  for (i in seq_len(draws)) {
    for (cpt in seq_len(npt)) {
      for (mg in seq_len(ng)) {
        r_results_draws[, , cpt, mg, i] <-
          matrix_power(r_draws[, , cpt, mg, i], pow)
      }
    }
  }
  # Summarize SE
  na <- nrow(z)
  for (mg in seq_len(ng)) {
    for (cpt in seq_len(npt)) {
      for (ca in seq_len(na)) {
        for (pa in seq_len(na)) {
          r_results_se[pa, ca, cpt, mg] <-
            stats::sd(r_results_draws[pa, ca, cpt, mg, ], na.rm = TRUE)
        }
      }
    }
  }

  # Create movement results data frame -----------------------------------------

  # Initialize
  movement_results <- matrix(0, nrow = prod(dim(r_results)), ncol = 6L)
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
          movement_results[row_ind, 5] <- r_results[pa, ca, cpt, mg]
          movement_results[row_ind, 6] <- r_results_se[pa, ca, cpt, mg]
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
  # Rename
  names(p_fit) <- NULL
  names(r_fit) <- NULL
  names(r_results) <- NULL
  names(r_results_se) <- NULL

  # Return list of results -----------------------------------------------------

  return(list(
    movement_rate = movement_results,
    p_fit = p_fit,
    r_fit = r_fit,
    r_results = r_results,
    r_results_se = r_results_se
  ))
}

#' Create Natural Mortality Results
#'
#' @param logit_exp_neg_m [numeric()] Scalar estimate.
#' @param logit_exp_neg_m_se [numeric()] Scalar estimate.
#' @param results_step [integer()] Results step as multiple of tag step.
#' @param estimate [logical()] Estimate natural mortality?
#'
#' @return [list()]
#'
create_mortality_results <- function (logit_exp_neg_m,
                                      logit_exp_neg_m_se,
                                      results_step,
                                      estimate) {
  if (estimate) {
    # Estimate
    m_fit <- -log(invlogit(logit_exp_neg_m))
    m_results <- m_fit * results_step
    # Standard error
    logit_exp_neg_m_draws <- stats::rnorm(
      1000,
      logit_exp_neg_m,
      logit_exp_neg_m_se
    )
    m_draws <- -log(invlogit(logit_exp_neg_m_draws))
    m_results_se <- stats::sd(m_draws * results_step, na.rm = TRUE)
    # Rename
    names(m_fit) <- NULL
    names(m_results) <- NULL
    names(m_results_se) <- NULL
  } else {
    m_fit <- NULL
    m_results <- NULL
    m_results_se <- NULL
  }
  # Data frame
  natural_mortality <- data.frame(Estimate = m_results, SE = m_results_se)
  # Return
  return(list(
    natural_mortality = natural_mortality,
    m_fit = m_fit,
    m_results = m_results,
    m_results_se = m_results_se
  ))
}


#' Create Fishing Rate Results
#'
#' @param vlogit_exp_neg_tf [numeric()] Vector estimates.
#' @param mlogit_exp_neg_tf_cov [matrix()] Covariances.
#' @param results_step [integer()] Results step as multiple of tag step.
#' @param nft [integer()] Number of fishing rate time steps.
#' @param nfa [integer()] Number of fishing rate areas.
#' @param estimate [logical()] Estimate fishing mortality rate(s)?
#'
#' @return [list()]
#'
create_fishing_results <- function (vlogit_exp_neg_tf,
                                    mlogit_exp_neg_tf_cov,
                                    results_step,
                                    nft,
                                    nfa,
                                    nfg,
                                    estimate) {
  if (estimate) {
    # Estimate
    vtf_fit <- -log(invlogit(vlogit_exp_neg_tf))
    tf_fit <- array(vtf_fit, dim = c(nfa, nft, nfg))
    f_fit <- aperm(tf_fit, c(2, 1, 3))
    f_results <- f_fit * results_step
    # Standard error
    mvlogit_exp_neg_tf_draws <- MASS::mvrnorm(
      n = 1000,
      mu = vlogit_exp_neg_tf,
      Sigma = mlogit_exp_neg_tf_cov
    )
    # Untransform
    mvtf_draws <- -log(invlogit(mvlogit_exp_neg_tf_draws))
    # Compute SE
    vtf_results_se <- apply(mvtf_draws * results_step, 2, stats::sd, na.rm = TRUE)
    tf_results_se <- array(vtf_results_se, dim = c(nfa, nft, nfg))
    f_results_se <- aperm(tf_results_se, c(2, 1, 3))
    # Rename
    names(f_fit) <- NULL
    names(f_results) <- NULL
    names(f_results_se) <- NULL
    # Initialize
    fishing_rate <- matrix(0, nrow = prod(dim(f_results)), ncol = 5L)
    row_ind <- 1L
    # Populate
    for (cfg in seq_len(nfg)) {
      for (cfa in seq_len(nfa)) {
        for (cft in seq_len(nft)) {
          fishing_rate[row_ind, 1] <- cfa
          fishing_rate[row_ind, 2] <- cft
          fishing_rate[row_ind, 3] <- cfg
          fishing_rate[row_ind, 4] <- f_results[cft, cfa, cfg]
          fishing_rate[row_ind, 5] <- f_results_se[cft, cfa, cfg]
          row_ind <- row_ind + 1L
        }
      }
    }
    # Set column names
    dim_names <- c("Fish_Area", "Fish_Step", "Fish_Group")
    res_names <- c("Estimate", "SE")
    colnames(fishing_rate) <- c(dim_names, res_names)
    # Convert to data frame
    fishing_rate <- as.data.frame(fishing_rate)
  } else {
    f_fit <- NULL
    f_results <- NULL
    f_results_se <- NULL
    # Data frame
    fishing_rate <- data.frame(
      Fish_Area = NULL,
      Fish_Step = NULL,
      Estimate = f_results,
      SE = f_results_se)
  }
  # Return
  return(list(
    fishing_rate = fishing_rate,
    f_fit = f_fit,
    f_results = f_results,
    f_results_se = f_results_se
  ))
}

#' Create Fishing Mortality Bias Results
#'
#' @param log_b [numeric()] Vector.
#' @param log_b_cov [matrix()] Covariances.
#' @param estimate [logical()] Estimate or return NULL.
#'
#' @return [list()]
#'
create_bias_results <- function (log_b,
                                 log_b_cov,
                                 estimate) {

  if (estimate) {
    # Estimate
    b_fit <- exp(log_b)
    b_results <- b_fit
    # Compute draws
    log_b_draws <- MASS::mvrnorm(
      n = 1000,
      mu = log_b,
      Sigma = log_b_cov
    )
    # Untransform
    b_draws <- exp(log_b_draws)
    # Compute SE
    b_results_se <- apply(b_draws, 2, stats::sd, na.rm = TRUE)
    # Rename
    names(b_fit) <- NULL
    names(b_results) <- NULL
    names(b_results_se) <- NULL
  } else {
    b_fit <- NULL
    b_results <- NULL
    b_results_se <- NULL
  }
  # Data frame
  fishing_bias <- data.frame(Estimate = b_results, SE = b_results_se)
  # Return
  return(list(
    fishing_bias = fishing_bias,
    b_fit = b_fit,
    b_results = b_results,
    b_results_se = b_results_se
  ))
}

#' Create Dispersion Results
#'
#' @param log_k [numeric()] Scalar estimate.
#' @param log_k_se [numeric()] Scalar estimate.
#' @param estimate [logical()] Estimate negative binomial dispersion?
#'
#' @return [list()]
#'
create_dispersion_results <- function (log_k,
                                       log_k_se,
                                       estimate) {
  if (estimate) {
    # Estimate
    k_fit <- exp(log_k)
    k_results <- k_fit
    # Compute SE
    log_k_draws <- stats::rnorm(1000, log_k, log_k_se)
    k_results_se <- stats::sd(exp(log_k_draws), na.rm = TRUE)
    # Rename
    names(k_fit) <- NULL
    names(k_results) <- NULL
    names(k_results_se) <- NULL
  } else {
    k_fit <- NULL
    k_results <- NULL
    k_results_se <- NULL
  }
  # Data frame
  dispersion <- data.frame(Estimate = k_results, SE = k_results_se)
  # Return
  return(list(
    dispersion = dispersion,
    k_fit = k_fit,
    k_results = k_results,
    k_results_se = k_results_se
  ))
}
