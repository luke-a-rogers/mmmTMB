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



