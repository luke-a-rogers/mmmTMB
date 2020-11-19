#' Check Data Argument for mmmFit()
#'
#' @param data [list()]
#'
#' @return
#'
#' @examples
#'
check_data <- function (data = NULL) {

  # Check list and elements ----------------------------------------------------

  checkmate::assert_list(data)
  checkmate::assert_matrix(data$x, mode = "integerish")
  checkmate::assert_matrix(data$y, mode = "integerish")
  checkmate::assert_matrix(data$z, mode = "integerish")
  checkmate::assert_matrix(data$l, mode = "double", null.ok = TRUE)
  checkmate::assert_matrix(data$w, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(data$x, lower = 0)
  checkmate::assert_numeric(data$y, lower = 0)
  checkmate::assert_numeric(data$z, lower = 0, upper = 1)
  checkmate::assert_numeric(data$l, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_numeric(data$w, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_numeric(data$h, lower = 0, len = 1)
  checkmate::assert_numeric(data$u, lower = 0, upper = 1, len = 1)
  # Optionally parameters
  checkmate::assert_matrix(data$f, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(data$f, lower = 0, null.ok = TRUE)
  checkmate::assert_numeric(data$m, lower = 0, len = 1, null.ok = TRUE)

  # Check f xor m null ---------------------------------------------------------

  if (is.null(data$f) && is.null(data$m)) {
    stop("matrix f or scalar m must be defined in data")
  }

  #---------------- Check required data ---------------------------------------#

  # Tag releases
  checkmate::assert_matrix(data$x, mode = "integerish", any.missing = FALSE)
  checkmate::assert_matrix(data$x, ncols = 4, null.ok = FALSE)
  checkmate::assert_numeric(data$x, lower = 0, null.ok = FALSE)
  checkmate::assert_true(colnames(data$x)[1] == "release_step")
  checkmate::assert_true(colnames(data$x)[2] == "release_area")
  checkmate::assert_true(colnames(data$x)[3] == "group")
  checkmate::assert_true(colnames(data$x)[4] == "count")
  # Tag recoveries
  checkmate::assert_matrix(data$y, mode = "integerish", any.missing = FALSE)
  checkmate::assert_matrix(data$y, ncols = 6, null.ok = FALSE)
  checkmate::assert_numeric(data$y, lower = 0, null.ok = FALSE)
  checkmate::assert_true(colnames(data$y)[1] == "release_step")
  checkmate::assert_true(colnames(data$y)[2] == "release_area")
  checkmate::assert_true(colnames(data$y)[3] == "group")
  checkmate::assert_true(colnames(data$y)[4] == "recover_step")
  checkmate::assert_true(colnames(data$y)[5] == "recover_area")
  checkmate::assert_true(colnames(data$y)[6] == "count")
  # Index matrix
  checkmate::assert_true(all(diag(data$z) == 0L))
  # Check values
  if (min(data$x[, "release_step"]) > 0) {
    cat("caution: x time step not indexed from zero")
  }

  # Check index matrix ---------------------------------------------------------

  checkmate::assert_matrix(data$z, mode = "integerish")
  checkmate::assert_matrix(data$z, any.missing = FALSE)
  checkmate::assert_numeric(data$z, lower = 0, upper = 1)
  checkmate::assert_true(nrow(data$z) == ncol(data$z))

  # Set default data -----------------------------------------------------------

  if (is.null(data$l)) l <- matrix(1, nrow = 1L, ncol = 1L)
  if (is.null(data$w)) w <- matrix(1, nrow = 1L, ncol = 1L)

  # Check data arguments -------------------------------------------------------

  # Tag reporting rate
  checkmate::assert_matrix(l, mode = "double")
  checkmate::assert_matrix(l, any.missing = FALSE)
  checkmate::assert_numeric(l, lower = 0)
  # checkmate::assert_true(ncol(l) == 1 || ncol(l) == na)
  # Fishing rate weighting
  checkmate::assert_matrix(w, mode = "double")
  checkmate::assert_matrix(w, any.missing = FALSE)
  checkmate::assert_numeric(w, lower = 0)
  # checkmate::assert_true(ncol(w) == 1 || ncol(w) == na)
  checkmate::assert_true(all(colSums(w) == 1))

}

#' Check Parameters Argument for mmmFit()
#'
#' @param parameters [list()]
#'
#' @return
#'
#' @examples
#'
check_parameters <- function (parameters) {

  # Check parameters argument --------------------------------------------------

  checkmate::assert_list(parameters, null.ok = TRUE)
  # Optional initial values
  checkmate::assert_array(parameters$p, mode = "double", null.ok = TRUE)
  checkmate::assert_array(parameters$p, d = 3L, null.ok = TRUE)
  # Optionally data
  checkmate::assert_matrix(parameters$f, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(parameters$f, lower = 0, null.ok = TRUE)
  checkmate::assert_double(parameters$m, lower = 0, len = 1, null.ok = TRUE)

  #---------------- Check arguments for tmb_parameters ------------------------#

  # # Movement parameter array
  # checkmate::assert_array(p, mode = "double", any.missing = FALSE, d = 3)
  # checkmate::assert_true(dim(p)[1] == npt, na.ok = FALSE)
  # checkmate::assert_true(dim(p)[2] == np, na.ok = FALSE)
  # checkmate::assert_true(dim(p)[3] == ng, na.ok = FALSE)
  # # Fishing rate matrix
  # checkmate::assert_matrix(f, mode = "double", any.missing = FALSE)
  # checkmate::assert_numeric(f, lower = 0, null.ok = TRUE)
  # checkmate::assert_true(dim(f)[1] == nft, na.ok = FALSE)
  # checkmate::assert_true(dim(f)[2] == nfa, na.ok = FALSE)
}

#' Check Settings Argument for mmmFit()
#'
#' @param settings [list()]
#'
#' @return
#'
#' @examples
#'
check_settings <- function (settings) {

  # Check settings -------------------------------------------------------------

  # Check settings
  checkmate::assert_list(settings, null.ok = FALSE)
  # Check error family
  checkmate::assert_integerish(settings$error_family, len = 1)
  checkmate::assert_integerish(settings$error_family, lower = 0, upper = 1)
  checkmate::assert_integerish(settings$error_family, any.missing = FALSE)
  # Check max liberty
  checkmate::assert_integerish(settings$max_liberty, lower = 0, len = 1)
  checkmate::assert_integerish(settings$max_liberty, any.missing = FALSE)
  # Check time varying
  checkmate::assert_integerish(settings$time_varying, len = 1)
  checkmate::assert_integerish(settings$time_varying, lower = 0, upper = 1)
  checkmate::assert_integerish(settings$time_varying, any.missing = FALSE)
  # Check time process
  checkmate::assert_integerish(settings$time_process, len = 1)
  checkmate::assert_integerish(settings$time_process, lower = 0, upper = 1)
  checkmate::assert_integerish(settings$time_process, any.missing = FALSE)
  # Check cycle length
  checkmate::assert_integerish(settings$cycle_length, len = 1)
  checkmate::assert_integerish(settings$cycle_length, lower = 0)
  checkmate::assert_integerish(settings$cycle_length, any.missing = FALSE)
  # Check block length
  checkmate::assert_integerish(settings$block_length, len = 1)
  checkmate::assert_integerish(settings$block_length, lower = 0)
  checkmate::assert_integerish(settings$block_length, any.missing = FALSE)
  # Check results step
  checkmate::assert_integerish(settings$results_step, len = 1)
  checkmate::assert_integerish(settings$results_step, lower = 1)
  checkmate::assert_integerish(settings$results_step, any.missing = FALSE)
  # Check nlminb loops
  checkmate::assert_integerish(settings$nlminb_loops, len = 1)
  checkmate::assert_integerish(settings$nlminb_loops, lower = 0)
  checkmate::assert_integerish(settings$nlminb_loops, any.missing = FALSE)
  # Check newton iterations
  checkmate::assert_integerish(settings$newton_iters, len = 1)
  checkmate::assert_integerish(settings$newton_iters, lower = 0)
  checkmate::assert_integerish(settings$newton_iters, any.missing = FALSE)
  # Check OpenMP cores
  checkmate::assert_integerish(settings$openmp_cores, len = 1)
  checkmate::assert_integerish(settings$openmp_cores, lower = 1)
  checkmate::assert_integerish(settings$openmp_cores, any.missing = FALSE)
}

#' Get Convergence Diagnostics from TMB Model Fit
#'
#' @description Get convergence diagnostics from a TMB model fit.
#'   This function is from the sdmTMB' package at
#'   github.com/pbs-assess/sdmTMB/blob/master/R/fit.R
#'
#' @param sd_report [sdreport()] See [TMB::sdreport()].
#'
#' @return [list()]
#'
#' @examples
#'
get_convergence_diagnostics <- function(sd_report) {
  final_grads <- sd_report$gradient.fixed
  bad_eig <- FALSE
  if (!is.null(sd_report$pdHess)) {
    if (!sd_report$pdHess) {
      warning("The model may not have converged: ",
              "non-positive-definite Hessian matrix.", call. = FALSE)
    } else {
      eigval <- try(1 / eigen(sd_report$cov.fixed)$values, silent = TRUE)
      if (is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
        warning("The model may not have converged: ",
                "extreme or very small eigen values detected.", call. = FALSE)
        bad_eig <- TRUE
      }
      if (any(abs(final_grads) > 0.01))
        warning("The model may not have converged. ",
                "Maximum final gradient: ", max(abs(final_grads)), ".", call. = FALSE)
    }
  }
  invisible(list(final_grads = final_grads, bad_eig = bad_eig))
}
