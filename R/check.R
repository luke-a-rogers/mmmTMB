#' Check Data Argument for mmmFit()
#'
#' @param data [list()]
#'
#' @return NULL
#'
check_data <- function (data = NULL) {

  # Check data -----------------------------------------------------------------

  checkmate::assert_list(data)

  # Check x --------------------------------------------------------------------

  checkmate::assert_matrix(data$x, mode = "integerish", any.missing = FALSE)
  checkmate::assert_matrix(data$x, ncols = 4)
  checkmate::assert_numeric(data$x, lower = 0, finite = TRUE)
  checkmate::assert_choice(colnames(data$x)[1], c("release_step"))
  checkmate::assert_choice(colnames(data$x)[2], c("release_area"))
  checkmate::assert_choice(colnames(data$x)[3], c("group"))
  checkmate::assert_choice(colnames(data$x)[4], c("count"))

  # Check y --------------------------------------------------------------------

  checkmate::assert_matrix(data$y, mode = "integerish", any.missing = FALSE)
  checkmate::assert_matrix(data$y, ncols = 6)
  checkmate::assert_numeric(data$y, lower = 0, finite = TRUE)
  checkmate::assert_choice(colnames(data$y)[1], c("release_step"))
  checkmate::assert_choice(colnames(data$y)[2], c("release_area"))
  checkmate::assert_choice(colnames(data$y)[3], c("group"))
  checkmate::assert_choice(colnames(data$y)[4], c("recover_step"))
  checkmate::assert_choice(colnames(data$y)[5], c("recover_area"))
  checkmate::assert_choice(colnames(data$y)[6], c("count"))

  # Check z --------------------------------------------------------------------

  checkmate::assert_matrix(data$z, mode = "integerish", any.missing = FALSE)
  checkmate::assert_numeric(data$z, lower = 0, upper = 1)

  # Check l --------------------------------------------------------------------

  checkmate::assert_matrix(data$l, mode = "double", null.ok = TRUE)
  checkmate::assert_matrix(data$l, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_numeric(data$l, lower = 0, upper = 1, null.ok = TRUE)

  # Check w --------------------------------------------------------------------

  checkmate::assert_matrix(data$w, mode = "double", null.ok = TRUE)
  checkmate::assert_matrix(data$w, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_numeric(data$w, lower = 0, upper = 12, null.ok = TRUE)

  # Check f --------------------------------------------------------------------

  checkmate::assert_matrix(data$f, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(data$f, lower = 0, finite = TRUE, null.ok = TRUE)

  # Check m --------------------------------------------------------------------

  checkmate::assert_numeric(data$m, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_numeric(data$m,  finite = TRUE, null.ok = TRUE)

  # Check h --------------------------------------------------------------------

  checkmate::assert_numeric(data$h, lower = 0, len = 1, finite = TRUE)

  # Check u --------------------------------------------------------------------

  checkmate::assert_numeric(data$u, lower = 0, upper = 1, len = 1)

  # Define index limits --------------------------------------------------------

  nt <- max(c(data$x[, "release_step"], data$y[, "recover_step"])) + 1L
  na <- max(c(data$x[, "release_area"], data$y[, "recover_area"])) + 1L
  ng <- max(c(data$x[, "group"], data$y[, "group"])) + 1L # Indexed from zero
  np <- as.integer(sum(data$z))

  # Check dimensions -----------------------------------------------------------

  # z
  checkmate::assert_true(nrow(data$z) == na)
  checkmate::assert_true(ncol(data$z) == na)
  # l
  if (!is.null(data$l)) {
    checkmate::assert_true(nrow(data$l) > 0)
    checkmate::assert_true(ncol(data$l) == 1 || ncol(data$l) == na)
  }
  # w
  if (!is.null(data$w)) {
    checkmate::assert_true(nrow(data$w) > 0)
    checkmate::assert_true(ncol(data$w) == 1 || ncol(data$w) == na)
  }

  # Check values ---------------------------------------------------------------

  # z
  checkmate::assert_true(all(diag(data$z) == 0L))
  # w
  if (!is.null(data$w)) checkmate::assert_true(all(colMeans(data$w) == 1))

  # Check not f and m null -----------------------------------------------------

  if (is.null(data$f) && is.null(data$m)) {
    stop("matrix f or scalar m must be defined in the data")
  }
}

#' Check Parameters Argument for mmmFit()
#'
#' @param parameters [list()]
#'
#' @return NULL
#'
check_parameters <- function (parameters) {

  # Check parameters -----------------------------------------------------------

  checkmate::assert_list(parameters, null.ok = TRUE)

  # Check p --------------------------------------------------------------------

  checkmate::assert_array(parameters$p, mode = "double", null.ok = TRUE)
  checkmate::assert_array(parameters$p, d = 3L, null.ok = TRUE)

  # Check f --------------------------------------------------------------------

  checkmate::assert_matrix(parameters$f, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(parameters$f, lower = 0, null.ok = TRUE)
  checkmate::assert_numeric(parameters$f, finite = TRUE, null.ok = TRUE)

  # Check m --------------------------------------------------------------------

  checkmate::assert_numeric(parameters$m, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_numeric(parameters$m,  finite = TRUE, null.ok = TRUE)

  # Check b --------------------------------------------------------------------

  checkmate::assert_numeric(parameters$b, lower = 0, null.ok = TRUE)
  checkmate::assert_numeric(parameters$b, finite = TRUE, null.ok = TRUE)

  # Check k --------------------------------------------------------------------

  checkmate::assert_numeric(parameters$k, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_numeric(parameters$k,  finite = TRUE, null.ok = TRUE)

}

#' Check Settings Argument for mmmFit()
#'
#' @param settings [list()]
#'
#' @return NULL
#'
check_settings <- function (settings) {

  # Check settings -------------------------------------------------------------

  checkmate::assert_list(settings, null.ok = FALSE)

  # Check error family ---------------------------------------------------------

  checkmate::assert_integerish(settings$error_family, len = 1)
  checkmate::assert_integerish(settings$error_family, lower = 0, upper = 2)
  checkmate::assert_integerish(settings$error_family, any.missing = FALSE)

  # Check max liberty ----------------------------------------------------------

  checkmate::assert_integerish(settings$max_liberty, lower = 0, len = 1)
  checkmate::assert_integerish(settings$max_liberty, any.missing = FALSE)

  # Check time varying ---------------------------------------------------------

  checkmate::assert_integerish(settings$time_varying, len = 1)
  checkmate::assert_integerish(settings$time_varying, lower = 0, upper = 1)
  checkmate::assert_integerish(settings$time_varying, any.missing = FALSE)

  # Check time process ---------------------------------------------------------

  checkmate::assert_integerish(settings$time_process, len = 1)
  checkmate::assert_integerish(settings$time_process, lower = 0, upper = 1)
  checkmate::assert_integerish(settings$time_process, any.missing = FALSE)

  # Check cycle length ---------------------------------------------------------

  checkmate::assert_integerish(settings$cycle_length, len = 1)
  checkmate::assert_integerish(settings$cycle_length, lower = 0)
  checkmate::assert_integerish(settings$cycle_length, any.missing = FALSE)

  # Check block length ---------------------------------------------------------

  checkmate::assert_integerish(settings$block_length, len = 1)
  checkmate::assert_integerish(settings$block_length, lower = 0)
  checkmate::assert_integerish(settings$block_length, any.missing = FALSE)

  # Check results step ---------------------------------------------------------

  checkmate::assert_integerish(settings$results_step, len = 1)
  checkmate::assert_integerish(settings$results_step, lower = 1)
  checkmate::assert_integerish(settings$results_step, any.missing = FALSE)

  # Check nlminb loops ---------------------------------------------------------

  checkmate::assert_integerish(settings$nlminb_loops, len = 1)
  checkmate::assert_integerish(settings$nlminb_loops, lower = 0)
  checkmate::assert_integerish(settings$nlminb_loops, any.missing = FALSE)

  # Check newton iterations ----------------------------------------------------

  checkmate::assert_integerish(settings$newton_iters, len = 1)
  checkmate::assert_integerish(settings$newton_iters, lower = 0)
  checkmate::assert_integerish(settings$newton_iters, any.missing = FALSE)

  # Check OpenMP cores ---------------------------------------------------------

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
get_convergence_diagnostics <- function(sd_report) {
  final_grads <- sd_report$gradient.fixed
  bad_eig <- FALSE
  if (!is.null(sd_report$pdHess)) {
    if (!sd_report$pdHess) {
      warning("The model may not have converged: ",
              "non-positive-definite Hessian matrix.", call. = FALSE)
    } else {
      eigval <- try(1 / eigen(sd_report$cov.fixed)$values, silent = TRUE)
      if (methods::is(eigval, "try-error") || (min(eigval) < .Machine$double.eps * 10)) {
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
