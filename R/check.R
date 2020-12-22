#' Check Tag Release Data
#'
#' @param x [matrix()][mmmTags()] Tag release data
#'
#' @return \code{NULL}
#'
check_x <- function (x) {
  # Check x
  checkmate::assert_matrix(x, mode = "integerish", any.missing = FALSE)
  checkmate::assert_matrix(x, ncols = 4)
  checkmate::assert_numeric(x, lower = 0, finite = TRUE)
  checkmate::assert_choice(colnames(x)[1], c("release_step"))
  checkmate::assert_choice(colnames(x)[2], c("release_area"))
  checkmate::assert_choice(colnames(x)[3], c("group"))
  checkmate::assert_choice(colnames(x)[4], c("count"))
}

#' Check Tag Recovery Data
#'
#' @param y [matrix()][mmmTags()] Tag recovery data
#'
#' @return \code{NULL}
#'
check_y <- function (y) {
  # Check y
  checkmate::assert_matrix(y, mode = "integerish", any.missing = FALSE)
  checkmate::assert_matrix(y, ncols = 6)
  checkmate::assert_numeric(y, lower = 0, finite = TRUE)
  checkmate::assert_choice(colnames(y)[1], c("release_step"))
  checkmate::assert_choice(colnames(y)[2], c("release_area"))
  checkmate::assert_choice(colnames(y)[3], c("group"))
  checkmate::assert_choice(colnames(y)[4], c("recover_step"))
  checkmate::assert_choice(colnames(y)[5], c("recover_area"))
  checkmate::assert_choice(colnames(y)[6], c("count"))
}

#' Check Index Matrix
#'
#' @param z [matrix()][mmmIndex()] Square binary index matrix
#'
#' @return \code{NULL}
#'
check_z <- function (z) {
  # Check z
  checkmate::assert_matrix(z, mode = "integerish", any.missing = FALSE)
  checkmate::assert_numeric(z, lower = 0, upper = 1)
}

#' Check Instantaneous Tag Loss Rate
#'
#' @param h [numeric()] Instantaneous tag loss rate
#'
#' @return \code{NULL}
#'
check_h <- function (h) {
  # Check h
  checkmate::assert_numeric(h, lower = 0, len = 1, finite = TRUE)
}

#' Check Initial Tag Loss Rate (Proportion)
#'
#' @param u [numeric()] Initial tag loss rate (proportion)
#'
#' @return \code{NULL}
#'
check_u <- function (u) {
  # Check u
  checkmate::assert_numeric(u, lower = 0, upper = 1, len = 1)
}

#' Check Tag Reporting Rate Data
#'
#' @param z [array()] Tag reporting rate data
#' @param null_ok [logical()] May l be \code{NULL}?
#'
#' @return \code{NULL}
#'
check_l <- function (l, null_ok = TRUE) {
  # Check l
  checkmate::assert_array(l, mode = "double", d = 3, null.ok = null_ok)
  checkmate::assert_array(l, any.missing = FALSE, null.ok = null_ok)
  checkmate::assert_numeric(l, lower = 0, upper = 1, null.ok = null_ok)
}

#' Check Fishing Mortality Rate
#'
#' @param f [array()] Fishing mortality rate
#' @param null_ok [logical()] May f be \code{NULL}?
#'
#' @return \code{NULL}
#'
check_f <- function (f, null_ok = TRUE) {
  # Check f
  checkmate::assert_array(f, mode = "double", d = 3, null.ok = null_ok)
  checkmate::assert_numeric(f, lower = 0, finite = TRUE, null.ok = null_ok)
}

#' Check Natural Mortality Rate
#'
#' @param m [numeric()] Natural mortality rate
#' @param null_ok [logical()] May m be \code{NULL}?
#'
#' @return \code{NULL}
#'
check_m <- function (m, null_ok = TRUE) {
  # Check m
  checkmate::assert_numeric(m, lower = 0, len = 1, null.ok = null_ok)
  checkmate::assert_numeric(m,  finite = TRUE, null.ok = null_ok)
}

#' Check Movement Parameters
#'
#' @param p [array()] Movement parameters
#' @param null_ok [logical()] May p be \code{NULL}?
#'
#' @return \code{NULL}
#'
check_p <- function (p, null_ok = TRUE) {
  # Check p
  checkmate::assert_array(p, mode = "double", null.ok = null_ok)
  checkmate::assert_array(p, d = 3L, null.ok = null_ok)
}

#' Negative Binomial Dispersioin Parameter
#'
#' @param k [numeric()] Negative Binomial Dispersion. The variance is given
#'   by var = mu + mu / k (NB1) or var = mu + mu^2 / k (NB2).
#' @param null_ok [logical()] May k be \code{NULL}?
#'
#' @return \code{NULL}
#'
check_k <- function (k, null_ok = TRUE) {
  # Check k
  checkmate::assert_numeric(k, lower = 0, len = 1, null.ok = null_ok)
  checkmate::assert_numeric(k,  finite = TRUE, null.ok = null_ok)
}



#' Check Data Argument for mmmFit()
#'
#' @param data [list()]
#' @param exclude [character()] vector of element names in \code{data}
#'   to exclude from check.
#'
#' @return \code{NULL}
#'
check_data <- function (data = NULL, exclude = NULL) {

  # Check data -----------------------------------------------------------------

  checkmate::assert_list(data)

  # Check required elements ----------------------------------------------------

  if (!is.element("x", exclude)) check_x(data$x)
  if (!is.element("y", exclude)) check_y(data$y)
  if (!is.element("z", exclude)) check_z(data$z)
  if (!is.element("h", exclude)) check_h(data$h)
  if (!is.element("u", exclude)) check_u(data$u)

  # Check optional elements ----------------------------------------------------

  if (!is.element("l", exclude)) check_l(data$l)
  if (!is.element("f", exclude)) check_f(data$f)
  if (!is.element("m", exclude)) check_m(data$m)

  # Define index limits --------------------------------------------------------

  # TODO: Continue from here
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
#' @param exclude [character()] vector of element names in \code{parameters}
#'   to exclude from check.
#'
#' @return \code{NULL}
#'
check_parameters <- function (parameters, exclude = NULL) {

  # Check parameters -----------------------------------------------------------

  checkmate::assert_list(parameters, null.ok = TRUE)

  # Check required elements ----------------------------------------------------

  if (!is.element("p", exclude)) check_p(parameters$p)

  # Check optional elements ----------------------------------------------------

  if (!is.element("f", exclude)) check_f(parameters$f)
  if (!is.element("m", exclude)) check_m(parameters$m)
  if (!is.element("k", exclude)) check_k(parameters$k)
}

#' Check Settings Argument for mmmFit()
#'
#' @param settings [list()]
#'
#' @return \code{NULL}
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
