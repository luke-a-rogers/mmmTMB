#' Get Convergence Diagnostics from TMB Model Fit
#'
#' @description Get convergence diagnostics from a TMB model fit. This function is from the
#' 'sdmTMB' package at github.com/pbs-assess/sdmTMB/blob/master/R/fit.R
#'
#' @param sd_report A TMB::sdreport() object
#'
#' @return A list
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
