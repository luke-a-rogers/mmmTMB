## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib mmmTMB, .registration = TRUE
## usethis namespace: end
NULL


#' Fit a Markov Movement Model to Mark-Recapture Data
#'
#' @description Estimate movement probabilities among areas by fitting a Markov
#' movement model to mark release and recovery data.
#'
#' Array dimensions are abbreviated for space and are given by:
#' \itemize{
#'   \item{\code{nt}: Count of time steps in the data}
#'   \item{\code{na}: Count of areas in the data}
#'   \item{\code{ng}: Count of group categories}
#' }
#'
#' @usage mmmTMB(data_list, ...)
#'
#' @param data_list [list()] A named list of data objects (see details)
#' @param structure_list [list()] A named list specifying model structure (see
#' details)
#' @param optimizer_list [list()] A named list of optimizer arguments (see
#' details)
#' @param parameter_list [list()] An optional named list of initial values
#' (see details)
#'
#' @details
#' The argument \code{data_list} must be a named \code{list()} of objects:
#' \itemize{
#'   \item{\code{released_3d} an \code{array()} with \code{dim = c(nt, na, ng)}}
#'   \item{\code{recovered_5d} an \code{array()} with
#'     \code{dim = c(nt, nt, na, na, ng)}}
#'   \item{\code{capture_rate_2d} an \code{array()} with \code{dim = c(na, nt)}}
#'   \item{\code{report_ratio_2d} an \code{array()} with \code{dim = c(na, nt)}}
#'   \item{\code{template_2d} an \code{array()} with \code{dim = c(na, na)}}
#'   \item{\code{tag_loss_rate} a \code{numeric()} with \code{length = 1}}
#'   \item{\code{imm_loss_ratio} a \code{numeric()} with \code{length = 1}}
#' }
#'
#' The argument \code{structure_list} must be a named \code{list()} of objects:
#' \itemize{
#'   \item{\code{recapture_delay} \code{integer()} Minimum delay between
#'   release and recovery}
#'   \item{\code{error_family} \code{integer()} One of 0: Poisson; 1: NB1;
#'   or 2: NB2}
#'   \item{\code{result_units} \code{integer()} Results time step as multiple
#'   of data time step}
#' }
#'
#' The argument \code{optimizer_list} must be a named \code{list()} of objects:
#' \itemize{
#'   \item{\code{newton_steps} \code{integer()} Number of Newton optimization
#'   steps}
#'   \item{\code{nlminb_loops} \code{integer()} Number of times to run
#'   [stats::nlminb()] optimization.}
#'   \item{\code{nlminb_control} \code{list()} see \code{mmmTMBcontrol()}}
#'   \item{\code{openmp_cores} \code{integer()} Number of cores for TMB}
#' }
#'
#'
#' The argument \code{parameter_list} must be \code{NULL} or a named
#' \code{list()} of objects:
#' \itemize{
#'   \item{TBD...}
#' }
#'
#' @return An object of class \code{mmmTMB}
#'
#' @export
#'
#' @examples
#'
mmmTMB <- function (data_list,
                    structure_list,
                    optimizer_list,
                    parameter_list = NULL) {

  #---------------- Fix certain arguments --------------------------------------#

  time_process <- 0

  #---------------- Start the clock -------------------------------------------#

  tictoc::tic("mmmTMB")

  #---------------- Unpack arguments ------------------------------------------#

  released_3d <- data_list$released_3d
  template_2d <- data_list$template_2d
  error_family <- structure_list$error_family
  nlminb_control <- optimizer_list$nlminb_control
  nlminb_loops <- optimizer_list$nlminb_loops
  newton_steps <- optimizer_list$newton_steps
  openmp_cores <- optimizer_list$openmp_cores


  #---------------- Check arguments -------------------------------------------#

  # Assign dimensions
  np <- sum(template_2d)
  nt <- dim(released_3d)[1]
  na <- dim(released_3d)[2]
  ng <- dim(released_3d)[3]

  #---------------- Define the number of cores --------------------------------#

  openmp_cores <- structure_list$openmp_cores
  if (!is.null(openmp_cores)) {TMB::openmp(n = openmp_cores)}

  #---------------- Define the movement index ---------------------------------#

  if (time_process == 0) {
    nv <- 1
    movement_index <- rep(0, nt)
  } # TODO: else if (time_process > 0) {}

  #---------------- Create the data list --------------------------------------#

  cat("creating tmb_data \n")
  tmb_data <- c(data_list,
                structure_list,
                movement_index = list(movement_index))

  #---------------- Create the parameter list ---------------------------------#

  cat("creating tmb_parameters \n")
  if (!is.null(parameter_list)) {
    tmb_parameters <- parameter_list
  } else {
    tmb_parameters <- list(
      movement_parameters_3d = array(-3, dim = c(nv, np, ng)),
      log_capture_bias_2d = array(0, dim = c(na, ng)),
      log_natural_mortality = -3,
      log_dispersion = 0.8
    )
  }

  #---------------- Define parameter structure --------------------------------#

  # Initialize tmb_map
  cat("creating tmb_map \n")
  tmb_map <- list()

  # TODO: Add capture_bias_map argument and conditional
  capture_bias_map <- factor(rep(1, na * ng))
  tmb_map <- c(tmb_map, list(log_capture_bias_2d = capture_bias_map))

  # Map off log_dispersion?
  if (error_family == 0) {
    tmb_map <- c(tmb_map, list(log_dispersion = as.factor(NA)))
  }

  #---------------- Define random effects -------------------------------------#

  # Initialize tmb_random
  tmb_random <- character(0)

  # TODO: Add random effects timevarying process

  #---------------- Create the simulation object ------------------------------#


  #---------------- Create the model object -----------------------------------#
  # TODO: Add previous fit option

  cat("creating tmb_obj \n")
  tmb_obj <- TMB::MakeADFun(data = tmb_data,
                            parameters = tmb_parameters,
                            map = tmb_map,
                            random = tmb_random,
                            DLL = "mmmTMB")

  #---------------- Optimize the objective function ---------------------------#

  cat("\nmodel mgc \n")
  {
    # Start the clock
    tictoc::tic("nlminb")

    # Perform initial optimization
    tmb_opt <- stats::nlminb(
      start = tmb_obj$par,
      objective = tmb_obj$fn,
      gradient = tmb_obj$gr,
      control = nlminb_control)

    # Iterate optimization
    if (nlminb_loops > 0) {
      cat("running extra nlminb loops \n")
      for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
        temp <- tmb_opt[c("iterations", "evaluations")]
        tmb_opt <- stats::nlminb(
          start = tmb_opt$par,
          objective = tmb_obj$fn,
          gradient = tmb_obj$gr,
          control = nlminb_control)
        tmb_opt[["iterations"]] <- tmb_opt[["iterations"]] + temp[["iterations"]]
        tmb_opt[["evaluations"]] <- tmb_opt[["evaluations"]] + temp[["evaluations"]]
      }
    }

    # Perform additional Newton steps
    if (newton_steps > 0) {
      cat("\nrunning newtonsteps \n")
      for (i in seq_len(newton_steps)) {
        g <- as.numeric(tmb_obj$gr(tmb_opt$par))
        h <- optimHess(tmb_opt$par, fn = tmb_obj$fn, gr = tmb_obj$gr)
        tmb_opt$par <- tmb_opt$par - solve(h, g)
        tmb_opt$objective <- tmb_obj$fn(tmb_opt$par)
      }
    }
    # Stop the clock
    tictoc::toc()
  }

  #---------------- Create sd report ------------------------------------------#

  cat("\nsd_report mgc \n")
  sd_report <- TMB::sdreport(tmb_obj)
  conv <- get_convergence_diagnostics(sd_report)
  mgc <- max(abs(conv$final_grads))

  #---------------- Create results --------------------------------------------#\
  # TODO: Compute results and SEs

  # Results list
  results_list <- list(movement_probabilities = NULL,
                       natural_mortality = NULL,
                       capture_bias = NULL,
                       dispersion = NULL)

  #---------------- Stop the clock --------------------------------------------#

  tictoc::toc()

  #---------------- Return an mmmTMB object -----------------------------------#

  structure(list(
    adfun      = tmb_obj, # report(),
    model      = tmb_opt, # convergence, objective,
    data       = data_list, # TODO: optionally place in simulated data
    structure  = structure_list,
    parameters = tmb_parameters,
    results    = results_list,
    map        = tmb_map,
    random     = tmb_random,
    sd_report  = sd_report, # value, sd, cov, par.fixed, cov.fixed, pdHess, gradient.fixed, env
    conv       = conv, # final_grads, bad_eig
    mgc        = mgc),
    class      = "mmmTMB")

}


#' Optimization control options
#'
#' Any arguments to pass to [stats::nlminb()].
#'
#' @param eval.max [numeric(1)] Maximum number of evaluations of the objective
#' function allowed.
#' @param iter.max [numeric(1)] Maximum number of iterations allowed.
#' @param ... Anything else. See the 'Control parameters' section of
#'   [stats::nlminb()].
#'
#' @export
#'
mmmTMBcontrol <- function(eval.max = 1e4, iter.max = 1e4, ...) {
  list(eval.max = eval.max, iter.max = iter.max, ...)
}
