## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib mmmTMB, .registration = TRUE
## usethis namespace: end
NULL

#' Fit a Markov Movement Model
#'
#' @description Estimate movement rates from single-mark single-recapture
#' tagging data.
#'
#' @param data [list()] Tag and other data. See Details.
#' @param parameters [list()] Initial parameter values. See Details.
#' @param map [list()] Values to override \code{map} defaults.
#'   See [TMB::MakeADFun()].
#' @param settings [list()] Values that help define the model.
#'   See [mmmSet()].
#' @param control [list()] Optimization control options to pass to
#'   [stats::nlminb()]. See [mmmControl()].
#'
#' @details The [list()] argument \code{data} must contain:
#' \itemize{
#'   \item \code{x} An integer matrix of tag release counts formatted
#'     as described in [mmmTags()].
#'   \item \code{y} An integer matrix of tag recovery counts formatted
#'     as described in [mmmTags()].
#'   \item \code{z} A square binary index matrix that indicates allowed
#'     direct movement between areas (from rows to columns) from one
#'     time step to the next. Ones off the diagonal indicate allowed
#'     direct movement. Zeros off the diagonal represent disallowed
#'     direct movement. Numbers on the diagonal are ignored because
#'     self-movement is always allowed. See [mmmIndex()].
#'   \item \code{h} The (positive) scalar instantaneous tag loss rate.
#'   \item \code{u} The scalar proportion of tags lost during release.
#' }
#' The list argument \code{data} may optionally contain:
#' \itemize{
#'   \item \code{l} A matrix of tag reporting rates (proportions).
#'     Defaults to one if omitted (full reporting). See [mmmRates()].
#'   \item \code{w} A matrix of fishing mortality rate weights. Useful
#'     when the tag time step is finer than the fishing rate step.
#'     Defaults to one (equal fishing rates across tag time steps
#'     within a fishing rate time step). See [mmmWeights()].
#'   \item \code{f} A matrix of fishing mortality rates. Estimated as
#'     parameter(s) if omitted. See \code{parameters} and \code{settings}.
#'     Also see [mmmRates()].
#'   \item \code{m} The scalar instantaneous natural mortality rate.
#'     Estimated if omitted. One or both of \code{f} and \code{m} must be
#'     provided. See \code{parameters}.
#' }
#' The list argument \code{parameters} may optionally contain the initial
#' values:
#' \itemize{
#'   \item \code{p} An array of movement parameters. See TBD.
#'   \item \code{f} A matrix of fishing mortality rates. See \code{settings}.
#'   \item \code{m} The scalar instantaneous natural mortality rate.
#'   \item \code{b} The vector fishing mortality bias by group.
#'   \item \code{k} The scalar negative binomial dispersion.
#' }
#'
#' @return A list of class \code{mmmFit}.
#' @export
#'
#' @author Luke A. Rogers
#'
#' @examples
#'
mmmFit <- function (data,
                    parameters = NULL,
                    map = NULL,
                    settings = mmmSet(),
                    control = mmmControl()) {

  # Start the clock ------------------------------------------------------------

  tictoc::tic("mmmFit")

  # Check arguments ------------------------------------------------------------

  check_data(data)
  check_parameters(parameters)
  check_settings(settings)

  # Assign data ----------------------------------------------------------------

  x <- data$x # Tag releases
  y <- data$y # Tag recoveries
  z <- data$z # Movement index
  h <- data$h # Instantaneous tag loss rate
  u <- data$u # Initial tag loss proportion

  # Assign optional data -------------------------------------------------------

  if (!is.null(data$l)) l <- data$l else l <- matrix(1, nrow = 1L, ncol = 1L)
  if (!is.null(data$w)) w <- data$w else w <- matrix(1, nrow = 1L, ncol = 1L)

  # Assign data / parameters ---------------------------------------------------

  if (!is.null(data$f)) f <- data$f else f <- parameters$f # May still be NULL
  if (!is.null(data$m)) m <- data$m else m <- parameters$m # May still be NULL

  # Assign parameters ----------------------------------------------------------

  if (!is.null(parameters$p)) p <- parameters$p else p <- NULL # Movement
  if (!is.null(parameters$b)) b <- parameters$b else b <- NULL # Fishing bias
  if (!is.null(parameters$k)) k <- parameters$k else k <- NULL # NB Dispersion

  # Assign settings ------------------------------------------------------------

  error_family <- settings$error_family
  min_liberty  <- settings$min_liberty
  max_liberty  <- settings$max_liberty
  time_varying <- settings$time_varying
  time_process <- settings$time_process
  cycle_length <- settings$cycle_length
  block_length <- settings$block_length
  results_step <- settings$results_step
  nlminb_loops <- settings$nlminb_loops
  newton_iters <- settings$newton_iters
  openmp_cores <- settings$openmp_cores

  # Set index limits -----------------------------------------------------------

  nt <- max(c(x[, "release_step"], y[, "recover_step"])) + 1L # From zero
  na <- max(c(x[, "release_area"], y[, "recover_area"])) + 1L # From zero
  ng <- max(c(x[, "group"], y[, "group"])) + 1L # Indexed from zero for TMB
  np <- as.integer(sum(z))

  # Create tag reporting rate indexes ------------------------------------------

  nlt <- nrow(l)
  nla <- ncol(l)
  vlt <- rep(c(seq_len(nlt) - 1L), ceiling(nt / nlt))[seq_len(nt)]
  vla <- rep(c(seq_len(nla) - 1L), ceiling(na / nla))[seq_len(na)]

  # Create fishing rate weighting indexes --------------------------------------

  nwt <- nrow(w)
  nwa <- ncol(w)
  vwt <- rep(c(seq_len(nwt) - 1L), ceiling(nt / nwt))[seq_len(nt)]
  vwa <- rep(c(seq_len(nwa) - 1L), ceiling(na / nwa))[seq_len(na)]

  # Create fishing rate indexes ------------------------------------------------

  if (!is.null(f)) nft <- nrow(f) else nft <- 1L
  if (!is.null(f)) nfa <- ncol(f) else nfa <- 1L
  vft <- rep(seq_len(nft) - 1L, each = ceiling(nt / nft))[seq_len(nt)]
  vfa <- rep(seq_len(nfa) - 1L, each = ceiling(na / nfa))[seq_len(na)]
  if (nt %% nft) cat("caution: nft does not divide nt evenly\n")
  if (!is.element(nfa, c(1, na))) stop("nfa must equal 1 or na\n")

  # Create movement parameter indexes ------------------------------------------

  if (time_varying) {
    if (!block_length && !cycle_length) {
      npt <- nt
      vpt <- c(seq_len(nt) - 1L) # Index from zero for TMB
    } else if (block_length && !cycle_length) {
      npt <- ceiling(nt / block_length)
      vpt <- rep(seq_len(npt) - 1L, each = block_length)[seq_len(nt)]
    } else if (cycle_length && !block_length) {
      npt <- cycle_length
      vpt <- rep(seq_len(npt) - 1L, ceiling(nt / npt))[seq_len(nt)]
    } else {
      npt <- cycle_length
      vpt_cycle <- rep(seq_len(npt) - 1L, ceiling(nt / npt))[seq_len(nt)]
      vpt <- rep(vpt_cycle, each = block_length)[seq_len(nt)]
    }
  } else {
    npt <- 1L
    vpt <- rep(0L, nt) # Index from zero for TMB
  }

  # Estimate parameter? --------------------------------------------------------

  estimate_p <- TRUE
  estimate_f <- ifelse(is.null(data$f), TRUE, FALSE)
  estimate_m <- ifelse(is.null(data$m), TRUE, FALSE)
  estimate_b <- ifelse(!estimate_f | ng > 1, TRUE, FALSE)
  estimate_k <- ifelse(error_family, TRUE, FALSE)

  # Assign default initial parameter values ------------------------------------

  if (is.null(p)) p <- array(0, dim = c(npt, np, ng)) # Movement parameters
  if (is.null(f)) f <- array(0.10, dim = c(nft, nfa)) # Fishing rates
  if (is.null(m)) m <- 0.10 # Natural mortality rate
  if (is.null(b)) b <- rep(1, ng) # Fishing bias
  if (is.null(k)) k <- 1L # Negative binomial dispersion parameter

  # Confirm parameter dimensions -----------------------------------------------

  if (!all(dim(p) == c(npt, np, ng))) stop("dim(p) must equal c(npt, np, ng)\n")
  if (!(ncol(f) == 1 || ncol(f) == na)) stop("ncol(f) must equal one or na\n")
  if (nrow(f) < 1) stop("nrow(f) must be greater than zero\n")
  if (length(m) != 1) stop("length(m) must be one")
  if (length(b) != ng) stop("length(b) must equal ng\n")
  if (length(k) != 1) stop("length(k) must be one")

  # Create tmb data ------------------------------------------------------------

  cat("\ncreating tmb_data \n")
  tmb_data <- list(
    tx = t(x), # Tag release matrix
    ty = t(y), # Tag recovery matrix
    tz = t(z), # Movement index matrix
    tl = t(l), # Tag reporting rate matrix
    tw = t(w), # Fish rate weighting if tag step finer than fish rate step
    h = h,     # instantaneous tag loss rate
    d = (1L - u), # Proportion of tags still attached following release
    error_family = error_family,
    time_varying = time_varying,
    time_process = time_process,
    min_liberty  = min_liberty,
    max_liberty  = max_liberty,
    # Index limits
    np = np, # Index limits: number of parameters
    nt = nt, # Index limits: number of time steps
    na = na, # Index limits: number of areas
    ng = ng, # Index limits: number of groups
    npt = npt, # Secondary index limits: number of parameter time steps
    nft = nft, # Secondary index limits: number of fishing rate time steps
    nfa = nfa, # Secondary index limits: number of fishing rate areas
    nlt = nlt, # Secondary index limits: number of reporting rate time steps
    nla = nla, # Secondary index limits: number of reporting rate areas
    nwt = nwt, # Secondary index limits: number of fish rate weighting steps
    nwa = nwa, # Secondary index limits: number of fish rate weighting areas
    # Index vectors
    vpt = vpt, # Index vector: time step for movement parameters
    vft = vft, # Index vector: time step for fishing rate
    vfa = vfa, # Index vector: area for fishing rate
    vlt = vlt, # Index vector: time step for tag reporting rate
    vla = vla, # Index vector: area for tag reporting rate
    vwt = vwt, # Index vector: time step for fishing rate weighting
    vwa = vwa  # Index vector: area for fishing rate weighting
  )

  # Create tmb parameters ------------------------------------------------------

  cat("creating tmb_parameters \n")
  tmb_parameters <- list(
    tp = aperm(p, c(2,1,3)), # Movement parameters
    logit_exp_neg_tf = logit(exp(-(t(f)))), # Fishing mortality
    logit_exp_neg_m = logit(exp(-(m))), # Natural mortality
    log_b = log(b), # Fishing bias
    log_k = log(k) # Negative binomial dispersion
  )

  # Create tmb map -------------------------------------------------------------

  cat("creating tmb_map \n")
  tmb_map <- list()
  # Default
  if (!is.null(data$f)) tmb_map$logit_exp_neg_tf <- factor(rep(NA, nfa * nft))
  if (!is.null(data$m)) tmb_map$logit_exp_neg_m <- as.factor(NA)
  if (!estimate_b) tmb_map$log_b <- as.factor(rep(NA, ng))
  if (error_family == 0) tmb_map$log_k <- as.factor(NA)
  # User defined
  if (!is.null(map$p)) tmb_map$tp <- aperm(map$p, c(2, 1, 3))
  if (!is.null(map$f)) tmb_map$logit_exp_neg_tf <- t(map$f)
  if (!is.null(map$m)) tmb_map$logit_exp_neg_m <- map$m
  if (!is.null(map$b)) tmb_map$log_b <- map$b

  # Set OpenMP cores -----------------------------------------------------------

  TMB::openmp(n = settings$openmp_cores)

  # Create ADFun object --------------------------------------------------------

  tictoc::tic("adfun")
  cat("creating adfun \n")
  adfun <- TMB::MakeADFun(data = tmb_data,
                          parameters = tmb_parameters,
                          map = tmb_map,
                          # random = tmb_random,
                          DLL = "mmmTMB")
  tictoc::toc()

  # Optimize the objective function --------------------------------------------

  tictoc::tic("nlminb")
  cat("\noptimizing objective\n")
  cat("\nmodel mgc \n")

  # Optimize
  model <- stats::nlminb(
    start = adfun$par,
    objective = adfun$fn,
    gradient = adfun$gr,
    control = control)

  # Iterate optimization
  if (nlminb_loops > 0) {
    cat("\nrunning extra nlminb loops \n")
    for (i in seq_len(nlminb_loops)) {
      cat(paste0("running nlminb loop #", i, "\n"))
      tmp <- model[c("iterations", "evaluations")]
      model <- stats::nlminb(
        start = model$par,
        objective = adfun$fn,
        gradient = adfun$gr,
        control = control)
      model[["iterations"]] <- model[["iterations"]] + tmp[["iterations"]]
      model[["evaluations"]] <- model[["evaluations"]] + tmp[["evaluations"]]
    }
  }

  # Perform Newton iterations
  if (newton_iters > 0) {
    cat("\nrunning newton iters \n")
    for (i in seq_len(newton_iters)) {
      cat(paste0("running newton iter #", i, "\n"))
      g <- as.numeric(adfun$gr(model$par))
      h <- optimHess(model$par, fn = adfun$fn, gr = adfun$gr)
      model$par <- model$par - solve(h, g)
      model$objective <- adfun$fn(model$par)
    }
  }
  tictoc::toc()

  # Create sd_report -----------------------------------------------------------

  tictoc::tic("sd_report")
  cat("\ncreating sd_report")
  cat("\nsd_report mgc \n")
  sd_report <- TMB::sdreport(adfun)
  convergence <- get_convergence_diagnostics(sd_report)
  convergence$mgc <- max(abs(convergence$final_grads))
  cat("positive definite hessian:", sd_report$pdHess, "\n")
  cat("mgc:", convergence$mgc, "\n")
  tictoc::toc()

  # Compute movement rate results ----------------------------------------------

  # Extract values
  vtp_fit <- subset_by_name(model$par, "tp")
  mtp_cov <- subset_by_name(sd_report$cov.fixed, "tp")
  # Compute
  movement_list <- create_movement_results(
    vtp = vtp_fit,
    mtp = mtp_cov,
    np = np,
    npt = npt,
    ng = ng,
    z = z,
    pow = results_step,
    draws = 1000
  )

  # Compute natural mortality results ------------------------------------------

  # Extract values
  logit_exp_neg_m_fit <- subset_by_name(model$par, "logit_exp_neg_m")
  logit_exp_neg_m_se <- subset_by_name(sd_report$cov.fixed, "logit_exp_neg_m")
  # Compute
  mortality_list <- create_mortality_results(
    logit_exp_neg_m = logit_exp_neg_m_fit,
    logit_exp_neg_m_se = logit_exp_neg_m_se,
    results_step = results_step,
    estimate = estimate_m
  )

  # Compute fishing mortality rate results -------------------------------------

  # Define source of covariance matrix
  cov_source <- sd_report$cov.fixed
  # Extract values
  vlogit_exp_neg_tf_fit <- subset_by_name(model$par, "logit_exp_neg_tf")
  mlogit_exp_neg_tf_cov <- subset_by_name(cov_source, "logit_exp_neg_tf")
  # Compute
  fishing_rate_list <- create_fishing_results(
    vlogit_exp_neg_tf = vlogit_exp_neg_tf_fit,
    mlogit_exp_neg_tf_cov = mlogit_exp_neg_tf_cov,
    results_step = results_step,
    nft = nft,
    nfa = nfa,
    estimate = estimate_f
  )

  # Compute fishing mortality bias results -------------------------------------

  # Extract values
  log_b <- subset_by_name(model$par, "log_b")
  log_b_cov <- subset_by_name(sd_report$cov.fixed, "log_b")
  # Compute
  fishing_bias_list <- create_bias_results(
    log_b = log_b,
    log_b_cov = log_b_cov,
    estimate = estimate_b
  )

  # Compute dispersion results -------------------------------------------------

  # Extract values
  log_k_fit <- subset_by_name(model$par, "log_k")
  log_k_se <- subset_by_name(sd_report$cov.fixed, "log_k")
  # Compute
  dispersion_list <- create_dispersion_results(
    log_k = log_k_fit,
    log_k_se = log_k_se,
    estimate = estimate_k
  )

  # Assemble results -----------------------------------------------------------

  results <- list(
    movement_rate = movement_list$movement_rate,
    natural_mortality = mortality_list$natural_mortality,
    fishing_rate = fishing_rate_list$fishing_rate,
    fishing_bias = fishing_bias_list$fishing_bias,
    dispersion = dispersion_list$dispersion
  )

  # Assemble inputs ------------------------------------------------------------

  inputs <- list(
    data       = data,
    parameters = parameters,
    settings   = settings,
    # random     = random,
    map        = map,
    control    = control
  )

  # Assemble parameters --------------------------------------------------------

  parameters <- list(
    p = movement_list$p_fit,
    f = fishing_rate_list$f_fit,
    m = mortality_list$m_fit,
    b = fishing_bias_list$b_fit,
    k = dispersion_list$k_fit
  )

  # Stop the clock -------------------------------------------------------------

  tictoc::toc()

  # Return mmmFit object -------------------------------------------------------

  cat("\nreturning mmmFit object\n")
  structure(list(
    inputs      = inputs,
    adfun       = adfun,
    model       = model,
    results     = results,
    sd_report   = sd_report,
    parameters  = parameters,
    convergence = convergence),
    class       = c("mmmFit"))
}

#' Settings for \code{mmmFit()}
#'
#' @param error_family Integer. \code{0} = Poisson; \code{1} = NB1;
#'   \code{2} = NB2.
#' @param max_liberty Integer. Max time steps at liberty before
#'   recapture.
#' @param time_varying Integer. \code{0} = None; \code{1} = Time varying
#' movement rates.
#' @param time_process Integer. \code{0} = None; \code{1} = Negative
#' binomial (NB1).
#' @param cycle_length Integer. Cycle length or \code{0} for no cycle.
#' @param block_length Integer. Block length for movement rate.
#' @param results_step Integer. Results time step as a multiple of the data
#'   time step.
#' @param nlminb_loops Integer. Number of times to run [stats::nlminb()]
#'   optimization.
#' @param newton_steps Integer. Number of Newton optimization steps.
#' @param openmp_cores Integer. Number of parallel cores for TMB.
#'
#' @return A list of settings for \code{mmmFit()}
#' @export
#'
#' @examples
#'
mmmSet <- function (error_family = c(nb1 = 1, poisson = 0),
                    min_liberty  = 1,
                    max_liberty  = 1000,
                    time_varying = 0,
                    time_process = c(none = 0, rw = 1),
                    cycle_length = 0,
                    block_length = 0,
                    results_step = 1,
                    nlminb_loops = 5,
                    newton_iters = 5,
                    openmp_cores = 1) {

  #--------------- Return a list of settings ----------------------------------#

  structure(list(
    error_family = error_family[1],
    min_liberty  = min_liberty,
    max_liberty  = max_liberty,
    time_varying = time_varying,
    time_process = time_process[1],
    cycle_length = cycle_length,
    block_length = block_length,
    results_step = results_step,
    nlminb_loops = nlminb_loops,
    newton_iters = newton_iters,
    openmp_cores = openmp_cores),
    class = "mmmSet")
}


#' Optimization control options
#'
#' Any arguments to pass to \code{stats::nlminb()}.
#'
#' @param eval.max Maximum number of evaluations of the objective
#' function allowed.
#' @param iter.max Maximum number of iterations allowed.
#' @param ... Anything else. See the 'Control parameters' section of
#'   \code{stats::nlminb()}.
#'
#' @export
#'
mmmControl <- function (eval.max = 1e4, iter.max = 1e4, ...) {
  list(eval.max = eval.max, iter.max = iter.max, ...)
}
