## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib mmmTMB, .registration = TRUE
## usethis namespace: end
NULL

#' Fit a Markovian Movement Model
#'
#' @description Estimate movement rates
#'
#' @param data A list of named data objects. See details.
#' @param parameters A list of initial parameter values. See details.
#' @param settings A list of values that define the model. See details.
#' @param random A character vector of parameters to treat as random effects.
#' @param map A list of values to override defaults.
#' @param nlminb_control A list. See \code{mmmControl()}.
#'
#' @details The following must be included in the \code{data} list: \itemize{
#'   \item \code{T} An array holding tagged release counts. See \code{arrayT()}.
#'   \item \code{R} An array holding recovered counts. See \code{arrayR()}.
#'   \item \code{I} A square binary matrix indicating which direct movement
#'   rates to estimate.} The following can be optionally included in the
#'   \code{data} list and are otherwise excluded from the model \itemize{ \item
#'   \code{lambda} An array of tag reporting rates. } The following can be
#'   optionally included in the \code{data} list to be treated as data, or
#'   omitted to be treated as parameters with default initial values, or
#'   included in the \code{parameters} list along with initial values \itemize{
#'   \item \code{M} An array of instantaneous natural mortality rates. \item
#'   \code{F} An array of instantaneous fishing mortality rates. \item \code{b}
#'   A vector of selectivity values. \item \code{h} A scalar instantaneous tag
#'   loss rate. \item \code{c} A scalar mean initial tag loss rate.} The
#'   following can be included in the \code{settings} list, or omitted to be
#'   assigned default values \itemize{ \item \code{family} An integer specifying
#'   the error family 0: Poisson; 1: NB1; 2: NB2. \item \code{timestep} An
#'   integer specifying the results time step as a multiple of the data time
#'   step. \item \code{nlminb_loops} An integer specifying the number of
#'   \code{nlminb} loops. \item \code{newton_steps} An integer specifying the
#'   number of Newton steps. \item \code{openmp_cores} An integer specifying the
#'   number of \code{OpenMP} parallel cores.}
#'
#' @return A list of class \code{mmmFit} and \code{mmmTMB}.
#' @export
#'
#' @author Luke A. Rogers
#'
#' @examples
#'
mmmFit <- function(data,
                   parameters = NULL,
                   map = NULL,
                   random = NULL,
                   settings = mmmSet(),
                   control = mmmControl()) {

  #---------------- Start the clock -------------------------------------------#

  tictoc::tic("mmmFit")

  #---------------- Check data argument ---------------------------------------#

  cat("checking arguments \n")
  # TODO

  #---------------- Check parameters argument ---------------------------------#
  # TODO

  #---------------- Check settings argument -----------------------------------#

  # Error family
  if (is.null(settings$family)) {
    family <- 1L
  } else {
    stopifnot(
      is.numeric(settings$family),
      length(settings$family) == 1,
      is.element(settings$family, c(0L, 1L, 2L)))
    family <- settings$family
  }
  # Results timestep
  if (is.null(settings$timestep)) {
    timestep <- 1L
  } else {
    stopifnot(
      is.numeric(settings$timestep),
      length(settings$timestep) == 1,
      settings$timestep > 0,
      settings$timestep == as.integer(settings$timestep))
    timestep <- settings$timestep
  }
  # Number nlminb_loops
  if (is.null(settings$nlminb_loops)) {
    nlminb_loops <- 5L
  } else {
    stopifnot(
      is.numeric(settings$nlminb_loops),
      length(settings$nlminb_loops) == 1,
      settings$nlminb_loops >= 0,
      settings$nlminb_loops == as.integer(settings$nlminb_loops))
    nlminb_loops <- settings$nlminb_loops
  }
  # Number newton_steps
  if (is.null(settings$newton_steps)) {
    newton_steps <- 5L
  } else {
    stopifnot(
      is.numeric(settings$newton_steps),
      length(settings$newton_steps) == 1,
      settings$newton_steps >= 0,
      settings$newton_steps == as.integer(settings$newton_steps))
    newton_steps <- settings$newton_steps
  }
  # Number openmp_cores
  if (is.null(settings$openmp_cores)) {
    openmp_cores <- as.integer(parallel::detectCores() / 2)
  } else {
    stopifnot(
      is.numeric(settings$openmp_cores),
      length(settings$openmp_cores) == 1,
      settings$openmp_cores > 0,
      settings$openmp_cores <= parallel::detectCores(),
      settings$openmp_cores == as.integer(settings$openmp_cores))
    openmp_cores <- settings$openmp_cores
  }

  #---------------- Unpack arguments ------------------------------------------#

  cat("unpacking arguments \n")
  # Arguments tags, mT, mR, liberty
  if (is.element("tags", names(data))) {
    # Check
    checkmate::assert_list(data$tags)
    checkmate::assert_true(all(is.element(c("mT", "mR"), names(data$tags))))
    checkmate::assert_matrix(data$tags$mT, mode = "integer", ncols = 4)
    checkmate::assert_matrix(data$tags$mR, mode = "integer", ncols = 6)
    checkmate::assert_class(data$tags, "mmmTags")
    # Assign
    cat("using mT and mR from tags data \n")
    mT <- data$tags$mT
    mR <- data$tags$mR
    cat("using time at liberty from tags data \n")
    liberty <- data$tags$steps_liberty
  } else {
    # Check
    checkmate::assert_matrix(data$mT, mode = "integerish", ncols = 4)
    checkmate::assert_matrix(data$mR, mode = "integerish", ncols = 6)
    checkmate::assert_numeric(settings$liberty, len = 2, null.ok = TRUE)
    # Assign
    mT <- data$mT
    mR <- data$mR
    if (is.null(settings$liberty)) {
      cat("using default time at liberty \n")
      liberty <- c(0, Inf)
    } else {
      liberty <- numeric(2)
      cat("using time at liberty from settings \n")
      liberty[1] <- max(0, settings$liberty[1], na.rm = TRUE)
      liberty[2] <- max(liberty[1], settings$liberty[2], na.rm = TRUE)
    }
  }
  cat(paste0("liberty = c(", liberty[1], ",", liberty[2], ") \n"))
  # Argument mI
  checkmate::assert_matrix(data$mI, mode = "integerish")
  mI <- data$mI
  # Argument mL: tag reporting rates
  if (!is.null(data$mL)) {
    mL <- data$mL
  } else {
    mL <- array(1L, dim = c(nt, na))
  }





  # DONE: mT, mR, and liberty, mI, mL

  # TODO


  # diag(I) <- 0L

  #---------------- Check dimensions ------------------------------------------#

  cat("checking argument dimensions \n")
  # TODO:


  #---------------- Create index vectors --------------------------------------#

  cat("creating index vectors \n")
  # TODO: Use settings input for these

  # Index limits
  np <- sum(mI)
  nt <- max(c(mT$release_step, mR$recover_step)) + 1L # Index starts at zero
  na <- max(c(mT$release_area, mR$recover_area)) + 1L # Index starts at zero
  ng <- max(c(mT$group, mR$group)) + 1L # Index starts at zero
  # Secondary index limits
  if (settings$time_varying) {
    if (!settings$block_length && !settings$cycle_length) {
      npt <- nt
      vpt <- c(seq_len(npt) - 1L) # Index from zero for C++
    } else if (settings$block_length && !settings$cycle_length) {
      npt <- ceiling(nt / settings$block_length)
      vpt <- rep(seq_len(npt) - 1L, each = settings$block_length)[seq_len(nt)]
    } else if (settings$cycle_length && !settings$block_length) {
      npt <- settings$cycle_length
      vpt <- rep(seq_len(npt) - 1L, ceiling(nt / npt))[seq_len(nt)]
    } else {
      npt <- settings$cycle_length
      vpt_cycle <- rep(seq_len(npt) - 1L, ceiling(nt / npt))[seq_len(nt)]
      vpt <- rep(vpt_cycle, each = settings$block_length)[seq_len(nt)]
    }
  } else {
    npt <- 1L
    vpt <- rep(0L, nt) # Index from zero for C++
  }

  nft <-  # Count of fishing rate time steps
  nfa <-  # Count of fishing rate areas
  # Index vectors



  #---------------- Create the tmb_data ---------------------------------------#

  cat("creating tmb_data \n")
  tmb_data <- list(
    mT = mT,
    mR = mR,
    mI = mI,
    mL = mL,
    family = family,
    # Index vectors
    np = np, # Index limits: number of parameters
    nt = nt, # Index limits: number of time steps
    na = na, # Index limits: number of areas
    ng = ng, # Index limits: number of groups
    npt = , # Secondary index limits: number of parameter time steps
    nft = , # Secondary index limits: number of fishing rate time steps
    nfa = , # Secondary index limits: number of fishing rate areas
    vpt = , #
    vft = , #
    vfa = , #
  )

  #---------------- Create parameter values -----------------------------------#

  # Array movement parameters
  if (!is.null(parameters$aP)) {
    aP <- parameters$aP
  } else {
    aP <- array(0, dim = c(npt, np, ng))
  }
  # Matrix log fishing mortality rate
  if (!is.null(data$mF)) {
    mlF <- log(data$mF)
  } else {
    mlF <- array(-3, dim = c(nft, nfa)) # nfg
  }
  # Vector log fishing selectivity and tag reporting bias
  if (!is.null(data$mB)) {
    mlB <- log(data$mB)
  } else {
    mlB <- array(0, dim = c(nba, nbg))
  }
  # Scalar log tag loss rate
  if (!is.null(data$sH)) {
    slH <- log(data$sH)
  } else {
    slH <- -3L
  }
  # Scalar log initial tag loss rate
  if (!is.null(data$sC)) {
    slC <- log(data$sC)
  } else {
    slC <- -2L
  }
  # Scalar log natural mortality
  if (!is.null(data$sM)) {
    slM <- log(data$sM)
  } else {
    slM <- -2L
  }
  # Scalar log negative binomial dispersal
  slD <- 0L

  #---------------- Create the tmb_parameters ---------------------------------#

  cat("creating tmb_parameters \n")
  tmb_parameters <- list(
    aP = aP, # Array: movement parameters
    mlF = mlF, # Matrix: log fishing mortality rates
    mlB = mlB, # Matrix: log fishery selectivity and tag reporting bias
    slH = slH, # Scalar: log tag loss rate
    slC = slC, # Scalar: log initial tag loss rate
    slM = slM, # Scalar: natural mortality
    slD = slD # Scalar: negative binomial dispersion
  )

  #---------------- Create the tmb_map ----------------------------------------#

  cat("creating tmb_map \n")
  tmb_map <- list()
  # Default
  if (!is.null(data$mF)) { tmb_map <- c(tmb_map, mlF = NA) }
  if (!is.null(data$mB)) { tmb_map <- c(tmb_map, mlB = NA) }
  if (!is.null(data$sH)) { tmb_map <- c(tmb_map, slH = NA) }
  if (!is.null(data$sC)) { tmb_map <- c(tmb_map, slC = NA) }
  if (!is.null(data$sM)) { tmb_map <- c(tmb_map, slM = NA) }
  if (family == 0) { tmb_map <- c(tmb_map, slD = NA)}
  # User defined
  if (!is.null(map$mF)) { tmb_map$mlF <- map$mF }
  if (!is.null(map$mB)) { tmb_map$mlB <- map$mB }
  if (!is.null(map$sH)) { tmb_map$slH <- map$sH }
  if (!is.null(map$sC)) { tmb_map$slC <- map$sC }
  if (!is.null(map$sM)) { tmb_map$slM <- map$sM }

  #---------------- Create the tmb_random -------------------------------------#

  cat("creating tmb_random \n")
  if (is.null(random)) {
    tmb_random <- character()
  } else {
    tmb_random <- random
  }

  #---------------- Define the number of cores --------------------------------#

  cat("setting openmp_cores")
  if (!is.null(settings$openmp_cores)) {TMB::openmp(n = settings$openmp_cores)}

  #---------------- Create the ADFun object -----------------------------------#

  tictoc::tic("adfun")
  cat("creating adfun \n")
  adfun <- TMB::MakeADFun(data = tmb_data,
                          parameters = tmb_parameters,
                          map = tmb_map,
                          random = tmb_random,
                          DLL = "mmmTMB")
  tictoc::toc()

  #---------------- Optimize the objective function ---------------------------#

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
    for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
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

  # Perform Newton steps
  if (newton_steps > 0) {
    cat("\nrunning newton steps \n")
    for (i in seq_len(newton_steps)) {
      cat(paste0("running newton step #", i, "\n"))
      g <- as.numeric(adfun$gr(model$par))
      h <- optimHess(model$par, fn = adfun$fn, gr = adfun$gr)
      model$par <- model$par - solve(h, g)
      model$objective <- adfun$fn(model$par)
    }
  }
  tictoc::toc()

  #---------------- Create sd_report ------------------------------------------#

  tictoc::tic("sd_report")
  cat("creating sd_report")
  cat("\nsd_report mgc \n")
  sd_report <- TMB::sdreport(adfun)
  conv_list <- get_convergence_diagnostics(sd_report)
  mgc <- max(abs(conv_list$final_grads))
  tictoc::toc()

  #---------------- Compute goodness of fit -----------------------------------#

  tictoc::tic("goodness")
  cat("computing goodness of fit\n")
  # TODO





  tictoc::toc()

  #---------------- Compute results -------------------------------------------#

  tictoc::tic("results")
  cat("computing results")
  # TODO
  # Use result_steps here





  tictoc::toc()

  #---------------- Stop the clock --------------------------------------------#

  tictoc::toc()

  #---------------- Return an mmmTMB object -----------------------------------#

  cat("\nreturning mmmTMB object\n")
  structure(list(
    data        = data_list,
    parameters  = parameters_list,
    settings    = settings,
    random      = random,
    map         = map_list,
    control     = control,
    initial     = initial_list,
    results     = results_list,
    sd_report   = sd_report,
    convergence = conv_list,
    goodness    = goodness_list,
    adfun       = adfun,
    model       = model,
    mgc         = mgc),
    class       = c("mmmFit", "mmmTMB"))
}




#' Fit a Markov Movement Model to Mark-Recapture Data
#'
#' @description Estimate movement probabilities among areas by fitting a Markov
#' movement model to mark release and recovery data.
#'
#' Array dimensions are abbreviated for space and are given by:
#' \itemize{
#'   \item{\code{nt}: Count of time units in the data}
#'   \item{\code{na}: Count of areas in the data}
#'   \item{\code{ng}: Count of class categories}
#' }
#'
#' @usage mmmTMB(data_list, ...)
#'
#' @param released_3d [array()] Release counts by release time, release area,
#' and class; \code{dim = c(nt, na, ng)}
#' @param recovered_5d [array()] Recovery counts by release time, recovery
#' time, release area, recovery area, and class;
#' \code{dim = c(nt, nt, na, na, ng)}
#' @param capture_rate_2d [array()] Instantaneous capture rate per time unit by
#' time and area; \code{dim = c(nt, na)}
#' @param report_ratio_2d [array()] Reporting ratio by time and area;
#' \code{dim = c(nt, na)}
#' @param tag_loss_rate [numeric()] Instantaneous tag loss rate per time unit,
#' assumed constant across times, areas, and classes
#' @param imm_loss_ratio [numeric()] Ratio of tags lost during release, assumed
#' constant across times, areas, and classes
#' @param template_2d [array()] Template indicating movement parameters to
#' estimate; \code{dim = c(na, na)}. See Details.
#' @param recapture_delay [integer()] Number of time units after release when
#' recovery becomes possible, \code{>= 1}.
#' @param error_family [integer()] One of 0: Poisson; 1: NB1; or 2: NB2
#' @param result_units [integer()] Results time unit as a multiple of the data
#' time unit. See details.
#' @param time_process [integer()] Time process for movement
#' rates. 0: Constant; 1: Time varying
#' @param time_pattern [integer()] Timevarying pattern for movement
#' probabilities. One of 0: Constant; 1: Stepped; 2: Cyclical. See Details.
#' @param pattern_size [integer()] Number of levels in the \code{time_pattern}.
#' @param newton_steps [integer()] Number of Newton optimization steps
#' @param nlminb_loops [integer()] Number of times to run [stats::nlminb()]
#' optimization.
#' @param openmp_cores [integer()] Number of cores for TMB
#' @param simulate_value [logical()] Simulate \code{recovered_5d}? See Details.
#' @param capture_map_2d [array()] Custom \code{TMB} map for
#' \code{log_capture_bias_2d} parameters, \code{dim = c(na, ng)}. See Details.
#' @param structure_list [list()] Input data and model structure constants as
#' a list. See Details.
#' @param parameter_list [list()] Input (log) initial or simulation parameter
#' values. See Details.
#' @param optimizer_list [list()] Input optimizer and \code{OpenMP} values. See
#' Details.
#' @param nlminb_control [list()] See [mmmControl()]
#'
#' @details
#' Describe \code{template_2d}...
#' Describe \code{result_units}...
#' Describe \code{time_pattern}...
#' Describe \code{capture_map_2d}...
#' Describe \code{structure_list}...
#' Describe \code{parameter_list}...
#' Describe \code{optimizer_list}...
#' Describe \code{simulate_value}...
#'
#' @return An [list()] of class \code{mmmTMB} with elements:
#' \itemize{
#'   \item{\code{adfun}: The TMB AD model object}
#'   \item{\code{model}: Output from [nlminb()]}
#'   \item{\code{structure}: A [list()] of data and constants}
#'   \item{\code{parameter}: A [list()] of initial (log) parameter values}
#'   \item{\code{optimizer}: A [list()] of optimizer constants}
#'   \item{\code{results}: A [list()] of results data frames}
#'   \item{\code{map}: A [list()] of parameter map [factor()]s}
#'   \item{\code{map}: A [list()] of random effects}
#'   \item{\code{sd_report}: A [list()] generated by [sdreport()]}
#'   \item{\code{conv}: A [list()] of convergence diagnostics}
#'   \item{\code{mgc}: The maximum absolute gradient component}
#' }
#'
#' @export
#'
#' @examples
#' # Fit
#' fit_list <- mmmTMB(
#'   released_3d = sim_released_3d,
#'   recovered_5d = sim_recovered_5d,
#'   capture_rate_2d = sim_capture_rate_2d,
#'   report_ratio_2d = sim_report_ratio_2d,
#'   tag_loss_rate = 0.02,
#'   imm_loss_ratio = 0.1,
#'   template_2d = sim_template_2d,
#'   openmp_cores = floor(parallel::detectCores() / 2))
#'
#'
mmmTMB <- function (released_3d, # Data
                    recovered_5d,
                    capture_rate_2d,
                    report_ratio_2d,
                    tag_loss_rate,
                    imm_loss_ratio,
                    template_2d, # Structure
                    recapture_delay = 1,
                    error_family = 1,
                    result_units = 1,
                    time_process = 0,
                    time_pattern = 0,
                    pattern_size = 0,
                    newton_steps = 0, # Optimizer
                    nlminb_loops = 0,
                    openmp_cores = NULL,
                    simulate_value = NULL,
                    capture_map_2d = NULL,
                    structure_list = NULL, # Optional lists
                    parameter_list = NULL,
                    optimizer_list = NULL,
                    nlminb_control = mmmControl()) {


  #---------------- Start the clock -------------------------------------------#

  tictoc::tic("mmmTMB")

  #---------------- Unpack arguments ------------------------------------------#

  if (!is.null(structure_list)) {
    released_3d <- structure_list$released_3d
    template_2d <- structure_list$template_2d
    error_family <- structure_list$error_family
  }
  if (!is.null(optimizer_list)) {
    newton_steps <- optimizer_list$newton_steps
    nlminb_loops <- optimizer_list$nlminb_loops
    openmp_cores <- optimizer_list$openmp_cores
  }

  #---------------- Check arguments: Part 1 of 2 ------------------------------#

  # TODO: check that structure_list has movement_index


  #---------------- Assign dimensions -----------------------------------------#

  np <- sum(template_2d)
  nt <- dim(released_3d)[1]
  na <- dim(released_3d)[2]
  ng <- dim(released_3d)[3]

  #---------------- Check arguments: Part 2 of 2 ------------------------------#



  #---------------- Define the number of cores --------------------------------#

  if (!is.null(openmp_cores)) {TMB::openmp(n = openmp_cores)}

  #---------------- Define the movement index ---------------------------------#

  if (time_process == 0) {
    nv <- 1
    movement_index <- rep(0, nt)
  } else {
    if (time_pattern == 0) {
      nv <- 1
      movement_index <- rep(0, nt)
    } else if (time_pattern == 1) {
      nv <- pattern_size
      movement_index <- rep(seq_len(nv) - 1, each = ceiling(nt/nv))[seq_len(nt)]
    } else if (time_pattern == 2) {
      nv <- pattern_size
      movement_index <- rep((seq_len(nv) - 1), ceiling(nt/nv))[seq_len(nt)]
    } else {
      warning("time_pattern not implemented")
    }
  }

  #---------------- Create the data list --------------------------------------#

  cat("creating tmb_data \n")
  if (!is.null(structure_list)) {
    tmb_data <- structure_list
  } else {
    tmb_data <- list(released_3d = released_3d,
                     recovered_5d = recovered_5d,
                     capture_rate_2d = capture_rate_2d,
                     report_ratio_2d = report_ratio_2d,
                     tag_loss_rate = tag_loss_rate,
                     imm_loss_ratio = imm_loss_ratio,
                     template_2d = template_2d,
                     recapture_delay = recapture_delay,
                     error_family = error_family,
                     result_units = result_units,
                     time_process = time_process,
                     movement_index = movement_index)
  }

  #---------------- Create the parameter list ---------------------------------#

  cat("creating tmb_parameters \n")
  if (!is.null(parameter_list)) {
    tmb_parameters <- parameter_list
  } else {
    tmb_parameters <- list(
      movement_parameters_3d = array(0, dim = c(nv, np, ng)),
      log_capture_bias_2d = array(0, dim = c(na, ng)),
      log_natural_mortality = 0,
      log_dispersion = 0
      # TODO: Add parameters for time-varying
    )
  }

  #---------------- Define map list -------------------------------------------#

  cat("creating tmb_map \n")
  tmb_map <- list()

  # Augment by capture bias map
  if (!is.null(capture_map_2d)) {
    tmb_map <- c(tmb_map, list(log_capture_bias_2d = as.factor(capture_map_2d)))
  } else {
    tmb_map <- c(tmb_map, list(log_capture_bias_2d = factor(rep(1, na * ng))))
  }

  # Map off log_dispersion?
  if (error_family == 0) {
    tmb_map <- c(tmb_map, list(log_dispersion = as.factor(NA)))
  }

  #---------------- Define random effects -------------------------------------#

  # Initialize tmb_random
  tmb_random <- character(0)

  #---------------- Create the simulation object ------------------------------#

  # Use isTRUE(simulate_value) to accommodate NULL

  #---------------- Create the model object -----------------------------------#
  # TODO: Add previous fit option

  cat("creating tmb_obj \n")
  tictoc::tic("creating tmb_obj")
  tmb_obj <- TMB::MakeADFun(data = tmb_data,
                            parameters = tmb_parameters,
                            map = tmb_map,
                            random = tmb_random,
                            DLL = "mmmTMB")
  tictoc::toc()

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
      cat("\nrunning extra nlminb loops \n")
      for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
        cat(paste0("running nlminb loop #", i, "\n"))
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

    # TODO: Troubleshoot Newton steps: output fed into next input correctly?
    # Perform additional Newton steps
    if (newton_steps > 0) {
      cat("\nrunning newton steps \n")
      for (i in seq_len(newton_steps)) {
        cat(paste0("running newton step #", i, "\n"))
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

  #---------------- Create optimizer list -------------------------------------#

  if (is.null(optimizer_list)) {
    optimizer_list <- list(newton_steps = newton_steps,
                           nlminb_loops = nlminb_loops,
                           openmp_cores = openmp_cores)
  }

  #---------------- Create movement probability results -----------------------#

  tictoc::tic("computing movement estimates and std errs")
  pars <- subset_by_name(tmb_opt$par, "movement_parameters_3d")
  if (time_process == 0) {
    covs <- subset_by_name(sd_report$cov.fixed, "movement_parameters_3d")
  } else {
    # warning("Parameters are a random effect: Figure out where to access covs")
    covs <- subset_by_name(sd_report$cov.fixed, "movement_parameters_3d")
  }
  mpr_df <- create_movement_probability_results(
    pars = pars,
    covs = covs,
    dims = c(nv, np, nt, na, ng),
    tp_2d = template_2d,
    result_units = result_units,
    n_draws = 1000)
  tictoc::toc()

  #---------------- Create natural mortality results --------------------------#

  nmr_mat <- summary(sd_report)["natural_mortality_results", , drop = FALSE]
  rownames(nmr_mat) <- NULL
  colnames(nmr_mat) <- c("Estimate", "SE")
  nmr_df <- as.data.frame(nmr_mat)

  #---------------- Create capture bias results -------------------------------#

  cbr_ind <- which(rownames(summary(sd_report)) == "capture_bias_2d")
  cbr_mat <- matrix(0, nrow = na * ng, ncol = 4)
  cbr_mat[, 1] <- rep(seq_len(ng), each = na)
  cbr_mat[, 2] <- rep(seq_len(na), ng)
  cbr_mat[, 3:4] <- summary(sd_report)[cbr_ind, ]
  colnames(cbr_mat) <- c("Class", "Area", "Estimate", "SE")
  cbr_df <- as.data.frame(cbr_mat)

  #---------------- Create dispersion results ---------------------------------#

  dsp_mat <- summary(sd_report)["dispersion", , drop = FALSE]
  rownames(dsp_mat) <- NULL
  colnames(dsp_mat) <- c("Estimate", "SE")
  dsp_df <- as.data.frame(dsp_mat)

  #---------------- Create results list ---------------------------------------#

  results_list <- list(movement_probability_results = mpr_df,
                       natural_mortality_results = nmr_df,
                       capture_bias = cbr_df,
                       dispersion = dsp_df)

  #---------------- Compute AIC -----------------------------------------------#

  num_doubles_or_na <- sum(unlist(lapply(
    tmb_map,
    function(x) length(x) - length(unique(x)) + length(which(is.na(unique(x))))
  )))
  k <- length(unlist(parameter_list)) - num_doubles_or_na
  nll <- tmb_opt$objective
  aic <- 2 * k + 2 * nll

  #---------------- Stop the clock --------------------------------------------#

  tictoc::toc()

  #---------------- Return an mmmTMB object -----------------------------------#

  structure(list(
    adfun      = tmb_obj, # report(),
    model      = tmb_opt, # convergence, objective,
    structure  = tmb_data, # TODO: optionally place in simulated data
    parameter  = tmb_parameters, # Initial values
    optimizer  = optimizer_list,
    results    = results_list,
    map        = tmb_map,
    random     = tmb_random,
    sd_report  = sd_report, # value, sd, cov, par.fixed, cov.fixed, pdHess, gradient.fixed, env
    conv       = conv, # final_grads, bad_eig
    mgc        = mgc,
    aic        = aic),
    class      = "mmmTMB")
}


#' Settings for \code{mmmFit()}
#'
#' @param error_family Character. One of \code{"nb1"} or \code{"poisson"}
#' @param time_varying Logical. Should movement rates vary through time?
#' @param time_process Character. One of \code{"none"} or \code{"rw"}
#' @param cycle_length Integer. Cycle length or \code{0} for no cycle.
#' @param block_length Integer. Block length for movement rate.
#' @param fish_rate_by Character. Estimate F as a single value (\code{"none"}),
#'   by time steps as blocks (\code{"block"}), by areas (\code{"areas"}), or by
#'   both blocks and areas (\code{"both"}).
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
mmmSet <- function (error_family = c("nb1", "poisson"),
                    time_varying = 0,
                    time_process = c("none", "rw"),
                    cycle_length = 0,
                    block_length = 1,
                    fish_rate_by = c("none", "block", "area", "both"),
                    results_step = 1,
                    nlminb_loops = 5,
                    newton_steps = 5,
                    openmp_cores = NULL) {

  #--------------- Return a list of settings ----------------------------------#

  structure(list(
    error_family = error_family[1],
    time_varying = time_varying,
    time_process = time_process[1],
    cycle_length = cycle_length,
    block_length = block_length,
    fish_rate_by = fish_rate_by[1],
    results_step = results_step,
    nlminb_loops = nlminb_loops,
    newton_steps = newton_steps,
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
