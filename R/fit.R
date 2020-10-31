## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib mmmTMB, .registration = TRUE
## usethis namespace: end
NULL

#' Fit a Markovian Movement Model
#'
#' @description Estimate movement rates
#'
#' @param data A list of named data objects. See Details.
#' @param parameters A list of initial parameter values. See Details.
#' @param random A character vector of parameters to estimate as random effects.
#' @param map A list of \code{map} values to override defaults.
#'   See \code{TMB::MakeADFun()}.
#' @param settings A list of values that define the model.
#'   See \code{mmmSet()}.
#' @param control A list of optimization control options to pass to
#'   \code{stats::nlminb()}. See \code{mmmControl()}.
#'
#' @details The list argument \code{data} must contain
#' \itemize{
#'   \item{\code{tags} a list of class \code{mmmTags} (recommended: see
#'     \code{mmmTags()}) OR both \code{mT} and \code{mR} integer matrices
#'     of tag release and recovery counts, respectively, formatted as
#'     described in \code{mmmTags()}}
#'   \item{\code{mI} a square binary index matrix representing movement
#'     between areas (from rows to columns). Ones represent movement that
#'     is allowed from one time step to the next. Zeros off the diagonal
#'     represent disallowed movement. Numbers on the diagonal are ignored
#'     because self-movement is always allowed.}
#' }
#' The list argument \code{data} may optionally contain
#' \itemize{
#'   \item{\code{mL} a matrix of tag reporting rates (proportions)
#'     in which the time step is given by the row and the area is
#'     given by the column. Defaults to ones (full reporting).}
#'   \item{mW}
#' }
#'
#'
#' @return A list of class \code{mmmFit}.
#' @export
#'
#' @author Luke A. Rogers
#'
#' @examples
#'
mmmFit <- function(data,
                   parameters = NULL,
                   random = NULL,
                   map = NULL,
                   settings = mmmSet(),
                   control = mmmControl()) {

  # TODO: Convert mI to transpose

  #---------------- Start the clock -------------------------------------------#

  tictoc::tic("mmmFit")

  #---------------- Check data argument ---------------------------------------#

  cat("checking data arguments \n")
  checkmate::assert_list(data, null.ok = FALSE)
  checkmate::assert_class(data$tags, "mmmTags", null.ok = TRUE)
  checkmate::assert_matrix(data$mT, mode = "integerish", null.ok = TRUE)
  checkmate::assert_matrix(data$mR, mode = "integerish", null.ok = TRUE)
  checkmate::assert_matrix(data$mI, mode = "integerish", null.ok = TRUE)
  checkmate::assert_matrix(data$mL, mode = "double", null.ok = TRUE)
  checkmate::assert_matrix(data$mW, mode = "double", null.ok = TRUE)
  # Optionally parameters
  checkmate::assert_matrix(data$mF, mode = "double", null.ok = TRUE)
  checkmate::assert_double(data$sM, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(data$sH, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(data$sC, lower = 0, len = 1, null.ok = TRUE)

  #---------------- Check parameters argument ---------------------------------#

  cat("checking parameter arguments \n")
  checkmate::assert_list(parameters, null.ok = TRUE)
  # Optionally data
  checkmate::assert_array(parameters$aP, mode = "double", null.ok = TRUE)
  checkmate::assert_matrix(parameters$mF, mode = "double", null.ok = TRUE)
  checkmate::assert_double(parameters$sM, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(parameters$sH, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(parameters$sC, lower = 0, len = 1, null.ok = TRUE)
  # Optionally specified initial values
  checkmate::assert_numeric(parameters$vB, lower = 0, null.ok = TRUE)

  #---------------- Check settings argument -----------------------------------#

  cat("checking settings arguments \n")
  checkmate::assert_list(settings, null.ok = FALSE)
  # Optionally specified settings
  checkmate::assert_string(settings$error_family, null.ok = TRUE)
  checkmate::assert_integerish(settings$time_varying, null.ok = FALSE)
  checkmate::assert_string(settings$time_process, null.ok = FALSE)
  checkmate::assert_integerish(settings$cycle_length, null.ok = FALSE)
  checkmate::assert_integerish(settings$block_length, null.ok = FALSE)
  checkmate::assert_string(settings$fish_rate_by, null.ok = FALSE)
  checkmate::assert_integerish(settings$results_step, null.ok = FALSE)
  checkmate::assert_integerish(settings$nlminb_loops, null.ok = FALSE)
  checkmate::assert_integerish(settings$newton_steps, null.ok = FALSE)
  checkmate::assert_integerish(settings$openmp_cores, null.ok = TRUE)

  #---------------- Set default settings --------------------------------------#

  # Error family
  if (is.null(settings$error_family)) { error_family <- 1L }
  else if (is.na(settings$error_family)) { error_family <- 1L }
  else if (settings$error_family == "poisson") { error_family <- 0L }
  else { error_family <- 1L }
  # Span liberty
  if (is.null(settings$span_liberty)) { span_liberty <- c(1, Inf) }
  else if (is.na(settings$span_liberty)) { span_liberty <- c(1, Inf) }
  else { span_liberty <- settings$span_liberty }
  # Time varying
  if (is.null(settings$time_varying)) { time_varying <- 0L }
  else if (is.na(settings$time_varying)) { time_varying <- 0L }
  else if (settings$time_varying == 1L) { time_varying <- 1L }
  else { time_varying <- 0L }
  # Time process
  if (is.null(settings$time_process)) { time_process <- 0L }
  else if (is.na(settings$time_process)) { time_process <- 0L }
  else if (settings$time_process == "rw") { time_process <- 1L }
  else { time_process <- 0L }
  # Cycle length
  if (is.null(settings$cycle_length)) { cycle_length <- 0L }
  else if (is.na(settings$cycle_length)) { cycle_length <- 0L }
  else if (settings$cycle_length > 0L) { cycle_length <- settings$cycle_length }
  else { cycle_length <- 0L }
  # Block length
  if (is.null(settings$block_length)) { block_length <- 0L }
  else if (is.na(settings$block_length)) { block_length <- 0L }
  else if (settings$block_length > 0L) { block_length <- settings$block_length }
  else { block_length <- 0L }
  # Fish rate by
  if (is.null(settings$fish_rate_by)) { fish_rate_by <- "none" }
  else if (is.na(settings$fish_rate_by)) { fish_rate_by <- "none" }
  else if (settings$fish_rate_by == "none") { fish_rate_by <- "none" }
  else if (settings$fish_rate_by == "block") { fish_rate_by <- "block" }
  else if (settings$fish_rate_by == "area") { fish_rate_by <- "area" }
  else if (settings$fish_rate_by == "both") { fish_rate_by <- "both" }
  else { fish_rate_by <- "none" }
  # Results step
  if (is.null(settings$results_step)) { results_step <- 1L }
  else if (is.na(settings$results_step)) { results_step <- 1L }
  else if (settings$results_step > 1L) { results_step <- settings$results_step }
  else { results_step <- 1L }
  # Nlminb loops
  if (is.null(settings$nlminb_loops)) { nlminb_loops <- 5L }
  else if (is.na(settings$nlminb_loops)) { nlminb_loops <- 5L }
  else if (settings$nlminb_loops > 0L) { nlminb_loops <- settings$nlminb_loops }
  else { nlminb_loops <- 5L }
  # Newton steps
  if (is.null(settings$newton_steps)) { newton_steps <- 5L }
  else if (is.na(settings$newton_steps)) { newton_steps <- 5L }
  else if (settings$newton_steps > 0L) { newton_steps <- settings$newton_steps }
  else { newton_steps <- 5L }
  # OpenMP cores
  if (is.null(settings$openmp_cores)) { openmp_cores <- 1L }
  else if (is.na(settings$openmp_cores)) { openmp_cores <- 1L }
  else if (settings$openmp_cores > 1L) { openmp_cores <- settings$openmp_cores }
  else { openmp_cores <- 1L }

  #---------------- Unpack required data --------------------------------------#

  # Tag matrices
  if (!is.null(data$tags)) {
    mT <- data$tags$mT
    mR <- data$tags$mR
  } else {
    mT <- data$mT
    mR <- data$mR
  }

  #---------------- Check required data ---------------------------------------#

  # Tag releases
  checkmate::assert_matrix(mT, mode = "integerish", any.missing = FALSE)
  checkmate::assert_matrix(mT, ncols = 4, null.ok = FALSE)
  checkmate::assert_true(all(mT >= 0))
  checkmate::assert_true(colnames(mT)[1] == "release_step")
  checkmate::assert_true(colnames(mT)[2] == "release_area")
  checkmate::assert_true(colnames(mT)[3] == "group")
  checkmate::assert_true(colnames(mT)[4] == "count")
  # Tag recoveries
  checkmate::assert_matrix(mR, mode = "integerish", any.missing = FALSE)
  checkmate::assert_matrix(mR, ncols = 6, null.ok = FALSE)
  checkmate::assert_true(all(mR >= 0))
  checkmate::assert_true(colnames(mR)[1] == "release_step")
  checkmate::assert_true(colnames(mR)[2] == "release_area")
  checkmate::assert_true(colnames(mR)[3] == "recover_step")
  checkmate::assert_true(colnames(mR)[4] == "recover_area")
  checkmate::assert_true(colnames(mR)[5] == "group")
  checkmate::assert_true(colnames(mR)[6] == "count")
  # Check values
  if (min(mT$release_step) > 0) cat("caution: time step not indexed from zero")

  #---------------- Create index limits ---------------------------------------#

  cat("creating index limits \n")
  # Set index limits
  nt <- max(c(mT$release_step, mR$recover_step)) + 1L # Index from zero for C++
  na <- max(c(mT$release_area, mR$recover_area)) + 1L # Index from zero for C++
  ng <- max(c(mT$group, mR$group)) + 1L # Index from zero for C++

  #---------------- Unpack index matrix ---------------------------------------#

  if (!is.null(data$mI)) {
    mI <- data$mI
    diag(mI) <- 0
  } else {
    mI <- matrix(1L, nrow = na, ncol = na)
    diag(mI) <- 0L
  }

  #---------------- Check index matrix ----------------------------------------#

  checkmate::assert_matrix(mI, mode = "integerish", null.ok = FALSE)
  checkmate::assert_matrix(mI, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_true(nrow(mI) == ncol(mI))
  checkmate::assert_true(all(is.element(mI, 0L:1L)))

  #---------------- Create parameter index limit ------------------------------#

  # Set parameter index limit
  np <- sum(mI)

  #---------------- Create parameter indexes ----------------------------------#

  # Set for movement parameters
  if (time_varying) {
    if (!block_length && !cycle_length) {
      npt <- nt
      vpt <- c(seq_len(nt) - 1L) # Index from zero for C++
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
    vpt <- rep(0L, nt) # Index from zero for C++
  }

  #---------------- Create fishing rate indexes -------------------------------#

  # Set for fishing rate time step
  if (fish_rate_by == "none" || fish_rate_by == "area") {
    nft <- 1L
    vft <- rep(0L, nt) # Index from zero for C++
  } else {
    if (block_length) {
      nft <- ceiling(nt / block_length)
      vft <- rep(seq_len(nft) - 1L, each = block_length)[seq_len(nt)]
    } else {
      nft <- nt
      vft <- c(seq_len(nt) - 1L) # Index from zero for C++
    }
  }
  # Set for fishing rate areas
  if (fish_rate_by == "none" || fish_rate_by == "block") {
    nfa <- 1L
    vfa <- rep(0L, na) # Index from zero for C++
  } else {
    nfa <- na
    vfa <- c(seq_len(na) - 1L) # Index from zero for C++
  }

  #---------------- Unpack arguments for tmb_data -----------------------------#

  # Tag reporting rate
  if (!is.null(data$mL)) { mL <- data$mL }
  else { mL <- matrix(1, nrow = 1L, ncol = 1L) }
  # Fishing rate weighting
  if (!is.null(data$mW)) { mW <- data$mW }
  else { mW <- matrix(1, nrow = 1L, ncol = 1L) }

  #---------------- Check arguments for tmb_data ------------------------------#

  # Tag reporting rate
  checkmate::assert_matrix(mL, mode = "double", null.ok = FALSE)
  checkmate::assert_matrix(mL, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_true(ncol(mL) == 1 || ncol(mL) == na)
  checkmate::assert_true(all(mL >= 0))
  # Fishing rate weighting
  checkmate::assert_matrix(mW, mode = "double", null.ok = FALSE)
  checkmate::assert_matrix(mW, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_true(ncol(mW) == 1 || ncol(mW) == na)
  checkmate::assert_true(all(mW >= 0))
  checkmate::assert_true(all(colSums(mW) == 1))

  #---------------- Create tag reporting rate indexes -------------------------#

  # Set for time step
  if (nrow(mL) == 1L) {
    nlt <- 1L
    vlt <- rep(0L, nt)
  } else {
    nlt <- nt
    vlt <- c(seq_len(nt) - 1L)
  }
  # Set for areas
  if (ncol(mL) == 1L) {
    nla <- 1L
    vla <- rep(0L, na)
  } else {
    nla <- na
    vla <- c(seq_len(na) - 1L)
  }

  #---------------- Create fishing rate weighting indexes ---------------------#

  # Set for time step
  if (nrow(mW) == 1L) {
    nwt <- 1L
    vwt <- rep(0L, nt)
  } else {
    nwt <- nrow(mW)
    vwt <- rep(c(seq_len(nwt) - 1L), ceiling(nt / nwt))[seq_len(nt)]
  }
  # Set for areas
  if (ncol(mW) == 1L) {
    nwa <- 1L
    vwa <- rep(0L, na)
  } else {
    nwa <- na
    vwa <- c(seq_len(na) - 1L)
  }

  #---------------- Update time span at liberty -------------------------------#

  if (!is.null(data$tags)) { init_liberty <- data$tags$init_liberty }
  else { init_liberty <- 1L }
  if (span_liberty[1] < init_liberty) {
    span_liberty[1] <- init_liberty
    cat(paste0("using initial time at liberty from mmmTags(): ", init_liberty))
  }

  #---------------- Unpack arguments for tmb_parameters -----------------------#

  # Array movement parameters
  if (!is.null(parameters$aP)) { aP <- parameters$aP }
  else { aP <- array(0, dim = c(npt, np, ng)) }
  # Matrix log fishing mortality rate
  if (!is.null(data$mF)) { lmF <- log(data$mF) }
  else if (!is.null(parameters$mF)) { lmF <- log(parameters$mF) }
  else { lmF <- array(-3, dim = c(nft, nfa)) }
  # Vector log fishing selectivity and tag reporting bias
  if (!is.null(parameters$vB)) { lvB <- log(parameters$vB) }
  else { lvB <- numeric(length = ng) }
  # Scalar log tag loss rate
  if (!is.null(data$sH)) { lsH <- log(data$sH) }
  if (!is.null(parameters$sH)) { lsH <- log(parameters$sH) }
  else { lsH <- -3L }
  # Scalar log initial tag loss rate
  if (!is.null(data$sC)) { lsC <- log(data$sC) }
  if (!is.null(parameters$sC)) { lsC <- log(parameters$sC) }
  else { lsC <- -2L }
  # Scalar log natural mortality
  if (!is.null(data$sM)) { lsM <- log(data$sM) }
  if (!is.null(parameters$sM)) { lsM <- log(parameters$sM) }
  else { lsM <- -2L }
  # Scalar log negative binomial dispersal
  lsD <- 0L

  #---------------- Check arguments for tmb_parameters ------------------------#

  # Movement parameter array
  checkmate::assert_array(aP, mode = "double", any.missing = FALSE, d = 3)
  checkmate::assert_true(dim(aP)[1] == npt, na.ok = FALSE)
  checkmate::assert_true(dim(aP)[2] == np, na.ok = FALSE)
  checkmate::assert_true(dim(aP)[3] == ng, na.ok = FALSE)
  # Fishing rate matrix
  checkmate::assert_matrix(mF, mode = "double", any.missing = FALSE)
  checkmate::assert_true(dim(mF)[1] == nft, na.ok = FALSE)
  checkmate::assert_true(dim(mF)[2] == nfa, na.ok = FALSE)

  #---------------- Create tmb_data -------------------------------------------#

  cat("creating tmb_data \n")
  tmb_data <- list(
    mT = mT,
    mR = mR,
    mI = t(mI),
    mL = mL,
    mW = mW,
    error_family = error_family,
    span_liberty = span_liberty,
    time_varying = time_varying,
    time_process = time_process,
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
    vwa = vwa # Index vector: area for fishing rate weighting
  )

  #---------------- Create tmb_parameters -------------------------------------#

  cat("creating tmb_parameters \n")
  tmb_parameters <- list(
    aP = aP, # Array: movement parameters
    lmF = lmF, # Matrix: log fishing mortality rates
    lvB = lvB, # Vector: log fishery selectivity and tag reporting bias
    lsM = lsM, # Scalar: log natural mortality
    lsH = lsH, # Scalar: log tag loss rate
    lsC = lsC, # Scalar: log initial tag loss rate
    lsD = lsD # Scalar: negative binomial dispersion
  )

  #---------------- Create tmb_random -----------------------------------------#

  cat("creating tmb_random \n")
  if (is.null(random)) { tmb_random <- character() }
  else { tmb_random <- random }

  #---------------- Create tmb_map --------------------------------------------#

  cat("creating tmb_map \n")
  tmb_map <- list()
  # Default
  if (!is.null(data$mF)) { tmb_map <- c(tmb_map, lmF = NA) }
  if (!is.null(data$sH)) { tmb_map <- c(tmb_map, lsH = NA) }
  if (!is.null(data$sC)) { tmb_map <- c(tmb_map, lsC = NA) }
  if (!is.null(data$sM)) { tmb_map <- c(tmb_map, lsM = NA) }
  if (error_family == 0) { tmb_map <- c(tmb_map, lsD = NA)}
  # User defined
  if (!is.null(map$mF)) { tmb_map$lmF <- map$mF }
  if (!is.null(map$vB)) { tmb_map$lvB <- map$vB }
  if (!is.null(map$sH)) { tmb_map$lsH <- map$sH }
  if (!is.null(map$sC)) { tmb_map$lsC <- map$sC }
  if (!is.null(map$sM)) { tmb_map$lsM <- map$sM }

  #---------------- Set the number of OpenMP cores ----------------------------#

  cat(paste0("using ", openmp_cores, " openmp cores"))
  TMB::openmp(n = openmp_cores)

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
  # TODO: Update results
  # Use result_steps here

  #---------------- Create movement probability results -----------------------#

  tictoc::tic("computing movement estimates and std errs")
  pars <- subset_by_name(model$par, "movement_parameters_3d")
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

  tictoc::toc()

  #----------------------------------------------------------------------------#

  tictoc::toc()

  #---------------- Stop the clock --------------------------------------------#

  tictoc::toc()

  #---------------- Return an mmmFit object -----------------------------------#

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
    class       = c("mmmFit"))
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
#' @param span_liberty Integer. Min and max time steps at liberty before
#'   recapture. The initial time at liberty \code{span_liberty[1]} should
#'   agree in duration with \code{days_liberty} from \code{mmmTags()}.
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
                    span_liberty = c(1, Inf),
                    time_varying = 0,
                    time_process = c("none", "rw"),
                    cycle_length = 0,
                    block_length = 0,
                    fish_rate_by = c("none", "block", "area", "both"),
                    results_step = 1,
                    nlminb_loops = 5,
                    newton_steps = 5,
                    openmp_cores = NULL) {

  #--------------- Return a list of settings ----------------------------------#

  structure(list(
    error_family = error_family[1],
    span_liberty = span_liberty[1:2],
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
