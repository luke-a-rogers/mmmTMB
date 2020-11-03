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
#'   \item{mF}
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
  checkmate::assert_numeric(data$mT, lower = 0, null.ok = TRUE)
  checkmate::assert_numeric(data$mR, lower = 0, null.ok = TRUE)
  checkmate::assert_numeric(data$mI, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_numeric(data$mL, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_numeric(data$mW, lower = 0, upper = 0, null.ok = TRUE)
  # Optionally parameters
  checkmate::assert_matrix(data$mF, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(data$mF, lower = 0, null.ok = TRUE)
  checkmate::assert_double(data$sM, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(data$sH, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(data$sC, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(data$sC, upper = 1, null.ok = TRUE)

  #---------------- Check parameters argument ---------------------------------#

  cat("checking parameter arguments \n")
  checkmate::assert_list(parameters, null.ok = TRUE)
  # Optionally data
  checkmate::assert_array(parameters$aP, mode = "double", null.ok = TRUE)
  checkmate::assert_matrix(parameters$mF, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(parameters$mF, lower = 0, null.ok = TRUE)
  checkmate::assert_double(parameters$sM, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(parameters$sH, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(parameters$sC, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(parameters$sC, upper = 1, null.ok = TRUE)
  # Optionally specified initial values
  checkmate::assert_numeric(parameters$vB, lower = 0, null.ok = TRUE)

  #---------------- Check settings argument -----------------------------------#

  cat("checking settings arguments \n")
  # Check settings
  checkmate::assert_list(settings, null.ok = FALSE)
  # Check error family setting
  checkmate::assert_integerish(settings$error_family, len = 1)
  checkmate::assert_integerish(settings$error_family, lower = 0, upper = 1)
  checkmate::assert_integerish(settings$error_family, any.missing = FALSE)
  # Check span liberty setting
  checkmate::assert_integerish(settings$span_liberty, len = 2)
  checkmate::assert_integerish(settings$span_liberty, any.missing = FALSE)
  checkmate::assert_integerish(settings$span_liberty, lower = 0)
  checkmate::assert_true(settings$span_liberty[1] < settings$time_varying[2])
  # Check time varying setting
  checkmate::assert_integerish(settings$time_varying, len = 1)
  checkmate::assert_integerish(settings$time_varying, lower = 0, upper = 1)
  checkmate::assert_integerish(settings$time_varying, any.missing = FALSE)
  # Check time process setting
  checkmate::assert_integerish(settings$time_process, len = 1)
  checkmate::assert_integerish(settings$time_process, lower = 0, upper = 1)
  checkmate::assert_integerish(settings$time_process, any.missing = FALSE)
  # Check cycle length setting
  checkmate::assert_integerish(settings$cycle_length, len = 1)
  checkmate::assert_integerish(settings$cycle_length, lower = 0)
  checkmate::assert_integerish(settings$cycle_length, any.missing = FALSE)
  # Check block length setting
  checkmate::assert_integerish(settings$block_length, len = 1)
  checkmate::assert_integerish(settings$block_length, lower = 0)
  checkmate::assert_integerish(settings$block_length, any.missing = FALSE)
  # Check fish rate by setting
  checkmate::assert_integerish(settings$fish_rate_by, len = 1)
  checkmate::assert_integerish(settings$fish_rate_by, lower = 0, upper = 3)
  checkmate::assert_integerish(settings$fish_rate_by, any.missing = FALSE)
  # Check results step setting
  checkmate::assert_integerish(settings$results_step, len = 1)
  checkmate::assert_integerish(settings$results_step, lower = 1)
  checkmate::assert_integerish(settings$results_step, any.missing = FALSE)
  # Check nlminb loops setting
  checkmate::assert_integerish(settings$nlminb_loops, len = 1)
  checkmate::assert_integerish(settings$nlminb_loops, lower = 0)
  checkmate::assert_integerish(settings$nlminb_loops, any.missing = FALSE)
  # Check newton steps setting
  checkmate::assert_integerish(settings$newton_steps, len = 1)
  checkmate::assert_integerish(settings$newton_steps, lower = 0)
  checkmate::assert_integerish(settings$newton_steps, any.missing = FALSE)
  # Check OpenMP cores setting
  checkmate::assert_integerish(settings$openmp_cores, len = 1)
  checkmate::assert_integerish(settings$openmp_cores, lower = 1)
  checkmate::assert_integerish(settings$openmp_cores, any.missing = FALSE)

  #---------------- Set default settings --------------------------------------#

  # Error family
  if (is.null(settings$error_family)) {
    error_family <- 1L
  } else {
    error_family <- settings$error_family
  }
  # Span liberty
  if (is.null(settings$span_liberty)) {
    span_liberty <- c(1, Inf)
  } else {
    span_liberty <- settings$span_liberty
  }
  # Time varying
  if (is.null(settings$time_varying)) {
    time_varying <- 0L
  } else {
    time_varying <- settings$time_varying
  }
  # Time process
  if (is.null(settings$time_process)) {
    time_process <- 0L
  } else {
    time_process <- settings$time_process
  }
  # Cycle length
  if (is.null(settings$cycle_length)) {
    cycle_length <- 0L
  } else {
    cycle_length <- settings$cycle_length
  }
  # Block length
  if (is.null(settings$block_length)) {
    block_length <- 0L
  } else {
    block_length <- settings$block_length
  }
  # Fish rate by
  if (is.null(settings$fish_rate_by)) {
    fish_rate_by <- 0L
  } else {
    fish_rate_by <- settings$fish_rate_by
  }
  # Results step
  if (is.null(settings$results_step)) {
    results_step <- 1L
  } else {
    results_step <- settings$results_step
  }
  # Number of nlminb loops
  if (is.null(settings$nlminb_loops)) {
    nlminb_loops <- 5L
  } else {
    nlminb_loops <- settings$nlminb_loops
  }
  # Newton steps
  if (is.null(settings$newton_steps)) {
    newton_steps <- 5L
  } else {
    newton_steps <- settings$newton_steps
  }
  # OpenMP cores
  if (is.null(settings$openmp_cores)) {
    openmp_cores <- 1L
  } else {
    openmp_cores <- settings$openmp_cores
  }

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
  checkmate::assert_numeric(mT, lower = 0, null.ok = FALSE)
  checkmate::assert_true(colnames(mT)[1] == "release_step")
  checkmate::assert_true(colnames(mT)[2] == "release_area")
  checkmate::assert_true(colnames(mT)[3] == "group")
  checkmate::assert_true(colnames(mT)[4] == "count")
  # Tag recoveries
  checkmate::assert_matrix(mR, mode = "integerish", any.missing = FALSE)
  checkmate::assert_matrix(mR, ncols = 6, null.ok = FALSE)
  checkmate::assert_numeric(mR, lower = 0, null.ok = FALSE)
  checkmate::assert_true(colnames(mR)[1] == "release_step")
  checkmate::assert_true(colnames(mR)[2] == "release_area")
  checkmate::assert_true(colnames(mR)[3] == "group")
  checkmate::assert_true(colnames(mR)[4] == "recover_step")
  checkmate::assert_true(colnames(mR)[5] == "recover_area")
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
  checkmate::assert_numeric(mI, lower = 0, upper = 1, null.ok = FALSE)
  checkmate::assert_true(nrow(mI) == ncol(mI))

  #---------------- Create parameter index limit ------------------------------#

  # Set parameter index limit
  np <- as.integer(sum(mI))

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

  if (!is.null(data$mF)) {
    nft <- nrow(data$mF)
    nfa <- ncol(data$mF)
    vft <- rep(seq_len(nft) - 1L, each = ceiling(nt / nft))[seq_len(nt)]
    vfa <- rep(seq_len(nfa) - 1L, each = ceiling(na / nfa))[seq_len(na)]
    if (nt %% nft) {
      cat("warning: nft does not divide nt evenly \n")
    }
    if (!is.element(nfa, c(1, na))) {
      cat("warning: nfa must equal 1 or na \n")
    }
    cat("using mF from data and ignoring fish_rate_by in settings \n")
  } else if (!is.null(parameters$mF)) {
    nft <- nrow(parameters$mF)
    nfa <- ncol(parameters$mF)
    vft <- rep(seq_len(nft) - 1L, each = ceiling(nt / nft))[seq_len(nt)]
    vfa <- rep(seq_len(nfa) - 1L, each = ceiling(na / nfa))[seq_len(na)]
    if (nt %% nft) {
      cat("warning: nft does not divide nt evenly \n")
    }
    if (!is.element(nfa, c(1, na))) {
      cat("warning: nfa must equal 1 or na \n")
    }
    cat("initial mF from parameters; ignoring fish_rate_by in settings \n")
  } else {
    # Set for fishing rate time step
    if (is.element(fish_rate_by, c(0L, 1L))) { # None or area
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
    if (is.element(fish_rate_by, c(0L, 2L))) { # None or block
      nfa <- 1L
      vfa <- rep(0L, na) # Index from zero for C++
    } else {
      nfa <- na
      vfa <- c(seq_len(na) - 1L) # Index from zero for C++
    }
    cat("estimating mF; using fish_rate_by from settings \n")
  }

  #---------------- Unpack arguments for tmb_data -----------------------------#

  # Tag reporting rate
  if (!is.null(data$mL)) {
    mL <- data$mL
  } else {
    mL <- matrix(1, nrow = 1L, ncol = 1L)
  }
  # Fishing rate weighting
  if (!is.null(data$mW)) {
    mW <- data$mW
  }
  else {
    mW <- matrix(1, nrow = 1L, ncol = 1L)
  }

  #---------------- Check arguments for tmb_data ------------------------------#

  # Tag reporting rate
  checkmate::assert_matrix(mL, mode = "double", null.ok = FALSE)
  checkmate::assert_matrix(mL, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(mL, lower = 0, null.ok = FALSE)
  checkmate::assert_true(ncol(mL) == 1 || ncol(mL) == na)
  # Fishing rate weighting
  checkmate::assert_matrix(mW, mode = "double", null.ok = FALSE)
  checkmate::assert_matrix(mW, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(mW, lower = 0, null.ok = FALSE)
  checkmate::assert_true(ncol(mW) == 1 || ncol(mW) == na)
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

  # Check data$tags for minimum steps at liberty
  if (!is.null(data$tags)) {
    init_liberty <- data$tags$init_liberty
  } else {
    init_liberty <- 1L
  }
  # Reconcile minimum steps at liberty
  if (span_liberty[1] < init_liberty) {
    span_liberty[1] <- init_liberty
    cat(paste0("using initial time at liberty from mmmTags(): ", init_liberty))
  }

  #---------------- Unpack arguments for tmb_parameters -----------------------#

  # Array movement parameters
  if (!is.null(parameters$aP)) {
    aP <- parameters$aP
  } else {
    aP <- array(0, dim = c(npt, np, ng))
  }
  # Matrix fishing mortality rate
  if (!is.null(data$mF)) {
    mF <- data$mF
  } else if (!is.null(parameters$mF)) {
    mF <- parameters$mF
  } else {
    mF <- array(3L, dim = c(nft, nfa))
  }
  # Scalar natural mortality
  if (!is.null(data$sM)) {
    sM <- data$sM
  } else if (!is.null(parameters$sM)) {
    sM <- parameters$sM
  } else {
    sM <- 3L
  }
  # Scalar tag loss rate
  if (!is.null(data$sH)) {
    sH <- data$sH
  } else if (!is.null(parameters$sH)) {
    sH <- parameters$sH
  } else {
    sH <- 5
  }
  # Scalar initial tag loss rate (proportion)
  if (!is.null(data$sC)) {
    sC <- data$sC
  } else if (!is.null(parameters$sC)) {
    sC <- parameters$sC
  } else {
    sC <- 0
  }
  # Vector fishing bias
  if (!is.null(parameters$vB)) {
    vB <- parameters$vB
  } else {
    vB <- rep(1, ng)
  }
  # Scalar negative binomial dispersal
  if (!is.null(parameters$sD)) {
    sD <- parameters$sD
  } else {
    sD <- 1L
  }
  #---------------- Check arguments for tmb_parameters ------------------------#

  # Movement parameter array
  checkmate::assert_array(aP, mode = "double", any.missing = FALSE, d = 3)
  checkmate::assert_true(dim(aP)[1] == npt, na.ok = FALSE)
  checkmate::assert_true(dim(aP)[2] == np, na.ok = FALSE)
  checkmate::assert_true(dim(aP)[3] == ng, na.ok = FALSE)
  # Fishing rate matrix
  checkmate::assert_matrix(mF, mode = "double", any.missing = FALSE)
  checkmate::assert_numeric(mF, lower = 0, null.ok = TRUE)
  checkmate::assert_true(dim(mF)[1] == nft, na.ok = FALSE)
  checkmate::assert_true(dim(mF)[2] == nfa, na.ok = FALSE)

  #---------------- Create tmb_data -------------------------------------------#

  cat("creating tmb_data \n")
  tmb_data <- list(
    tmT = t(mT),
    tmR = t(mR),
    tmI = t(mI),
    tmL = t(mL),
    tmW = t(mW),
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
    taP = aperm(aP, c(2,1,3)), # Movement parameters
    logit_exp_neg_tmF = logit(exp(-(t(mF)))), # Fishing mortality
    logit_exp_neg_sM = logit(exp(-(sM))), # Natural mortality
    logit_exp_neg_sH = logit(exp(-(sH))), # Tag loss rate
    logit_sA = logit(1 - sC), # Initial tag loss rate (proportion)
    log_vB = log(vB), # Fishing rate bias
    log_sD = log(sD) # Negative binomial dispersion
  )

  #---------------- Create tmb_random -----------------------------------------#

  cat("creating tmb_random \n")
  if (is.null(random)) {
    tmb_random <- character()
  } else {
    tmb_random <- random
  }

  #---------------- Create tmb_map --------------------------------------------#

  cat("creating tmb_map \n")
  tmb_map <- list()
  # Default
  if (!is.null(data$mF)) { tmb_map <- c(tmb_map, logit_exp_neg_tmF = NA) }
  if (!is.null(data$sM)) { tmb_map <- c(tmb_map, logit_exp_neg_sM = NA) }
  if (!is.null(data$sH)) { tmb_map <- c(tmb_map, logit_exp_neg_sH = NA) }
  if (!is.null(data$sC)) { tmb_map <- c(tmb_map, logit_sA = NA) }
  if (!is.null(data$vB)) { tmb_map <- c(tmb_map, log_vB = NA) }
  if (error_family == 0) { tmb_map <- c(tmb_map, log_sD = NA)}
  # User defined
  if (!is.null(map$mF)) { tmb_map$logit_exp_neg_tmF <- t(map$mF) }
  if (!is.null(map$sM)) { tmb_map$logit_exp_neg_sM <- t(map$sM) }
  if (!is.null(map$sH)) { tmb_map$logit_exp_neg_sH <- t(map$sH) }
  if (!is.null(map$sC)) { tmb_map$logit_sA <- map$sC }
  if (!is.null(map$vB)) { tmb_map$log_vB <- map$vB }

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
  # TODO: Continue from here

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
  # TODO: Fix before here

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

#' Settings for \code{mmmFit()}
#'
#' @param error_family Integer. \code{0} = Poisson; \code{1} = NB1.
#' @param span_liberty Integer. Min and max time steps at liberty before
#'   recapture. The initial time at liberty \code{span_liberty[1]} should
#'   agree in duration with \code{days_liberty} from \code{mmmTags()}.
#' @param time_varying Integer. \code{0} = None; \code{1} = Time varying
#' movement rates.
#' @param time_process Integer. \code{0} = None; \code{1} = Negative
#' binomial (NB1).
#' @param cycle_length Integer. Cycle length or \code{0} for no cycle.
#' @param block_length Integer. Block length for movement rate.
#' @param fish_rate_by Integer. Estimate F \code{0} = As a single value;
#'   \code{1} = By areas; \code{2} = By time steps as blocks;
#'   \code{3} = By both areas and blocks.
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
                    span_liberty = c(1, Inf),
                    time_varying = 0,
                    time_process = c(none = 0, rw = 1),
                    cycle_length = 0,
                    block_length = 0,
                    fish_rate_by = c(none = 0, area = 1, block = 2, both = 3),
                    results_step = 1,
                    nlminb_loops = 5,
                    newton_steps = 5,
                    openmp_cores = 1) {

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
