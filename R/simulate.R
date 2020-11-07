#' Simulate Tag Recovery Data
#'
#' @param data A list.
#' @param parameters A list.
#' @param settings A list. See \code{mmmSet()}.
#'
#' @return An object of class \code{mmmSim}.
#' @export
#'
#' @examples
#'
mmmSim <- function(data,
                   parameters,
                   settings = mmmSet()) {

  # Start the clock ------------------------------------------------------------

  tictoc::tic("mmmSim")

  # Check data arguments -------------------------------------------------------

  cat("checking data arguments \n")
  checkmate::assert_list(data, null.ok = FALSE)
  checkmate::assert_matrix(data$mT, mode = "integerish", null.ok = FALSE)
  checkmate::assert_matrix(data$mR, mode = "integerish", null.ok = TRUE)
  checkmate::assert_matrix(data$mI, mode = "integerish", null.ok = TRUE)
  checkmate::assert_matrix(data$mL, mode = "double", null.ok = TRUE)
  checkmate::assert_matrix(data$mW, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(data$mT, lower = 0, null.ok = FALSE)
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
  checkmate::assert_array(parameters$aK, mode = "double", null.ok = TRUE)
  checkmate::assert_matrix(parameters$mF, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(parameters$mF, lower = 0, null.ok = TRUE)
  checkmate::assert_double(parameters$sM, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(parameters$sH, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(parameters$sC, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(parameters$sC, upper = 1, null.ok = TRUE)
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

  # Assembly values ------------------------------------------------------------

  mT <- data$mT

  # Check required data --------------------------------------------------------

  # Tag releases
  checkmate::assert_matrix(mT, mode = "integerish", any.missing = FALSE)
  checkmate::assert_matrix(mT, ncols = 4, null.ok = FALSE)
  checkmate::assert_numeric(mT, lower = 0, null.ok = FALSE)
  checkmate::assert_true(colnames(mT)[1] == "release_step")
  checkmate::assert_true(colnames(mT)[2] == "release_area")
  checkmate::assert_true(colnames(mT)[3] == "group")
  checkmate::assert_true(colnames(mT)[4] == "count")
  # Check values
  if (min(mT$release_step) > 0) cat("caution: time step not indexed from zero")

  # Set index limits -----------------------------------------------------------

  nt <- max(mT$release_step) + 1L # Index from zero for C++
  na <- max(mT$release_area) + 1L # Index from zero for C++
  ng <- max(mT$group) + 1L # Index from zero for C++

  # Index matrix
  mI <- data$mI
  diag(mI) <- 0

  # Check index matrix ---------------------------------------------------------

  checkmate::assert_matrix(mI, mode = "integerish", null.ok = FALSE)
  checkmate::assert_matrix(mI, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(mI, lower = 0, upper = 1, null.ok = FALSE)
  checkmate::assert_true(nrow(mI) == ncol(mI))

  # Set index limit
  np <- as.integer(sum(mI))

  # Set parameter index values -------------------------------------------------

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

  # Create fishing rate indexes ------------------------------------------------

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
    stop("missing fishing rate matrix mF")
  }

  # Unpack arguments -----------------------------------------------------------

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

  # Check arguments ------------------------------------------------------------

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

  # Create tag reporting rate indexes ------------------------------------------

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

  # Create fishing rate weighting indexes --------------------------------------

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

  # Unpack arguments -----------------------------------------------------------

  # Array movement parameters
  if (!is.null(parameters$aK)) {
    aK <- parameters$aK
  } else if (!is.null(parameters$aP)) {
    aK <- create_movement_rates(aP, m = mI)
  } else {
    stop("data or parameter must contain array aP")
  }
  # Matrix fishing mortality rate
  if (!is.null(data$mF)) {
    mF <- data$mF
  } else if (!is.null(parameters$mF)) {
    mF <- parameters$mF
  } else {
    stop("data or parameter must contain matrix mF")
  }
  # Scalar natural mortality
  if (!is.null(data$sM)) {
    sM <- data$sM
  } else if (!is.null(parameters$sM)) {
    sM <- parameters$sM
  } else {
    stop("data or parameter must contain scalar sM")
  }
  # Scalar tag loss rate
  if (!is.null(data$sH)) {
    sH <- data$sH
  } else if (!is.null(parameters$sH)) {
    sH <- parameters$sH
  } else {
    sH <- 0L
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
  # Scalar negative binomial dispersion
  if (!is.null(parameters$sD)) {
    sD <- parameters$sD
  } else {
    sD <- 1L
  }
  # Check arguments ------------------------------------------------------------

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

  # Compute transposes ---------------------------------------------------------

  # Matrices
  tmT <- t(mT)
  tmI <- t(mI)
  tmL <- t(mL)
  tmW <- t(mW)
  tmF <- t(mF)
  # Other transformations
  sA <- (1 - sC)

  # Populate aN ----------------------------------------------------------------

  for (i in seq_len(ncol(tmT))) {
    aN[tmT[2, i], tmT[1, i], tmT[3, i], tmT[2, i], tmT[1, i]] = sA * tmT[4, i]
  }

  # Compute aS -----------------------------------------------------------------

  for (mg in seq_len(ng)) {
    for (ct in seq_len(nt)) {
      for (ca in seq_len(na)) {
        aS[ca, ct, mg] <- exp(-(vB[mg] * tmF[vfa[ca], vft[ct]] *
                                  tmW[vwa[ca], vwt[ct]] + sM + sH))
      }
    }
  }

  # Compute aR -----------------------------------------------------------------

  for (mt in seq_len(nt)) {
    for (ma in seq_len(na)) {
      for (mg in seq_len(ng)) {
        # Were tags released in this stratum?
        if (aN[ma, mt, mg, ma, mt] > 0) {
          # Project the abundance array forward
          if (mt < nt) {
            for (ct in c((mt + 1):nt)) {
              for (ca in seq_len(na)) {
                for (pa in seq_len(na)) {
                  aN[ca, ct, mg, ma, mt] <- aN[ca, ct, mg, ma, mt] +
                    aN[pa, ct - 1, mg, ma, mt] * aS[pa, ct - 1, mg] *
                    aK[pa, ca, vpt[ct], mg]
                }
              }
            }
          }
          # Compute predicted recoveries
          for (rt in c(mt:nt)) {
            # Was the recapture time within the span at liberty?
            if (rt > (mt + span_liberty[1] - 1L)) {
              if (rt < (mt + span_liberty[2] + 1L)) {
                for (ra in seq_len(na)) {
                  # Compute the recovery array
                  sR <- aN[ra, rt, mg, ma, mt] *
                    (1L - exp(-vB[mg] * tmF[vfa[ra], vft[rt]] *
                                     tmW[vwa[ra], vwt[rt]])) *
                    tmL[vla[ra], vlt[rt]]
                  if (error_family) {
                    if (sR > 0) {
                      sR <- rnbinom(1L, mu = sR, size = sD / sR)
                    } else {
                      sR <- 0L
                    }
                  } else {
                    sR <- rpois(1L, sR)
                  }
                  aR[ra, rt, mg, ma, mt] <- round(sR)
                }
              }
            }
          }
        } # End if()
      }
    }
  }

  # Compute mR -----------------------------------------------------------------

  # Initialize
  mR <- matrix(0L, nrow = sum(as.logical(aR)), ncol = 6)
  release_names <- c("release_step", "release_area", "group_index")
  recover_names <- c("recover_step", "recover_area", "count")
  colnames(mR) <- c(release_names, recover_names)
  # Populate
  ind <- 1
  for (ra in seq_len(na)) {
    for (rt in seq_len(nt)) {
      for (mg in seq_len(ng)) {
        for (ma in seq_len(na)) {
          for (mt in seq_len(nt)) {
            if (aR[ra, rt, mg, ma, mt] > 0L) {
              sR <- aR[ra, rt, mg, ma, mt]
              mR[ind, ] <- c(c(mt, ma, mg, rt, ra) - 1L, sR)
              ind <- ind + 1
            }
          }
        }
      }
    }
  }

  # Assemble lists -------------------------------------------------------------

  # Data
  data_list <- list(
    mT = mT,
    mR = mR,
    mI = mI,
    mL = mL,
    mW = mW
  )
  # Parameters
  parameter_list <- list(
    aP = aP,
    aK = aK,
    mF = mF,
    sM = sM,
    sH = sH,
    sC = sC,
    vB = vB,
    sD = sD
  )
  # Settings
  settings_list <- list()


  # Return list ----------------------------------------------------------------

  cat("\nreturning mmmSim object\n")
  structure(list(
    data        = data_list,
    parameters  = parameters_list,
    settings    = settings_list),
    class       = c("mmmSim"))
}
