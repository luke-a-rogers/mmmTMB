

mmmSim <- function(data,
                   parameters) {


  # Check arguments ------------------------------------------------------------



  # Assembly values ------------------------------------------------------------

  # Tag releases
  mT <- data$mT
  # Set index limits
  nt <- max(mT$release_step) + 1L # Index from zero for C++
  na <- max(mT$release_area) + 1L # Index from zero for C++
  ng <- max(mT$group) + 1L # Index from zero for C++
  # Index matrix
  mI <- data$mI
  diag(mI) <- 0
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


  # In: aP, sC
  # Out: aK, sA,



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
    mT = mR,
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
