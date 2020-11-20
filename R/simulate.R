#' Simulate Tag Recovery Data
#'
#' @param data [list()] Input data. See details.
#' @param parameters [list()] Optional parameters to be treated as data.
#'   See details.
#' @param settings [list()] Values that help define the simulation model.
#'   See \code{mmmSet()}.
#'
#' @details TBD
#'
#' @return An object of class \code{mmmSim}.
#' @export
#'
#' @examples
#' sim_x <- create_release_matrix()
#' sim_f <- as.matrix(0.04)
#' sim_m <- 0.1
#' d <- list(x = sim_x, z = sim_z, r = sim_r, f = sim_f, m = sim_m)
#' s1 <- mmmSim(d)
#'
mmmSim <- function(data,
                   parameters = NULL,
                   settings = mmmSet()) {

  # Start the clock ------------------------------------------------------------

  tictoc::tic("mmmSim")

  # Check data arguments -------------------------------------------------------

  cat("checking data arguments \n")
  checkmate::assert_list(data, null.ok = FALSE)
  checkmate::assert_matrix(data$x, mode = "integerish", null.ok = FALSE)
  checkmate::assert_matrix(data$z, mode = "integerish", null.ok = TRUE)
  checkmate::assert_matrix(data$l, mode = "double", null.ok = TRUE)
  checkmate::assert_matrix(data$w, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(data$x, lower = 0, null.ok = FALSE)
  checkmate::assert_numeric(data$z, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_numeric(data$l, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_numeric(data$w, lower = 0, upper = 0, null.ok = TRUE)
  # Optionally parameters
  checkmate::assert_matrix(data$f, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(data$f, lower = 0, null.ok = TRUE)
  checkmate::assert_double(data$m, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(data$h, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(data$u, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(data$u, upper = 1, null.ok = TRUE)
  # Also optionally parameters
  checkmate::assert_array(data$p, mode = "double", null.ok = TRUE)
  checkmate::assert_array(data$r, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(data$b, lower = 0, null.ok = TRUE)

  # Check parameters argument --------------------------------------------------

  cat("checking parameter arguments \n")
  checkmate::assert_list(parameters, null.ok = TRUE)
  # Optionally data
  checkmate::assert_array(parameters$p, mode = "double", null.ok = TRUE)
  checkmate::assert_array(parameters$r, mode = "double", null.ok = TRUE)
  checkmate::assert_matrix(parameters$f, mode = "double", null.ok = TRUE)
  checkmate::assert_numeric(parameters$f, lower = 0, null.ok = TRUE)
  checkmate::assert_double(parameters$m, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(parameters$h, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(parameters$u, lower = 0, len = 1, null.ok = TRUE)
  checkmate::assert_double(parameters$u, upper = 1, null.ok = TRUE)
  checkmate::assert_numeric(parameters$b, lower = 0, null.ok = TRUE)

  # Check settings argument ----------------------------------------------------

  cat("checking settings arguments \n")
  # Check settings
  checkmate::assert_list(settings, null.ok = FALSE)
  # Check error family setting
  checkmate::assert_integerish(settings$error_family, len = 1)
  checkmate::assert_integerish(settings$error_family, lower = 0, upper = 1)
  checkmate::assert_integerish(settings$error_family, any.missing = FALSE)
  # Check min liberty setting
  checkmate::assert_integerish(settings$min_liberty, len = 1)
  checkmate::assert_integerish(settings$min_liberty, any.missing = FALSE)
  checkmate::assert_integerish(settings$min_liberty, lower = 0)
  # Check max liberty setting
  checkmate::assert_integerish(settings$max_liberty, len = 1)
  checkmate::assert_integerish(settings$max_liberty, any.missing = FALSE)
  checkmate::assert_integerish(settings$max_liberty, lower = 0)
  checkmate::assert_true(settings$min_liberty < settings$max_liberty)
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

  # Assembly values ------------------------------------------------------------

  x <- data$x

  # Check required data --------------------------------------------------------

  # Tag releases
  checkmate::assert_matrix(x, mode = "integerish", any.missing = FALSE)
  checkmate::assert_matrix(x, ncols = 4, null.ok = FALSE)
  checkmate::assert_numeric(x, lower = 0, null.ok = FALSE)
  checkmate::assert_true(colnames(x)[1] == "release_step")
  checkmate::assert_true(colnames(x)[2] == "release_area")
  checkmate::assert_true(colnames(x)[3] == "group")
  checkmate::assert_true(colnames(x)[4] == "count")
  # Check values
  if (min(x[, "release_step"]) > 0) cat("caution: time step not indexed from zero")

  # Set index limits -----------------------------------------------------------

  nt <- max(x[, "release_step"]) + 1L # Index from zero for C++
  na <- max(x[, "release_area"]) + 1L # Index from zero for C++
  ng <- max(x[, "group"]) + 1L # Index from zero for C++

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

  # Define index matrix --------------------------------------------------------

  z <- data$z
  diag(z) <- 0

  # Check index matrix ---------------------------------------------------------

  checkmate::assert_matrix(z, mode = "integerish", null.ok = FALSE)
  checkmate::assert_matrix(z, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(z, lower = 0, upper = 1, null.ok = FALSE)
  checkmate::assert_true(nrow(z) == ncol(z))

  # Set index limit ------------------------------------------------------------

  np <- as.integer(sum(z))

  # Set parameter index values -------------------------------------------------

  if (time_varying) {
    if (!block_length && !cycle_length) {
      npt <- nt
      vpt <- c(seq_len(nt)) # Index from one for R
    } else if (block_length && !cycle_length) {
      npt <- ceiling(nt / block_length)
      vpt <- rep(seq_len(npt), each = block_length)[seq_len(nt)]
    } else if (cycle_length && !block_length) {
      npt <- cycle_length
      vpt <- rep(seq_len(npt), ceiling(nt / npt))[seq_len(nt)]
    } else {
      npt <- cycle_length
      vpt_cycle <- rep(seq_len(npt), ceiling(nt / npt))[seq_len(nt)]
      vpt <- rep(vpt_cycle, each = block_length)[seq_len(nt)]
    }
  } else {
    npt <- 1L
    vpt <- rep(1L, nt) # Index from one for R
  }

  # Create fishing rate indexes ------------------------------------------------

  if (!is.null(data$f)) {
    nft <- nrow(data$f)
    nfa <- ncol(data$f)
    vft <- rep(seq_len(nft), each = ceiling(nt / nft))[seq_len(nt)]
    vfa <- rep(seq_len(nfa), each = ceiling(na / nfa))[seq_len(na)]
    if (nt %% nft) {
      cat("warning: nft does not divide nt evenly \n")
    }
    if (!is.element(nfa, c(1, na))) {
      cat("warning: nfa must equal 1 or na \n")
    }
    cat("using f from data and ignoring fish_rate_by in settings \n")
  } else if (!is.null(parameters$f)) {
    nft <- nrow(parameters$f)
    nfa <- ncol(parameters$f)
    vft <- rep(seq_len(nft), each = ceiling(nt / nft))[seq_len(nt)]
    vfa <- rep(seq_len(nfa), each = ceiling(na / nfa))[seq_len(na)]
    if (nt %% nft) {
      cat("warning: nft does not divide nt evenly \n")
    }
    if (!is.element(nfa, c(1, na))) {
      cat("warning: nfa must equal 1 or na \n")
    }
    cat("initial f from parameters; ignoring fish_rate_by in settings \n")
  } else {
    stop("missing fishing rate matrix f")
  }

  # Unpack arguments -----------------------------------------------------------

  # Tag reporting rate
  if (!is.null(data$l)) {
    l <- data$l
  } else {
    l <- matrix(1, nrow = 1L, ncol = 1L)
  }
  # Fishing rate weighting
  if (!is.null(data$w)) {
    w <- data$w
  } else {
    w <- matrix(1, nrow = 1L, ncol = 1L)
  }

  # Check arguments ------------------------------------------------------------

  # Tag reporting rate
  checkmate::assert_matrix(l, mode = "double", null.ok = FALSE)
  checkmate::assert_matrix(l, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(l, lower = 0, null.ok = FALSE)
  checkmate::assert_true(ncol(l) == 1 || ncol(l) == na)
  # Fishing rate weighting
  checkmate::assert_matrix(w, mode = "double", null.ok = FALSE)
  checkmate::assert_matrix(w, any.missing = FALSE, null.ok = FALSE)
  checkmate::assert_numeric(w, lower = 0, null.ok = FALSE)
  checkmate::assert_true(ncol(w) == 1 || ncol(w) == na)
  checkmate::assert_true(all(colSums(w) == 1))

  # Create tag reporting rate indexes ------------------------------------------

  # Set for time step
  if (nrow(l) == 1L) {
    nlt <- 1L
    vlt <- rep(1L, nt) # Index from one for R
  } else {
    nlt <- nt
    vlt <- c(seq_len(nt))
  }
  # Set for areas
  if (ncol(l) == 1L) {
    nla <- 1L
    vla <- rep(1L, na)
  } else {
    nla <- na
    vla <- c(seq_len(na))
  }

  # Create fishing rate weighting indexes --------------------------------------

  # Set for time step
  if (nrow(w) == 1L) {
    nwt <- 1L
    vwt <- rep(1L, nt) # Index from one for R
  } else {
    nwt <- nrow(w)
    vwt <- rep(c(seq_len(nwt)), ceiling(nt / nwt))[seq_len(nt)]
  }
  # Set for areas
  if (ncol(w) == 1L) {
    nwa <- 1L
    vwa <- rep(1L, na)
  } else {
    nwa <- na
    vwa <- c(seq_len(na))
  }

  # Unpack arguments -----------------------------------------------------------

  # Array movement parameters
  if (!is.null(data$r)) {
    r <- data$r
    p <- create_movement_parameters(r, z = z)
  } else if (!is.null(data$p)) {
    p <- data$p
    r <- create_movement_rates(p, z = z)
  } else if (!is.null(parameters$r)) {
    r <- parameters$r
    p <- create_movement_parameters(r, z = z)
  } else if (!is.null(parameters$p)) {
    p <- parameters$p
    r <- create_movement_rates(parameters$p, z = z)
  } else {
    stop("data or parameter must contain array r or array p")
  }
  # Matrix fishing mortality rate
  if (!is.null(data$f)) {
    f <- data$f
  } else if (!is.null(parameters$f)) {
    f <- parameters$f
  } else {
    stop("data or parameter must contain matrix f")
  }
  # Scalar natural mortality
  if (!is.null(data$m)) {
    m <- data$m
  } else if (!is.null(parameters$m)) {
    m <- parameters$m
  } else {
    stop("data or parameter must contain scalar m")
  }
  # Scalar tag loss rate
  if (!is.null(data$h)) {
    h <- data$h
  } else if (!is.null(parameters$h)) {
    h <- parameters$h
  } else {
    h <- 0L
  }
  # Scalar initial tag loss rate (proportion)
  if (!is.null(data$u)) {
    u <- data$u
  } else if (!is.null(parameters$u)) {
    u <- parameters$u
  } else {
    u <- 0
  }
  # Vector fishing bias
  if (!is.null(data$b)) {
    b <- data$b
  } else if (!is.null(parameters$b)) {
    b <- parameters$b
  } else {
    b <- rep(1, ng)
  }
  # Scalar negative binomial dispersion
  if (!is.null(data$k)) {
    k <- data$k
  } else if (!is.null(parameters$k)) {
    k <- parameters$k
  } else {
    k <- 1L
  }

  # Check arguments ------------------------------------------------------------

  # Movement rate array
  checkmate::assert_array(r, mode = "double", any.missing = FALSE, d = 4)
  checkmate::assert_numeric(r, lower = 0, upper = 1, len = na * na * npt * ng)
  checkmate::assert_true(dim(r)[1] == na, na.ok = FALSE)
  checkmate::assert_true(dim(r)[2] == na, na.ok = FALSE)
  checkmate::assert_true(dim(r)[3] == npt, na.ok = FALSE)
  checkmate::assert_true(dim(r)[4] == ng, na.ok = FALSE)
  # Fishing rate matrix
  checkmate::assert_matrix(f, mode = "double", any.missing = FALSE)
  checkmate::assert_numeric(f, lower = 0, null.ok = FALSE)
  checkmate::assert_true(dim(f)[1] == nft, na.ok = FALSE)
  checkmate::assert_true(dim(f)[2] == nfa, na.ok = FALSE)

  # Compute transposes ---------------------------------------------------------

  # Matrices
  tx <- t(x)
  tz <- t(z)
  tl <- t(l)
  tw <- t(w)
  tf <- t(f)
  # Other transformations
  d <- (1 - u)

  # Initialize arrays ----------------------------------------------------------

  N <- array(0, dim = c(na, nt, ng, na, nt))
  S <- array(0, dim = c(na, nt, ng))
  Y <- array(0, dim = c(na, nt, ng, na, nt))

  # Populate N ----------------------------------------------------------------

  for (i in seq_len(ncol(tx))) {
    N[tx[2, i] + 1L,
       tx[1, i] + 1L,
       tx[3, i] + 1L,
       tx[2, i] + 1L,
       tx[1, i] + 1L] = d * tx[4, i]
  }

  # Compute S -----------------------------------------------------------------

  for (mg in seq_len(ng)) {
    for (ct in seq_len(nt)) {
      for (ca in seq_len(na)) {
        S[ca, ct, mg] <- exp(-(b[mg] * tf[vfa[ca], vft[ct]] *
                                  tw[vwa[ca], vwt[ct]] + m + h))
      }
    }
  }

  # Compute Y -----------------------------------------------------------------

  cat("only poisson error family implemented")
  for (mt in seq_len(nt)) {
    for (ma in seq_len(na)) {
      for (mg in seq_len(ng)) {
        # Were tags released in this stratum?
        if (N[ma, mt, mg, ma, mt] > 0) {
          # Project the abundance array forward
          if (mt < nt) {
            for (ct in c((mt + 1):nt)) {
              for (ca in seq_len(na)) {
                for (pa in seq_len(na)) {
                  N[ca, ct, mg, ma, mt] <- N[ca, ct, mg, ma, mt] +
                    N[pa, ct - 1, mg, ma, mt] * S[pa, ct - 1, mg] *
                    r[pa, ca, vpt[ct], mg]
                }
              }
            }
          }
          # Compute rt_min and rt_max from min and max liberty
          rt_min <- mt + min_liberty
          rt_max <- mt + max_liberty + 1L
          if (rt_max > nt) rt_max <- nt
          if (rt_min > rt_max) rt_min <- rt_max
          # Compute predicted recoveries
          for (rt in c(rt_min:rt_max)) {
            for (ra in seq_len(na)) {
              # Compute the recovery array
              Y_elem <- N[ra, rt, mg, ma, mt] *
                (1L - exp(-b[mg] * tf[vfa[ra], vft[rt]] *
                            tw[vwa[ra], vwt[rt]])) *
                tl[vla[ra], vlt[rt]]
                Y_elem <- rpois(1L, Y_elem)
              Y[ra, rt, mg, ma, mt] <- round(Y_elem)
            }
          }
        } # End if()
      }
    }
  }

  # Compute y -----------------------------------------------------------------

  # Initialize
  y <- matrix(0L, nrow = sum(as.logical(Y)), ncol = 6)
  release_names <- c("release_step", "release_area", "group")
  recover_names <- c("recover_step", "recover_area", "count")
  colnames(y) <- c(release_names, recover_names)
  # Populate
  ind <- 1
  for (ra in seq_len(na)) {
    for (rt in seq_len(nt)) {
      for (mg in seq_len(ng)) {
        for (ma in seq_len(na)) {
          for (mt in seq_len(nt)) {
            if (Y[ra, rt, mg, ma, mt] > 0L) {
              Y_elem <- Y[ra, rt, mg, ma, mt]
              y[ind, ] <- c(c(mt, ma, mg, rt, ra) - 1L, Y_elem)
              ind <- ind + 1
            }
          }
        }
      }
    }
  }

  # Assemble lists -------------------------------------------------------------

  # Simulation list
  simulation_list <- list(
    p = p,
    r = r,
    x = x,
    y = y,
    z = z,
    l = l,
    w = w,
    f = f,
    m = m,
    h = h,
    u = u,
    b = b,
    k = k
  )
  # Settings
  settings_list <- settings

  # Return list ----------------------------------------------------------------

  cat("\nreturning mmmSim object\n")
  structure(list(
    simulation  = simulation_list,
    settings    = settings_list),
    class       = c("mmmSim"))
}

#' Create Tag Release Matrix
#'
#' @param n [integer()] Mean tag release count per stratum
#' @param steps [integer()] Number of time steps
#' @param areas [integer()] Number of areas
#' @param groups [integer()] Number of groups
#' @param errors [logical()] Poisson distributed tag release counts?
#'
#' @return An integer matrix
#' @export
#'
#' @examples
#'
#' m1 <- create_release_matrix()
#' m2 <- create_release_matrix(groups = 3)
#' m3 <- create_release_matrix(errors = TRUE)
#'
create_release_matrix <- function (n = 1000,
                                   steps = 30,
                                   areas = 3L,
                                   groups = 1L,
                                   errors = FALSE) {

  # Check arguments
  checkmate::assert_integerish(n, lower = 0, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(steps, lower = 10, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(areas, lower = 2, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(groups, lower = 1, len = 1, any.missing = FALSE)
  checkmate::assert_logical(errors, len = 1, any.missing = FALSE)
  # Create vectors
  t <- rep(c(seq_len(steps) - 1L), each = areas * groups)
  a <- rep(rep(c(seq_len(areas) - 1L), each = groups), steps)
  g <- rep(c(seq_len(groups) - 1L), steps * areas)
  # Create counts
  if (errors) {
    count <- rpois(n = steps * areas * groups, lambda = n)
  } else {
    count <- rep(n, steps * areas * groups)
  }
  # Create matrix
  x <- as.matrix(data.frame(
    release_step = t,
    release_area = a,
    group = g,
    count = count))
  mode(x) <- c("integer")
  # Return
  return(x)
}
