#' Simulate Tag Recovery Data
#'
#' @param data [list()] Input data. See details.
#' @param parameters [list()] Parameters to be treated as data.
#'   See details.
#' @param settings [list()] Values that help define the simulation model.
#'   See \code{mmmSet()}.
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
#'   \item \code{f} A matrix of fishing mortality rates. Must be included
#'     in \code{data} or \code{parameters}. See [mmmRates()].
#'   \item \code{m} The scalar instantaneous natural mortality rate.
#'     Must be included in \code{data} or \code{parameters}.
#' }
#' The [list()] argument \code{parameters} must contain:
#' \itemize{
#'   \item \code{p} An array of movement parameters. See
#'     [create_movement_parameters()].
#' }
#' The list argument \code{parameters} may optionally contain:
#' \itemize{
#'   \item \code{f} A matrix of fishing mortality rates. Must be included
#'     in \code{data} or \code{parameters}. See [mmmRates()].
#'   \item \code{m} The scalar instantaneous natural mortality rate.
#'     Must be included in \code{data} or \code{parameters}.
#' }
#'
#' @return An object of class \code{mmmSim}.
#' @export
#'
#' @examples
#' # Data
#' sim_x <- create_release_matrix()
#' sim_z <- mmmIndex(3, pattern = 1)
#' sim_f <- as.matrix(0.04)
#' sim_m <- 0.1
#' sim_h <- 0.02
#' sim_u <- 0.1
#' # Parameters
#' v <- c(0.9, 0.1, 0.0, 0.1, 0.8, 0.3, 0.0, 0.1, 0.7)
#' sim_r <- array(v, dim = c(3,3,1,1))
#' sim_p <-create_movement_parameters(sim_r, sim_z)
#' # Simulation
#' data <- list(
#'   x = sim_x,
#'   z = sim_z,
#'   f = sim_f,
#'   m = sim_m,
#'   h = sim_h,
#'   u = sim_u)
#' parameters <- list(p = sim_p)
#' s1 <- mmmSim(data, parameters)
#'
mmmSim <- function(data,
                   parameters = NULL,
                   settings = mmmSet()) {

  # Start the clock ------------------------------------------------------------

  tictoc::tic("mmmSim")

  # Insert dummy element -------------------------------------------------------

  data$y <- matrix(0, nrow = 1, ncol = 6)
  rel_cols <- c("release_step", "release_area", "group")
  rec_cols <- c("recover_step", "recover_area", "count")
  colnames(data$y) <- c(rel_cols, rec_cols)

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

  # Set index limits -----------------------------------------------------------

  nt <- max(x[, "release_step"]) + 1L # Index from zero for C++
  na <- max(x[, "release_area"]) + 1L # Index from zero for C++
  ng <- max(x[, "group"]) + 1L # Index from zero for C++
  np <- as.integer(sum(z))

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

  # Create tag reporting rate indexes ------------------------------------------

  nlt <- nrow(l)
  nla <- ncol(l)
  vlt <- rep(c(seq_len(nlt)), ceiling(nt / nlt))[seq_len(nt)]
  vla <- rep(c(seq_len(nla)), ceiling(na / nla))[seq_len(na)]

  # Create fishing rate weighting indexes --------------------------------------

  nwt <- nrow(w)
  nwa <- ncol(w)
  vwt <- rep(c(seq_len(nwt)), ceiling(nt / nwt))[seq_len(nt)]
  vwa <- rep(c(seq_len(nwa)), ceiling(na / nwa))[seq_len(na)]

  # Create fishing rate indexes ------------------------------------------------

  if (!is.null(f)) nft <- nrow(f) else nft <- 1L
  if (!is.null(f)) nfa <- ncol(f) else nfa <- 1L
  vft <- rep(seq_len(nft), each = ceiling(nt / nft))[seq_len(nt)]
  vfa <- rep(seq_len(nfa), each = ceiling(na / nfa))[seq_len(na)]
  if (nt %% nft) cat("caution: nft does not divide nt evenly\n")
  if (!is.element(nfa, c(1, na))) stop("nfa must equal 1 or na\n")

  # Create movement parameter indexes ------------------------------------------

  if (time_varying) {
    if (!block_length && !cycle_length) {
      npt <- nt
      vpt <- c(seq_len(nt)) # Index from zero for TMB
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
    vpt <- rep(1L, nt) # Index from zero for TMB
  }

  # Confirm arguments ----------------------------------------------------------

  if (is.null(p)) stop("data or parameters must contain p\n")
  if (is.null(f)) stop("data or parameters must contain f\n")
  if (is.null(m)) stop("data or parameters must contain m\n")
  if (is.null(b)) b <- rep(1, ng)

  # Create rates ---------------------------------------------------------------

  r <- create_movement_rates(p, z)

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
        S[ca, ct, mg] <-
          exp(-b[mg] * tf[vfa[ca], vft[ct]] * tw[vwa[ca], vwt[ct]] - m - h)
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
    u = u
    # b = b,
    # k = k
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
