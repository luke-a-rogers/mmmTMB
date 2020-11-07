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
  m <- as.matrix(data.frame(
    release_step = t,
    release_area = a,
    group = g,
    count = count))
  mode(m) <- c("integer")
  # Return
  return(m)
}

#' Create Movement Index Matrix
#'
#' @param areas [integer()] Number of areas
#' @param pattern [integer()] One of \code{0}: full movement, or \code{1}:
#'   immediate neighbours only
#' @param allow [matrix()] Two columns. Each row gives a pair of area
#'   indexes specifying allowed directional movement from the first
#'   column to the second from one time step to the next. Indexed from 1.
#' @disallow [matrix()] Same as \code{allow}, but disallows movement
#'   between specified pairs.
#'
#' @return A square binary matrix with zero diagonal
#' @export
#'
#' @examples
#'
#' # Neighbours only (default)
#' m1 <- create_movement_index(6)
#'
#' # Full movement
#' m2 <- create_movement_index(6, 0)
#'
#' # Neighbours plus 1-6, 2-5, and 5-2, but not 6-1
#' ind <- matrix(c(1, 6, 2, 5, 5, 2), ncol = 2, byrow = TRUE)
#' m3 <- create_movement_index(6, 1, ind)
#'
#' # Full minus 1-6
#' dind <- matrix(c(1, 6), ncol = 2, byrow = TRUE)
#' m4 <- create_movement_index(6, 0, disallow = dind)
#'
create_movement_index <- function (areas,
                                   pattern = NULL,
                                   allow = NULL,
                                   disallow = NULL) {

  # Check arguments
  checkmate::assert_integerish(areas, lower = 2, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(pattern, len = 1, null.ok = TRUE)
  checkmate::assert_integerish(pattern, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_integerish(pattern, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_matrix(allow, mode = "integerish", null.ok = TRUE)
  checkmate::assert_matrix(allow, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_matrix(allow, ncols = 2, null.ok = TRUE)
  # Initialize matrix
  m <- matrix(0L, nrow = areas, ncol = areas)
  # Default
  if (is.null(pattern) & is.null(allow)) {
    # Immediate neighbours only
    for (i in seq_len(areas - 1L)) {
      m[i, i + 1L] <- 1L
      m[i + 1L, i] <- 1L
    }
  }
  # Add pattern
  if (!is.null(pattern)) {
    if (pattern) {
      for (i in seq_len(areas - 1L)) {
        m[i, i + 1L] <- 1L
        m[i + 1L, i] <- 1L
      }
    } else {
      m <- matrix(1L, nrow = areas, ncol = areas)
      diag(m) <- 0L
    }
  }
  # Allow indexes
  if (!is.null(allow)) {
    m[allow] <- 1L
    diag(m) <- 0L
  }
  # Disallow indexes
  if (!is.null(disallow)) {
    m[disallow] <- 0L
    diag(m) <- 0L
  }
  # Return
  return(m)
}
