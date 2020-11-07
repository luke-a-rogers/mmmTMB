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

