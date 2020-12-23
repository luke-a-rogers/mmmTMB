#' Prepare Tag Data
#'
#' @description Prepare tag release and recovery data prior to fitting
#' a Markov movement model via \code{mmmFit()}.
#'
#' @param x [data.frame()] Tag release data. See details.
#' @param y [data.frame()] Tag recovery data. See details.
#' @param cols [list()] Named list of character strings. See details.
#' @param groups [list()] Named list of atomic vectors. See details.
#' @param step [character()] One of \code{year}, \code{month}, or
#'   \code{day}.
#' @param limits [character()] Start and end dates as `"%Y-%M-%D"`
#'   character strings.
#' @param days_liberty [integer()] Minimum number of days at liberty.
#'
#' @details
#'
#' The \code{x} & \code{y} arguments must be data frames that contain
#' certain columns that have the same column name in both data frames.
#' Both \code{x} & \code{y} must include:
#' \itemize{
#'   \item [character()] column giving unique individual tag IDs.
#'   \item [Date()] column giving individual tag release dates.
#'   \item [integer()] column giving release areas indexed from one.
#' }
#' In addition, \code{y} must include:
#' \itemize{
#'   \item [Date()] column giving individual tag recovery dates.
#'   \item [integer()] column giving recovery areas indexed from one.
#' }
#' Optionally, both \code{x} & \code{y} can include:
#' \itemize{
#'   \item [atomic()] column giving values grouped by \code{groups}.
#' }
#' Note that [character()] columns in `"%Y-%M-%D"` format may be substituted
#' for [Date()] columns.
#'
#' The \code{cols} argument is a [list()] of named elements. The elements
#'   must be [character()] string names of columns in \code{x} or
#'   \code{y}. There must be one element for each of the required element
#'   names (except as noted):
#'   \itemize{
#'     \item \code{id}: \code{x} & \code{y} Individual tag identification
#'     \item \code{release_date}: \code{x} & \code{y} Release date
#'     \item \code{release_area}: \code{x} & \code{y} Release area
#'     \item \code{recover_date}: \code{y} Recovery date
#'     \item \code{recover_area}: \code{y} Recovery area
#'     \item \code{group}: (Optional) Grouping variable
#'   }
#'
#' The \code{groups} argument is a [list()] of named elements. The elements
#'   must be atomic vectors that include all values included in each named
#'   group. All element vectors must be of the same type. The names of the
#'   element vectors must be unique. For \code{groups} to have any effect,
#'   the element \code{group} in argument \code{cols} must name a column
#'   in both \code{x} & \code{y} that has values include in the element
#'   vectors in \code{groups}.
#'
#' @importFrom magrittr `%>%`
#'
#' @return A list object of class \code{mmmTags}.
#'
#' @export
#'
#' @examples
#'
#' x <- data.frame(
#'   Release_Date = c("2010-06-01", "2015-04-01", "2018-06-01", "2019-06-01"),
#'   Release_Area = c(1,1,1,1),
#'   Release_Group = c("F", "F", "F", "M"),
#'   ID = c("TAG001", "TAG002", "TAG003", "TAG004"))
#' y <- data.frame(
#'   Release_Date = c("2010-06-01", "2018-06-01", "2019-06-01"),
#'   Release_Area = c(1,1,1),
#'   Release_Group = c("F", "F", "M"),
#'   Recover_Date = c("2013-08-01", "2020-01-01", "2020-08-01"),
#'   Recover_Area = c(1,2,1),
#'   ID = c("TAG001", "TAG003", "TAG004"))
#' cols <- list(
#'   release_date = "Release_Date",
#'   release_area = "Release_Area",
#'   recover_date = "Recover_Date",
#'   recover_area = "Recover_Area",
#'   group = "Release_Group",
#'   id = "ID")
#' groups <- list(m = "M", f = "F")
#' step <- "year"
#' limits <- c("2010-01-01", "2020-12-31")
#' days_liberty <- 0L
#'
#' tags <- mmmTags(x, y, cols, groups, step, limits, days_liberty)
#'
mmmTags <- function (x,
                     y,
                     cols = NULL,
                     groups = NULL,
                     step = NULL,
                     limits = NULL,
                     days_liberty = NULL) {

  # Start the clock ------------------------------------------------------------

  tictoc::tic("mmmTags()")

  # Check arguments ------------------------------------------------------------

  # Arguments
  checkmate::assert_data_frame(x)
  checkmate::assert_data_frame(y)
  checkmate::assert_list(cols, unique = TRUE)
  checkmate::assert_list(groups, unique = TRUE, null.ok = TRUE)
  checkmate::assert_character(step, len = 1, null.ok = TRUE)
  checkmate::assert_choice(step, c("year", "month", "day"), null.ok = TRUE)
  checkmate::assert_integerish(days_liberty, lower = 0, len = 1, null.ok = TRUE)
  # Data frame x column names
  checkmate::assert_choice(cols$release_date, colnames(x))
  checkmate::assert_choice(cols$release_area, colnames(x))
  checkmate::assert_choice(cols$id, colnames(x))
  # Data frame y column names
  checkmate::assert_choice(cols$release_date, colnames(y))
  checkmate::assert_choice(cols$release_area, colnames(y))
  checkmate::assert_choice(cols$recover_date, colnames(y))
  checkmate::assert_choice(cols$recover_area, colnames(y))
  checkmate::assert_choice(cols$id, colnames(y))

  # Set defaults ---------------------------------------------------------------

  # Groups
  if (is.null(groups) | is.null(cols$group)) {
    groups <- list(nogroups = 1L)
    cols$group <- "nogroup"
    x$nogroup <- 1L
    y$nogroup <- 1L
    cat("no groups: if desired specify in groups and cols\n")
  }
  # Time step
  if (is.null(step)) {
    step <- "year"
    cat("using default time step: year")
  }
  # Days at liberty
  if (is.null(days_liberty)) {
    if (step == "year") {
      days_liberty <- 365L
    } else if (step == "month") {
      days_liberty <- 31L
    } else {
      days_liberty <- 1L
    }
  }

  # Check optional columns -----------------------------------------------------

  checkmate::assert_choice(cols$group, colnames(x))
  checkmate::assert_choice(cols$group, colnames(y))

  # Rename columns -------------------------------------------------------------

  # Rename columns x
  colnames(x)[which(colnames(x) == cols$release_date)] <- "release_date"
  colnames(x)[which(colnames(x) == cols$release_area)] <- "release_area"
  colnames(x)[which(colnames(x) == cols$group)] <- "group_raw"
  colnames(x)[which(colnames(x) == cols$id)] <- "id"
  # Rename columns y
  colnames(y)[which(colnames(y) == cols$release_date)] <- "release_date"
  colnames(y)[which(colnames(y) == cols$release_area)] <- "release_area"
  colnames(y)[which(colnames(y) == cols$group)] <- "group_raw"
  colnames(y)[which(colnames(y) == cols$recover_date)] <- "recover_date"
  colnames(y)[which(colnames(y) == cols$recover_area)] <- "recover_area"
  colnames(y)[which(colnames(y) == cols$id)] <- "id"

  # Define x and y -------------------------------------------------------------

  # Define columns
  release_cols <- c("release_date", "release_area")
  recover_cols <- c("recover_date", "recover_area")
  specify_cols <- c("group_raw", "id")
  colnames_x <- c(release_cols, specify_cols)
  colnames_y <- c(release_cols, recover_cols, specify_cols)
  # Define x and y
  x <- x[, colnames_x]
  y <- y[, colnames_y]
  # Drop NA values
  x <- tidyr::drop_na(x)
  y <- tidyr::drop_na(y)
  # Convert to dates
  x <- dplyr::mutate(x, release_date = lubridate::as_date(release_date))
  y <- dplyr::mutate(y, release_date = lubridate::as_date(release_date))
  y <- dplyr::mutate(y, recover_date = lubridate::as_date(recover_date))

  # Check column class or type -------------------------------------------------

  # Data frame x
  checkmate::assert_date(x$release_date, any.missing = FALSE)
  checkmate::assert_integerish(x$release_area, lower = 1, any.missing = FALSE)
  # checkmate::assert_integerish(x$group_raw, lower = 1, any.missing = FALSE)
  # checkmate::assert_character(x$id, any.missing = FALSE)

  # Data frame y
  checkmate::assert_date(y$release_date, any.missing = FALSE)
  checkmate::assert_date(y$recover_date, any.missing = FALSE)
  checkmate::assert_integerish(y$release_area, lower = 1, any.missing = FALSE)
  checkmate::assert_integerish(y$recover_area, lower = 1, any.missing = FALSE)
  # checkmate::assert_integerish(y$group_raw, lower = 1, any.missing = FALSE)
  # checkmate::assert_character(y$id, any.missing = FALSE)

  # Convert limits to dates ----------------------------------------------------

  if (is.null(limits)) {
    limits <- c(min(x$release_date), max(c(x$release_date, y$recover_date)))
  } else {
    limits <- lubridate::as_date(limits)
  }

  # Filter events by dates -----------------------------------------------------

  # Release
  x <- x %>% dplyr::filter(
    release_date >= limits[1],
    release_date <= limits[2])
  # Recover
  y <- y %>% dplyr::filter(
    release_date >= limits[1],
    release_date <= limits[2],
    recover_date >= limits[1],
    recover_date <= limits[2],
    release_date <= recover_date)

  # Arrange columns ------------------------------------------------------------

  x <- x %>% dplyr::arrange(id, release_date, release_area, .keep_all = TRUE)
  y <- y %>% dplyr::arrange(
    id,
    release_date,
    release_area,
    recover_date,
    recover_area,
    .keep_all = TRUE)

  # Remove duplicate events ----------------------------------------------------

  # Count the rows
  x_rows <- nrow(x)
  y_rows <- nrow(y)
  # Remove duplicates
  x <- dplyr::distinct(x, id, release_date, release_area, .keep_all = TRUE)
  y <- dplyr::distinct(y, id, release_date, release_area, .keep_all = TRUE)
  # Report losses
  cat("removed", x_rows - nrow(x), "duplicate release tags\n")
  cat("removed", y_rows - nrow(y), "duplicate recovery tags\n")

  # Filter x and y by days at liberty ------------------------------------------

  d <- dplyr::filter(y, recover_date - release_date < days_liberty)
  x <- dplyr::anti_join(x, d, by = colnames_x)[, colnames_x]
  y <- dplyr::inner_join(x, y, by = colnames_x)[, colnames_y]

  # Index areas from zero ------------------------------------------------------

  x <- dplyr::mutate(x, release_area = as.integer(release_area - 1L))
  y <- dplyr::mutate(y, release_area = as.integer(release_area - 1L))
  y <- dplyr::mutate(y, recover_area = as.integer(recover_area - 1L))

  # Convert groups to integer --------------------------------------------------

  x <- dplyr::mutate(x, group = mmmGroup(group_raw, groups))
  y <- dplyr::mutate(y, group = mmmGroup(group_raw, groups))

  # Convert date to step -------------------------------------------------------

  if (step == "year") {
    # Initial limit
    yr <- lubridate::year(limits[1])
    # x release step
    x <- dplyr::mutate(x, release_step = lubridate::year(release_date) - yr)
    x <- dplyr::mutate(x, release_step = as.integer(release_step))
    # y release step
    y <- dplyr::mutate(y, release_step = lubridate::year(release_date) - yr)
    y <- dplyr::mutate(y, release_step = as.integer(release_step))
    # y recover step
    y <- dplyr::mutate(y, recover_step = lubridate::year(recover_date) - yr)
    y <- dplyr::mutate(y, recover_step = as.integer(recover_step))
    # Steps at liberty
    steps_liberty <- round(exp(log(days_liberty) - log(365.25)))
  } else if (step == "month") {
    # Initial limits
    yr <- lubridate::year(limits[1])
    mo <- lubridate::month(limits[1])
    # x release step
    x <- dplyr::mutate(x, rel_yr = lubridate::year(release_date))
    x <- dplyr::mutate(x, rel_mo = lubridate::month(release_date))
    x <- dplyr::mutate(x, release_step = 12 * (rel_yr - yr) + (rel_mo - mo))
    x <- dplyr::mutate(x, release_step = as.integer(release_step))
    # y release step
    y <- dplyr::mutate(y, rel_yr = lubridate::year(release_date))
    y <- dplyr::mutate(y, rel_mo = lubridate::month(release_date))
    y <- dplyr::mutate(y, release_step = 12 * (rel_yr - yr) + (rel_mo - mo))
    y <- dplyr::mutate(y, release_step = as.integer(release_step))
    # y recover step
    y <- dplyr::mutate(y, rec_yr = lubridate::year(recover_date))
    y <- dplyr::mutate(y, rec_mo = lubridate::month(recover_date))
    y <- dplyr::mutate(y, recover_step = 12 * (rec_yr - yr) + (rec_mo - mo))
    y <- dplyr::mutate(y, recover_step = as.integer(recover_step))
    # Steps at liberty
    steps_liberty <- round(exp(log(days_liberty) - log(30)))
  } else if (step == "day") {
    # Initial limit
    d <- limits[1]
    # Steps
    x <- dplyr::mutate(x, release_step = as.integer(release_date - d))
    y <- dplyr::mutate(y, release_step = as.integer(release_date - d))
    y <- dplyr::mutate(y, recover_step = as.integer(recover_date - d))
    # Steps at liberty
    steps_liberty <- round(days_liberty)
  } else {
    stop("step must be one of 'year', 'month', or 'day'.")
  }
  # Print time step
  cat("using", step, "as the time step\n")

  # Select columns -------------------------------------------------------------

  x <- x %>% dplyr::select(release_step, release_area, group)
  y <- y %>% dplyr::select(
    release_step,
    release_area,
    group,
    recover_step,
    recover_area)

  # Create x -------------------------------------------------------------------

  x <- x %>%
    dplyr::group_by(release_step, release_area, group) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(release_step, release_area, group) %>%
    tidyr::drop_na() %>%
    as.matrix()
  # Confirm mode
  if (typeof(x) != "integer") {
    mode(x) <- "integer"
  }
  cat("x has", nrow(x), "rows and", sum(x[, "count"]), "tag releases\n")

  # Create y ------------------------------------------------------------------

  y <- y %>%
    dplyr::group_by(
      release_step,
      release_area,
      group,
      recover_step,
      recover_area) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(
      release_step,
      release_area,
      group,
      recover_step,
      recover_area) %>%
    tidyr::drop_na() %>%
    as.matrix()
  # Confirm mode
  if (typeof(y) != "integer") {
    mode(y) <- "integer"
  }
  cat("y has", nrow(y), "rows and", sum(y[, "count"]), "tag recoveries\n")

  # Stop the clock -------------------------------------------------------------

  tictoc::toc()

  # Return list of class mmmTags -----------------------------------------------

  return(structure(list(
    x = x,
    y = y,
    cols = cols,
    groups = groups,
    step = step,
    limits = limits,
    days_liberty = days_liberty,
    steps_liberty = steps_liberty),
    class = "mmmTags"))
}

#' Prepare Rate Data
#'
#' @description Prepare fishing mortality or tag reporting rate data prior
#' to fitting a Markov movement model via \code{mmmFit()}.
#'
#' @param x [data.frame()] Dates, areas, and rates data. See details.
#' @param cols [list()] Named list of character strings. See details.
#' @param reps [integer()] Number of tag time steps per rate time step (rows).
#' @param lims [integer()] First and last date (implemented as integer for now).
#' @param areas [list()] Named list of area values in index order. See details.
#' @param inst [logical()] Is the rate an instantaneous rate?
#'
#' @details The rate time step may differ from tag time step, as long as
#' the tag time step is a multiple of the rate time step.
#'
#' The \code{x} argument must be data frames that contains:
#'
#' \itemize{
#'   \item [integer()] column giving the date (integers for now).
#'   \item [atomic()] column giving the areas (see \code{areas} argument).
#'   \item [numeric()] column giving the rate for the given area and date.
#' }
#'
#' The \code{cols} argument is a [list()] of named elements. The elements
#'   must be [character()] string names of columns in \code{x}. There must
#'   be one element for each required element name:
#'   \itemize{
#'     \item \code{date}: Date as an integer (usually year, month, or day).
#'     \item \code{area}: Area as an atomic value.
#'     \item \code{rate}: Numeric rate corresponding to the date and area.
#'   }
#'
#' The \code{areas} argument is a [list()] of named [atomic()] elements.
#' The elements must correspond to elements in the area column of \code{x}.
#' The elements must be in the same order in which the areas are to be
#' indexed. The \code{areas} argument is optional if areas are already
#' indexed sequentially from one in \code{x}.
#'
#'
#' @return A matrix of class \code{mmmRates}.
#' @export
#'
#' @examples
#'
#' # One tag time step per row
#' d <- rep(2010:2016, 3)
#' a <- rep(c(1,2,3), each = 7)
#' r <- rep(seq(0.1, 0.7, length.out = 7), 3)
#' x <- data.frame(d = d, a = a, r = r)
#' cols <- list(date = "d", area = "a", rate = "r")
#' r1 <- mmmRates(x, cols)
#'
#' # Twelve tag time steps per row (e.g. monthly tags, yearly rates)
#' d <- rep(2010:2016, 3)
#' a <- rep(c(1,2,3), each = 7)
#' r <- rep(seq(0.1, 0.7, length.out = 7), 3)
#' x <- data.frame(d = d, a = a, r = r)
#' cols <- list(date = "d", area = "a", rate = "r")
#' reps <- 12
#' lims <- c(2010, 2016)
#' r2 <- mmmRates(x, cols, reps, lims, inst = TRUE)
#'
mmmRates <- function (x,
                      cols,
                      reps = 1,
                      lims = NULL,
                      areas = NULL,
                      inst = FALSE) {

  # Check arguments ------------------------------------------------------------

  # Check arguments
  checkmate::assert_data_frame(x)
  checkmate::assert_list(cols, types = "character", unique = TRUE)
  checkmate::assert_integerish(lims, len = 2, lower = 1, null.ok = TRUE)
  checkmate::assert_integerish(lims, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_list(areas, null.ok = TRUE)
  checkmate::assert_integerish(reps, lower = 1, len = 1, null.ok = FALSE)
  checkmate::assert_logical(inst, len = 1, null.ok = FALSE)
  # Check columns
  checkmate::assert_choice(cols$date, colnames(x))
  checkmate::assert_choice(cols$area, colnames(x))
  checkmate::assert_choice(cols$rate, colnames(x))

  # Rename columns -------------------------------------------------------------

  colnames(x)[which(colnames(x) == cols$date)] <- "date"
  colnames(x)[which(colnames(x) == cols$area)] <- "area"
  colnames(x)[which(colnames(x) == cols$rate)] <- "rate"

  # Define ---------------------------------------------------------------------

  # Define x
  colnames_x <- c("date", "area", "rate")
  x <- x[, colnames_x]
  # Date limits
  if (is.null(lims)) {
    lims <- c(min(x$date), max(x$date))
  }
  # Date indexes
  x <- dplyr::filter(x, date >= lims[1], date <= lims[2])
  x <- dplyr::mutate(x, date = date - lims[1] + 1L)

  # Check dates are sequential -------------------------------------------------

  # TODO: Check dates are sequential
  x <- dplyr::arrange(x, area, date)
  cat("caution: check that dates are sequential not implemented")

  # Convert areas to index -----------------------------------------------------

  # Indexed from one
  if (!is.null(areas)) {
    x <- dplyr::mutate(x, area = mmmGroup(area, areas) + 1L)
  }

  # Optional transformations ---------------------------------------------------

  # Deprecated
  # # Replicate rows
  # x <- x %>%
  #   dplyr::slice(rep(dplyr::row_number(), each = reps)) %>%
  #   dplyr::group_by(area) %>%
  #   dplyr::mutate(step = dplyr::row_number()) %>%
  #   dplyr::ungroup()

  # Convert instantaneous rate
  if (inst) {
    x$rate <- x$rate / reps
  }

  # Convert to step by area matrix ---------------------------------------------

  x <- x %>%
    dplyr::arrange(area) %>%
    tidyr::pivot_wider(names_from = area, values_from = rate) %>%
    dplyr::arrange(date) %>%
    dplyr::select(-date) %>%
    as.matrix()

  # Return ---------------------------------------------------------------------

  return(structure(x, class = "mmmRates"))
}

#' Create Monthly Weighting for Fishing Mortality Rates
#'
#' @description Create a weighting for the fishing mortality rate when the
#' fishing mortality rate is recorded at a larger time step than the tag
#' release and recovery time step. The weighting is based on tag recoveries
#' at the tag time step. Weights average to one for a given area.
#'
#' @param tags [mmmTags()] See [mmmTags()].
#' @param step [character()] Currently implemented for \code{"month"} only.
#'   Must equal the \code{step} element in \code{tags}.
#' @param nrows [integer()] Number of rows
#'
#' @return A matrix of monthly weighting for fishing mortality rates.
#'
#' @export
#'
mmmWeights <- function (tags, step = "month", nrows = 12L) {

  # Check arguments ------------------------------------------------------------

  checkmate::assert_class(tags, "mmmTags", null.ok = FALSE)
  checkmate::assert_choice(step, c("month"))
  checkmate::assert_choice(nrows, c(12))

  # Unpack arguments -----------------------------------------------------------

  y <- as.data.frame(tags$y)
  tag_step <- tags$step

  # Check time steps -----------------------------------------------------------

  if (step != tag_step) {
    stop("weight step must equal tag step")
  }

  # Create the data frame ------------------------------------------------------

  area <- rep(c(seq_len(max(y$recover_area) + 1L) - 1L), each = nrows)
  step <- rep(c(seq_len(nrows) - 1L), max(y$recover_area) + 1L)
  d <- data.frame(recover_area = area, weight_step = step)

  # Create the weights matrix --------------------------------------------------

  w <- y %>%
    dplyr::mutate(weight_step = recover_step %% 12L) %>%
    dplyr::select(recover_area, weight_step, count) %>%
    dplyr::arrange(recover_area, weight_step) %>%
    dplyr::group_by(recover_area, weight_step) %>%
    dplyr::mutate(weight_count = sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(recover_area, weight_step, weight_count) %>%
    dplyr::group_by(recover_area) %>%
    dplyr::mutate(weight = weight_count / sum(weight_count)) %>%
    dplyr::ungroup() %>%
    dplyr::full_join(d, by = c("recover_area", "weight_step")) %>%
    tidyr::replace_na(list(weight = 0)) %>%
    dplyr::select(-weight_count) %>%
    tidyr::pivot_wider(names_from = recover_area, values_from = weight) %>%
    dplyr::select(-weight_step) %>%
    as.matrix()

  # Convert from (sum to one) to (average to one) ------------------------------

  w <- w * 12

  # Return the weights matrix --------------------------------------------------

  return(structure(w, class = "mmmWeights"))
}


#' Create Movement Index Matrix
#'
#' @param n [integer()] Number of areas.
#' @param pattern [integer()] One of \code{0}: full movement, or \code{1}:
#'   direct movement between numerically sequential neighbours only.
#' @param allow [integer()] [matrix()] Each row indicates directional
#'   movement between a pair of areas indexed from one.
#' @param disallow [integer()] [matrix()] As for \code{allow}, but
#'   specified movement is disallowed.
#'
#' @return A square binary matrix with zero diagonal of class \code{mmmIndex}
#' @export
#'
#' @examples
#'
#' # Neighbours (default)
#' mmmIndex(6)
#'
#' # Full movement
#' mmmIndex(6, 0)
#'
#' # Neighbours plus 1-6, 2-5, and 5-2, but not 6-1
#' allow <- matrix(c(1,6,2,5,5,2), ncol = 2, byrow = TRUE)
#' mmmIndex(6, 1, allow)
#'
#' # Full minus 1-6
#' disallow <- matrix(c(1,6), ncol = 2, byrow = TRUE)
#' mmmIndex(6, 0, disallow = disallow)
#'
mmmIndex <- function (n, pattern = NULL, allow = NULL, disallow = NULL) {

  # Check arguments
  checkmate::assert_integerish(n, lower = 2, len = 1, any.missing = FALSE)
  checkmate::assert_integerish(pattern, len = 1, null.ok = TRUE)
  checkmate::assert_integerish(pattern, lower = 0, upper = 1, null.ok = TRUE)
  checkmate::assert_integerish(pattern, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_matrix(allow, mode = "integerish", null.ok = TRUE)
  checkmate::assert_matrix(allow, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_matrix(allow, ncols = 2, null.ok = TRUE)
  # Initialize matrix
  z <- matrix(0L, nrow = n, ncol = n)
  # Default
  if (is.null(pattern) & is.null(allow)) {
    # Immediate neighbours only
    for (i in seq_len(n - 1L)) {
      z[i, i + 1L] <- 1L
      z[i + 1L, i] <- 1L
    }
  }
  # Add pattern
  if (!is.null(pattern)) {
    if (pattern) {
      for (i in seq_len(n - 1L)) {
        z[i, i + 1L] <- 1L
        z[i + 1L, i] <- 1L
      }
    } else {
      z <- matrix(1L, nrow = n, ncol = n)
      diag(z) <- 0L
    }
  }
  # Allow indexes
  if (!is.null(allow)) {
    z[allow] <- 1L
    diag(z) <- 0L
  }
  # Disallow indexes
  if (!is.null(disallow)) {
    z[disallow] <- 0L
    diag(z) <- 0L
  }
  # Return
  return(structure(z, class = "mmmIndex"))
}

#' Map Values to Groups
#'
#' @param x [atomic()] Vector of values.
#' @param groups [list()] Named vectors of group elements.
#'
#' @return [integer()] Vector.
#'
#' @export
#'
#' @examples
#'
#' # Numeric
#' x <- c(1, 4, 2, 3, NA, 2, 7)
#' y <- list(a = 1:2, b = 3:4)
#' g1 <- mmmGroup(x, y)
#' # Character
#' x <- c("M", "F", "F", NA, "M", "N", "F")
#' y <- list(m = "M", f = "F")
#' g2 <- mmmGroup(x, y)
#' # Factor
#' x <- factor(c("L", "S", "M", NA, "M", "H"))
#' y <- list(s = "S", m = "M", l = "L")
#' g3 <- mmmGroup(x, y)
#'
mmmGroup <- function (x, groups) {

  #--------------- Check arguments --------------------------------------------#

  # TODO: Check that all names of x are unique
  # TODO: Check that all elements of x are vectors of the same type
  # TODO: Check that all elements of elements of x are unique
  checkmate::assert_vector(x)
  checkmate::assert_list(groups, unique = TRUE)

  #--------------- Create index -----------------------------------------------#

  s <- c(seq_along(groups) - 1L)
  x <- data.frame(value = x)
  d <- data.frame(value = unlist(groups), index = rep(s, lengths(groups)))

  #--------------- Return vector ----------------------------------------------#

  return(as.integer(dplyr::left_join(x, d, by = "value")$index))
}
