#' Prepare Tag Data
#'
#' @param x A data frame of named columns. See details.
#' @param y A data frame of named columns. See details.
#' @param cols A list of named character elements. See details.
#' @param groups A list of named vectors. See details.
#' @param time_step A character string. One of \code{year}, \code{month}, or
#'   \code{day}.
#' @param limits A character vector with two elements. "%Y-%M-%D".
#' @param days_liberty Integer scalar. Minimum number of days at liberty.
#'
#' @details TBD
#'
#' @importFrom magrittr `%>%`
#'
#' @return
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
#' time_step <- "year"
#' limits <- c("2010-01-01", "2020-12-31")
#' days_liberty <- 0L
#'
#' tags <- mmmTags(x, y, cols, groups, time_step, limits, days_liberty)
#'
#'
mmmTags <- function (x,
                     y,
                     cols = NULL,
                     groups = NULL,
                     step = NULL,
                     limits = NULL,
                     days_liberty = NULL) {

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
    cat("using no groups: groups or cols$group null")
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
  checkmate::assert_integerish(x$group_raw, lower = 1, any.missing = FALSE)
  checkmate::assert_character(x$id, len = 1, any.missing = FALSE)
  # Data frame y
  checkmate::assert_date(y$release_date, any.missing = FALSE)
  checkmate::assert_date(y$recover_date, any.missing = FALSE)
  checkmate::assert_integerish(y$release_area, lower = 1, any.missing = FALSE)
  checkmate::assert_integerish(y$recover_area, lower = 1, any.missing = FALSE)
  checkmate::assert_integerish(y$group_raw, lower = 1, any.missing = FALSE)
  checkmate::assert_character(y$id, len = 1, any.missing = FALSE)

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
  cat("removed", y_rows - nrow(y), "duplicate recover tags\n")

  # Filter x and y by days at liberty ------------------------------------------

  d <- dplyr::filter(y, recover_date - release_date < days_liberty)
  x <- dplyr::anti_join(x, d, by = colnames_x)[, colnames_x]
  y <- dplyr::inner_join(x, y, by = colnames_x)[, colnames_y]

  # Index areas from zero ------------------------------------------------------

  x <- dplyr::mutate(x, release_area = as.integer(release_area - 1L))
  y <- dplyr::mutate(y, release_area = as.integer(release_area - 1L))
  y <- dplyr::mutate(y, recover_area = as.integer(recover_area - 1L))

  # Convert groups to integer --------------------------------------------------

  x <- dplyr::mutate(x, group_index = mmmGroup(group, groups))
  y <- dplyr::mutate(y, group_index = mmmGroup(group, groups))


  # TODO: Continue from here
  # TODO: Check that group is identical between release and recover



  #--------------- Convert date to time step ----------------------------------#

  if (step == "year") {
    x <- x %>% dplyr::mutate(
      release_step = as.integer(
        lubridate::year(release_date) - lubridate::year(limits[1])))
    y <- y %>% dplyr::mutate(
      release_step = as.integer(
        lubridate::year(release_date) - lubridate::year(limits[1])),
      recover_step = as.integer(
        lubridate::year(recover_date) - lubridate::year(limits[1])))
    # Convert liberty
    steps_liberty <- round(exp(log(days_liberty) - log(365.25)))
    cat(paste0("days at liberty: ", days_liberty,"; "))
    cat(paste0("steps (years) at liberty: ", steps_liberty[1]), "\n")
  } else if (step == "month") {
    x <- x %>% dplyr::mutate(
      release_step = as.integer(
        12 * (lubridate::year(release_date) - lubridate::year(limits[1])) +
          lubridate::month(release_date) - lubridate::month(limits[1])))
    y <- y %>% dplyr::mutate(
      release_step = as.integer(
        12 * (lubridate::year(release_date) - lubridate::year(limits[1])) +
          lubridate::month(release_date) - lubridate::month(limits[1])),
      recover_step = as.integer(
        12 * (lubridate::year(recover_date) - lubridate::year(limits[1])) +
          lubridate::month(recover_date) - lubridate::month(limits[1])))
    # Convert liberty
    steps_liberty <- round(exp(log(days_liberty) - log(30)))
    cat(paste0("days at liberty: ", days_liberty[1], "; "))
    cat(paste0("steps (months) at liberty: ", steps_liberty[1]), "\n")
  } else if (step == "day") {
    x <- x %>% dplyr::mutate(
      release_step = as.integer(release_date - limits[1]))
    y <- y %>% dplyr::mutate(
      release_step = as.integer(release_date - limits[1]),
      recover_step = as.integer(recover_date - limits[1]))
    steps_liberty <- days_liberty
    cat(paste0("days at liberty: ", days_liberty), "; ")
    cat(paste0("steps (days) at liberty: ", steps_liberty[1]), "\n")
  } else {
    stop("step must be one of 'year', 'month', or 'day'.")
  }


  #--------------- Select columns ---------------------------------------------#

  x <- dplyr::select(x, release_step, release_area, group_index)
  y <- dplyr::select(y, release_step, release_area, group_index,
                     recover_step, recover_area)

  #--------------- Create mT --------------------------------------------------#

  mT <- dplyr::group_by(x, release_step, release_area, group_index) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(release_step, release_area, group_index) %>%
    as.matrix()
  if (typeof(mT) != "integer") { mode(mT) <- "integer" }

  #--------------- Create mR --------------------------------------------------#

  mR <- dplyr::group_by(y, release_step, release_area, group_index) %>%
    dplyr::group_by(recover_step, recover_area, .add = TRUE) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    dplyr::arrange(release_step, release_area, group_index,
                   recover_step, recover_area) %>%
    as.matrix()
  if (typeof(mR) != "integer") { mode(mR) <- "integer" }

  #--------------- Return a list ----------------------------------------------#

  # Include appropriate units at liberty
  return(structure(list(
    mT = mT,
    mR = mR,
    groups = groups,
    step = step,
    limits = limits,
    days_liberty = days_liberty,
    steps_liberty = steps_liberty),
    class = "mmmTags"))
}

#' Prepare Rate Data
#'
#' @param x A data frame of named columns. See details.
#' @param cols A list of named character elements. See details.
#' @param lims An integer vector of length 2.
#' @param areas A character vector of area names.
#' @param reps An integer. Number of time steps per row.
#' @param inst Logical. Is the rate an instantaneous rate?
#'
#' @return
#' @export
#'
#' @examples
#' date <- rep(2010:2016, 3)
#' area <- rep(c(1,2,3), each = 7)
#' rate <- rep(seq(0.1, 0.5, length.out = 7), 3)
#' x <- data.frame(date = date, area = area, rate = rate)
#' r1 <- mmmRates(x)
#'
mmmRates <- function (x,
                      cols = NULL,
                      lims = NULL,
                      areas = NULL,
                      reps = 1,
                      inst = FALSE) {

  #--------------- Check arguments --------------------------------------------#

  # TODO: Check all arguments
  checkmate::assert_data_frame(x)
  checkmate::assert_list(cols, unique = TRUE, null.ok = TRUE)
  checkmate::assert_integerish(lims, len = 2, lower = 0, null.ok = TRUE)
  checkmate::assert_integerish(lims, any.missing = FALSE, null.ok = TRUE)
  checkmate::assert_list(areas, null.ok = TRUE)
  checkmate::assert_integerish(reps, lower = 1, len = 1, null.ok = FALSE)
  checkmate::assert_logical(inst, len = 1, null.ok = FALSE)

  #--------------- Set default cols -------------------------------------------#

  if (is.null(cols)) {
    cols <- list(date = "date", area = "area", rate = "rate")
  }

  #--------------- Check required columns -------------------------------------#

  checkmate::assert_true(is.element(cols$date, colnames(x)))
  checkmate::assert_true(is.element(cols$area, colnames(x)))
  checkmate::assert_true(is.element(cols$rate, colnames(x)))

  #--------------- Rename columns ---------------------------------------------#

  colnames(x)[which(colnames(x) == cols$date)] <- "date"
  colnames(x)[which(colnames(x) == cols$area)] <- "area"
  colnames(x)[which(colnames(x) == cols$rate)] <- "rate"

  #--------------- Define x ---------------------------------------------------#

  colnames_x <- c("date", "area", "rate")
  x <- x[, colnames_x]

  #--------------- Define date limits -----------------------------------------#

  if (is.null(lims)) { lims <- c(min(x$date), max(x$date)) }

  #--------------- Convert dates to index -------------------------------------#

  x <- dplyr::filter(x, date >= lims[1], date <= lims[2]) %>%
    dplyr::mutate(date = date - lims[1] + 1L)

  #--------------- Check dates are sequential ---------------------------------#

  # TODO: Check dates are sequential

  #--------------- Convert areas to integers ----------------------------------#

  if (is.character(x$area)) {
    x <- dplyr::mutate(x, area = mmmGroup(area, areas) + 1L)
  }

  #--------------- Optionally replicate rows ----------------------------------#

  x <- dplyr::slice(x, rep(dplyr::row_number(), each = reps)) %>%
    dplyr::group_by(area) %>%
    dplyr::mutate(step = dplyr::row_number()) %>%
    dplyr::ungroup()

  #--------------- Optionally convert rate ------------------------------------#

  if (inst) { x$rate <- x$rate / reps }

  #--------------- Convert to step by area matrix -----------------------------#

  x <- dplyr::arrange(x, area) %>%
    tidyr::pivot_wider(names_from = area, values_from = rate) %>%
    dplyr::arrange(step) %>%
    dplyr::select(-date, -step) %>%
    as.matrix()
  return(x)
}

#' Create Monthly Weighting for Annual Fishing Mortality
#'
#' @param tags An \code{mmmTags} object. See \code{mmmTags()}.
#' @param rate_step String. Currenlty implemented for \code{"month"} only.
#'
#' @return A matrix of monthly weights for fishing mortality rates.
#'
#' @export
#'
#' @examples
#'
mmmWeights <- function (tags,
                        rate_step = "year") {

  #--------------- Check arguments --------------------------------------------#

  checkmate::assert_class(tags, "mmmTags", null.ok = FALSE)

  #---------------- Unpack arguments ------------------------------------------#

  mR <- as.data.frame(tags$mR)
  time_step <- tags$time_step

  #---------------- Create area-month data frame ------------------------------#

  recover_area <- rep(c(seq_len(max(mR$recover_area) + 1L) - 1L), each = 12)
  month_index <- rep(c(seq_len(12) - 1L), max(mR$recover_area) + 1L)
  d <- data.frame(recover_area = recover_area, month_index = month_index)

  #---------------- Create weights matrix -------------------------------------#

  if (time_step == "month" & rate_step == "year") {
    mW <- dplyr::mutate(mR, month_index = recover_step %% 12L) %>%
      dplyr::select(recover_area, month_index) %>%
      dplyr::group_by(recover_area, month_index) %>%
      dplyr::mutate(count = dplyr::n()) %>%
      dplyr::ungroup() %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      dplyr::group_by(recover_area) %>%
      dplyr::mutate(weight = count / sum(count)) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(recover_area, month_index) %>%
      dplyr::full_join(d, by = c("recover_area", "month_index")) %>%
      dplyr::arrange(recover_area, month_index) %>%
      tidyr::replace_na(list(weight = 0)) %>%
      dplyr::select(-count) %>%
      tidyr::pivot_wider(names_from = recover_area, values_from = weight) %>%
      dplyr::select(-month_index) %>%
      as.matrix()
  } else {
    mW <- NULL
    cat("warning: currently implemented for tags$time_step == 'month' only \n")
  }

  #--------------- Return weights matrix --------------------------------------#

  return(mW)
}


#' Create Movement Index Matrix
#'
#' @param n [integer()] Number of areas.
#' @param pattern [integer()] One of \code{0}: full movement, or \code{1}:
#'   direct movement between numerically sequential neighbours only
#' @param allow [integer()] [matrix()] Each row indicates directional
#'   movement between a pair of areas indexed from \code{1}.
#' @param disallow [integer()] [matrix()] As for \code{allow}, but
#'   specified movement is disallowed.
#'
#' @return A square binary matrix
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
  m <- matrix(0L, nrow = n, ncol = n)
  # Default
  if (is.null(pattern) & is.null(allow)) {
    # Immediate neighbours only
    for (i in seq_len(n - 1L)) {
      m[i, i + 1L] <- 1L
      m[i + 1L, i] <- 1L
    }
  }
  # Add pattern
  if (!is.null(pattern)) {
    if (pattern) {
      for (i in seq_len(n - 1L)) {
        m[i, i + 1L] <- 1L
        m[i + 1L, i] <- 1L
      }
    } else {
      m <- matrix(1L, nrow = n, ncol = n)
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

#' Map Values to Groups
#'
#' @param x A vector
#' @param groups A list of named group vectors
#'
#' @return An integer vector
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

  z <- c(seq_along(groups) - 1L)
  x <- data.frame(value = x)
  d <- data.frame(value = unlist(groups), index = rep(z, lengths(groups)))

  #--------------- Return vector ----------------------------------------------#

  return(as.integer(dplyr::left_join(x, d, by = "value")$index))
}
