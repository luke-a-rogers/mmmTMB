#' Prepare Tag Data
#'
#' @param x A data frame of named columns. See details.
#' @param y A data frame of named columns. See details.
#' @param cols A list of named character elements. See details.
#' @param groups A list of named vectors. See details.
#' @param time_step A character string. One of \code{year}, \code{month}, or
#'   \code{day}.
#' @param date_lims A character vector with two elements. "%Y-%M-%D".
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
#'   Release_Date = c("2010-06-01", "2015-04-01", "2020-06-01"),
#'   Release_Area = c(1,1,1),
#'   Release_Group = c("F","F","M"),
#'   ID = c("TAG001", "TAG002", "TAG003"))
#' y <- data.frame(
#'   Release_Date = c("2010-06-01", "2020-06-01"),
#'   Release_Area = c(1,1),
#'   Recover_Date = c("2013-08-01", "2020-08-01"),
#'   Recover_Area = c(1,1),
#'   Release_Group = c("F", "M"),
#'   ID = c("TAG001", "TAG003"))
#' cols <- list(
#'   release_date = "Release_Date",
#'   release_area = "Release_Area",
#'   recover_date = "Recover_Date",
#'   recover_area = "Recover_Area",
#'   group = "Release_Group",
#'   id = "ID")
#' groups <- list(m = "M", f = "F")
#' time_step <- "year"
#' date_lims <- c("2010-01-01", "2020-12-31")
#' days_liberty <- c(0L, Inf)
#'
#' tags <- mmmTags(x, y, cols, groups, time_step, date_lims, days_liberty)
#'
#'
mmmTags <- function (x,
                     y,
                     cols = NULL,
                     groups = NULL,
                     time_step = NULL,
                     date_lims = NULL,
                     days_liberty = c(0L, Inf)) {

  #--------------- Check arguments --------------------------------------------#

  checkmate::assert_data_frame(x)
  checkmate::assert_data_frame(y)
  checkmate::assert_list(cols, unique = TRUE)
  checkmate::assert_list(groups, unique = TRUE, null.ok = TRUE)
  checkmate::assert_character(time_step, len = 1, null.ok = TRUE)
  checkmate::assert_numeric(days_liberty, lower = 0, len = 2)

  #--------------- Check required columns -------------------------------------#

  checkmate::assert_true(is.element(cols$release_date, colnames(x)))
  checkmate::assert_true(is.element(cols$release_area, colnames(x)))
  checkmate::assert_true(is.element(cols$release_date, colnames(y)))
  checkmate::assert_true(is.element(cols$release_area, colnames(y)))
  checkmate::assert_true(is.element(cols$recover_date, colnames(y)))
  checkmate::assert_true(is.element(cols$recover_area, colnames(y)))
  checkmate::assert_true(is.element(cols$id, colnames(x)))
  checkmate::assert_true(is.element(cols$id, colnames(y)))

  #--------------- Check optional columns -------------------------------------#

  if (is.null(cols$group)) {
    cols$group <- "group"
    x$group <- 1L
    y$group <- 1L
  }
  checkmate::assert_true(is.element(cols$group, colnames(x)))
  checkmate::assert_true(is.element(cols$group, colnames(y)))

  #--------------- Rename columns ---------------------------------------------#

  # Rename columns x
  colnames(x)[which(colnames(x) == cols$release_date)] <- "release_date"
  colnames(x)[which(colnames(x) == cols$release_area)] <- "release_area"
  colnames(x)[which(colnames(x) == cols$group)] <- "group"
  colnames(x)[which(colnames(x) == cols$id)] <- "id"
  # Rename columns y
  colnames(y)[which(colnames(y) == cols$release_date)] <- "release_date"
  colnames(y)[which(colnames(y) == cols$release_area)] <- "release_area"
  colnames(y)[which(colnames(y) == cols$recover_date)] <- "recover_date"
  colnames(y)[which(colnames(y) == cols$recover_area)] <- "recover_area"
  colnames(y)[which(colnames(y) == cols$group)] <- "group"
  colnames(y)[which(colnames(y) == cols$id)] <- "id"

  #--------------- Define x and y ---------------------------------------------#

  release_cols <- c("release_date", "release_area")
  recover_cols <- c("recover_date", "recover_area")
  specify_cols <- c("group", "id")
  colnames_x <- c(release_cols, specify_cols)
  colnames_y <- c(release_cols, recover_cols, specify_cols)
  # Define
  x <- x[, colnames_x]
  y <- y[, colnames_y]

  #--------------- Check column classes ---------------------------------------#

  # TODO: Check column classes
  # TODO: Check that groups classes and group column classes match


  #--------------- Drop NA values ---------------------------------------------#

  x <- tidyr::drop_na(x)
  y <- tidyr::drop_na(y)

  #--------------- Convert columns to dates -----------------------------------#

  x <- dplyr::mutate(x, release_date = lubridate::as_date(release_date))
  y <- dplyr::mutate(y,
    release_date = lubridate::as_date(release_date),
    recover_date = lubridate::as_date(recover_date)
  )

  #--------------- Convert date limits to dates -------------------------------#

  if (is.null(date_lims)) {
    date_lims <- c(min(x$release_date), max(c(x$release_date, y$recover_date)))
  } else {
    date_lims <- lubridate::as_date(date_lims)
  }

  #--------------- Filter events by dates -------------------------------------#

  # Release
  x <- x %>% dplyr::filter(
    release_date >= date_lims[1],
    release_date <= date_lims[2])
  # Recover
  y <- y %>% dplyr::filter(
    release_date >= date_lims[1],
    release_date <= date_lims[2],
    recover_date >= date_lims[1],
    recover_date <= date_lims[2],
    release_date <= recover_date)

  #--------------- Remove duplicates from x and y -----------------------------#

  # The tag reporting rate accounts for data losses due to errors
  x_rows <- nrow(x)
  y_rows <- nrow(y)
  x <- dplyr::distinct(x, id, release_date, release_area, .keep_all = TRUE)
  y <- dplyr::distinct(y, id, release_date, release_area, .keep_all = TRUE)
  cat(paste0("duplicate release tags removed: ", x_rows - nrow(x)), "\n")
  cat(paste0("duplicate recover tags removed: ", y_rows - nrow(y)), "\n")

  #--------------- Filter y by x ----------------------------------------------#

  y <- dplyr::inner_join(x, y, by = colnames_x)[, colnames_y]

  #--------------- Filter y by days at liberty --------------------------------#

  # Time steps at liberty must also be tracked in the movement model
  y <- y %>% dplyr::filter(
    recover_date - release_date >= days_liberty[1],
    recover_date - release_date <= days_liberty[2])

  #--------------- Convert date to time step ----------------------------------#

  if (is.null(time_step) || time_step == "year") {
    x <- x %>% dplyr::mutate(
      release_step = as.integer(
        lubridate::year(release_date) - lubridate::year(date_lims[1])))
    y <- y %>% dplyr::mutate(
      release_step = as.integer(
        lubridate::year(release_date) - lubridate::year(date_lims[1])),
      recover_step = as.integer(
        lubridate::year(recover_date) - lubridate::year(date_lims[1])))
    # Convert liberty
    steps_liberty <- round(exp(log(days_liberty) - log(365.25)))
    cat(paste("days at liberty:", days_liberty[1], days_liberty[2]), "\n")
    cat(paste("steps at liberty:", steps_liberty[1], steps_liberty[2]), "\n")
    cat("use caution if days and steps (years) at liberty disagree.\n")
  } else if (time_step == "month") {
    x <- x %>% dplyr::mutate(
      release_step = as.integer(
        12 * (lubridate::year(release_date) - lubridate::year(date_lims[1])) +
          lubridate::month(release_date) - lubridate::month(date_lims[1])))
    y <- y %>% dplyr::mutate(
      release_step = as.integer(
        12 * (lubridate::year(release_date) - lubridate::year(date_lims[1])) +
          lubridate::month(release_date) - lubridate::month(date_lims[1])),
      recover_step = as.integer(
        12 * (lubridate::year(recover_date) - lubridate::year(date_lims[1])) +
          lubridate::month(recover_date) - lubridate::month(date_lims[1])))
    # Convert liberty
    steps_liberty <- round(exp(log(days_liberty) - log(30)))
    cat(paste("days at liberty:", days_liberty[1], days_liberty[2]), "\n")
    cat(paste("steps at liberty:", steps_liberty[1], steps_liberty[2]), "\n")
    cat("use caution if days and steps (years) at liberty disagree.\n")
  } else if (time_step == "day") {
    x <- x %>% dplyr::mutate(
      release_step = as.integer(release_date - date_lims[1]))
    y <- y %>% dplyr::mutate(
      release_step = as.integer(release_date - date_lims[1]),
      recover_step = as.integer(recover_date - date_lims[1]))
    steps_liberty <- days_liberty
    cat(paste("days at liberty:", days_liberty[1], days_liberty[2]), "\n")
  } else {
    stop("time_step must be one of 'year', 'month', or 'day'.")
  }

  #--------------- Index areas from zero --------------------------------------#

  x <- dplyr::mutate(x, release_area = as.integer(release_area - 1L))
  y <- dplyr::mutate(y, release_area = as.integer(release_area - 1L),
                     recover_area = as.integer(recover_area - 1L))

  #--------------- Convert group to integer -----------------------------------#

  x <- dplyr::mutate(x, group_index = mmmGroup(group, groups))
  y <- dplyr::mutate(y, group_index = mmmGroup(group, groups))

  #--------------- Drop NA values ---------------------------------------------#

  x <- tidyr::drop_na(x)
  y <- tidyr::drop_na(y)

  #--------------- Select columns ---------------------------------------------#

  x <- dplyr::select(x, release_step, release_area, group_index)
  y <- dplyr::select(y, release_step, release_area, recover_step, recover_area,
                     group_index)


  #--------------- Create mT --------------------------------------------------#

  mT <- dplyr::group_by(x, release_step, release_area, group_index) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    as.matrix()
  if (typeof(mT) != "integer") {mode(mT) <- "integer"}

  #--------------- Create mR --------------------------------------------------#

  mR <- dplyr::group_by(y, release_step, release_area, recover_step,
                        recover_area, group_index) %>%
    dplyr::mutate(count = dplyr::n()) %>%
    dplyr::ungroup() %>%
    dplyr::distinct(.keep_all = TRUE) %>%
    as.matrix()
  if (typeof(mR) != "integer") {mode(mR) <- "integer"}

  #--------------- Return a list ----------------------------------------------#

  # Include appropriate units at liberty
  return(structure(list(
    mT = mT,
    mR = mR,
    groups = groups,
    time_step = time_step,
    date_lims = date_lims,
    days_liberty = days_liberty,
    steps_liberty = steps_liberty),
    class = "mmmTags"))
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
