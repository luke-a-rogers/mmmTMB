#' Prepare Tag Data
#'
#' @param x A data frame of named columns. See details.
#' @param y A data frame of named columns. See details.
#' @param cols A list of named character elements. See details.
#' @param time_step A character string. One of \code{year}, \code{month}, or
#'   \code{day}.
#' @param date_lims A character vector with two elements. "%Y-%M-%D".
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
#'   Release_Group = c(1,1,1),
#'   ID = c("TAG001", "TAG002", "TAG003"))
#' y <- data.frame(
#'   Release_Date = c("2010-06-01", "2020-06-01"),
#'   Release_Area = c(1,1),
#'   Recover_Date = c("2013-08-01", "2020-08-01"),
#'   Recover_Area = c(1,1),
#'   Release_Group = c(1,1),
#'   ID = c("TAG001", "TAG003"))
#' cols = list(
#'   release_date = "Release_Date",
#'   release_area = "Release_Area",
#'   recover_date = "Recover_Date",
#'   recover_area = "Recover_Area",
#'   group = "Release_Group",
#'   id = "ID")
#' time_step <- "year"
#' date_lims <- c("2010-01-01", "2020-12-31")
#'
#'
#'
#'
mmmTags <- function (x,
                     y,
                     cols = NULL,
                     time_step = NULL,
                     date_lims = NULL,
                     days_liberty = c(0L, Inf)) {

  #--------------- Check arguments --------------------------------------------#

  checkmate::assert_data_frame(x)
  checkmate::assert_data_frame(y)
  checkmate::assert_list(cols, unique = TRUE)
  checkmate::assert_character(time_step, len = 1)
  checkmate::assert_character(date_lims, len = 2, null.ok = TRUE)
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
    release_date <= date_lims[2]
  )
  # Recover
  y <- y %>% dplyr::filter(
    release_date >= date_lims[1],
    release_date <= date_lims[2],
    recover_date >= date_lims[1],
    recover_date <= date_lims[2],
    release_date < recover_date
  )

  # TODO: Check all
  # - all x id unique
  # - y id / release_date / release_area in x (use inner_join)
  # - Decide: exclude releases from x when recoveries excluded by days_liberty?
  # (I think yes - other exclusions (e.g. data errors roll into reporting rate))


  #--------------- Convert date to time step ----------------------------------#

  # TODO: day, month, or year


  #--------------- Convert group to integer -----------------------------------#

  # TODO:
  # - Categorical character like sex
  # - Factor with levels
  # - Continuous with ranges


  #--------------- Create mT --------------------------------------------------#

  # TODO: Handle summation
  # TODO: Decrement indices for c++


  #--------------- Create mR --------------------------------------------------#

  # TODO: Handle summation
  # TODO: Decrement indices for c++

  #--------------- Return a list ----------------------------------------------#

  return(list(
    structure(mT = mT, class = "mT"),
    structure(mR = mR, class = "mR"))
  )
}


