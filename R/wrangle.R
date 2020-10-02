#' Array T From Data Frame
#'
#' @description Make array T from a data frame.
#'
#' @param x A data frame of tag releases. See details.
#' @param cols A list matching required names to column names in x. See details.
#' @param replace_na A number to replace NA values before they are removed.
#'
#' @details The data frame \code{x} must include columns for tag release area
#'   and tag release time. Currently, both should be integer valued and start
#'   indexing at one. If each row is corresponds to one tag release, no further
#'   columns are needed. If each row summarizes a number of tag releases in one
#'   area at one time, then \code{x} must include a column of release counts.
#'   Optionally, an integer column indicating individual class (e.g. sex, length
#'   class) can be included.
#'
#'   The list \code{cols} must include
#'   \itemize{
#'     \item \code{}
#'   }
#'
#'
#' @return An array of class \code{arrayT}
#' @export
#'
#' @examples
#'
arrayT <- function(x,
                   cols,
                   replace_na = NA) {

}

# arrayR <- function(x,
#                    cols,
#                    replace_na = NA) {
#
# }


