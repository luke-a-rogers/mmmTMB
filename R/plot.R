#' Plot Movement Probabilities from an 'mmmTMB' Object
#'
#' @param x An object of class \code{mmmTMB}
#' @param class [integer()] Release class
#' @param xlab [character()] The x label of the plot
#' @param ylab [character()] The y label of the plot
#' @param area_names [character()] Vector of area names
#' @param font_size_probs [numeric()] Font size for probabilities
#' @param font_nudge_probs [numeric()] Upward displacement for probabilities
#' @param font_size_stderr [numeric()] Font size for standard errors
#' @param font_nudge_stderr [numeric()] Downward displacement for standard errors
#' @param legend_name [character()] Legend title
#' @param legend_position [character()] Legend position
#'
#' @return An object with class \code{gg} and \code{ggplot}
#'
#' @export
#'
#' @examples
#'
plot.mmmTMB <- function (x = NULL,
                         class = 1,
                         xlab = "Area To",
                         ylab = "Area From",
                         area_names = NULL,
                         font_size_probs = 3,
                         font_nudge_probs = 0.15,
                         font_size_stderr = 2,
                         font_nudge_stderr = 0.15,
                         legend_name = "Estimate",
                         legend_position = "right") {

  #---------------- Check arguments -------------------------------------------#


  #---------------- Extract movement results data -----------------------------#

  if (is.element("mmmTMB", class(x))) {
    x <- x$results$movement_probability_results
    x_inds <- which(x$Class == class)
    x <- x[x_inds,]
  } else if (is.data.frame(x)) {
    # TODO: Check appropriate columns are present
    x_inds <- which(x$Class == class)
    x <- x[x_inds,]
  } else if (!is.null(x)) {
    warning("x not as expected")
  }

  #---------------- Convert columns to factor ---------------------------------#

  x$Area_From <- factor(
    x$Area_From,
    levels = sort(unique(x$Area_From), decreasing = TRUE))
  x$Area_To <- factor(
    x$Area_To,
    levels = sort(unique(x$Area_To)))

  #---------------- Construct geom object -------------------------------------#

  ggplot2::ggplot(
    data = x,
    mapping = ggplot2::aes(x = Area_To, y = Area_From, fill = Estimate),
  ) +
    ggplot2::geom_tile(
      color = "white",
      width = 0.975,
      height = 0.975) +
    # Use viridis?
    ggplot2::scale_fill_viridis_c(
      direction = 1,
      option = "plasma",
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    # Add text
    ggplot2::theme_minimal() +
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        label = round(Estimate, 2),
        col = as.factor(ifelse(Estimate >= 0.5, 0, 1))),
      fontface = "plain",
      nudge_y = font_nudge_probs,
      size = font_size_probs) +
    # Add SE
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        label = paste0("(", round(SE, 3), ")"),
        col = as.factor(ifelse(Estimate >= 0.5, 0, 1))),
      fontface = "plain",
      nudge_y = -font_nudge_stderr,
      size = font_size_stderr) +
    ggplot2::scale_color_grey(start = 0, end = 1) +
    ggplot2::guides(col = FALSE) +
    ggplot2::xlab(xlab) +
    ggplot2::ylab(ylab) +
    ggplot2::labs(fill = "Estimate") +
    ggplot2::scale_x_discrete(position = "bottom", labels = area_names) +
    ggplot2::scale_y_discrete(position = "left", labels = rev(area_names)) +
    # Add theme
    ggplot2::theme(
      legend.position = legend_position,
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 12, face = "plain"),
      axis.title.y = ggplot2::element_text(size = 12, face = "plain"),
      axis.text.y = ggplot2::element_text(
        size = 10,
        face = "plain",
        color = "black",
        margin = ggplot2::margin(t = 0, r = -2, b = 0, l = 5)),
      axis.text.x = ggplot2::element_text(
        size = 10,
        face = "plain",
        color = "black",
        margin = ggplot2::margin(t = -2)))
}
