#' Plot Movement Rate Heatmap Matrix
#'
#' @param x An object of class \code{mmmFit}
#' @param group [integer()] Release group
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
mmmHeatmap <- function (x = NULL,
                        group = 1,
                        xlab = "Area To",
                        ylab = "Area From",
                        area_names = NULL,
                        font_size_probs = 3,
                        font_nudge_probs = 0.15,
                        font_size_stderr = 2,
                        font_nudge_stderr = 0.15,
                        legend_name = "Estimate",
                        legend_position = "right") {

  # Check arguments ------------------------------------------------------------


  # Define x -------------------------------------------------------------------

  if (is.element("mmmFit", class(x))) {
    x <- x$results$movement_rate
    x_inds <- which(x$Group == group)
    x <- x[x_inds,]
  } else if (is.data.frame(x)) {
    # TODO: Check appropriate columns are present
    x_inds <- which(x$Group == group)
    x <- x[x_inds,]
  } else {
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
        label = sensibly_round(Estimate, 2),
        col = as.factor(ifelse(Estimate >= 0.5, 0, 1))),
      fontface = "plain",
      nudge_y = font_nudge_probs,
      size = font_size_probs) +
    # Add SE
    ggplot2::geom_text(
      mapping = ggplot2::aes(
        label = paste0("(", sensibly_round(SE, 3), ")"),
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

#' Lineplot Movement Rates from an 'mmmFit' Object
#'
#' @param x An object of class \code{mmmTMB}
#' @param x_name [character()] Axis name
#' @param y_name [character()] Axis name
#' @param x_breaks [numeric()] Axis breaks
#' @param x_labels [character()] Axis labels
#' @param x_labels_angle [numeric()] Axis labels angle
#' @param x_labels_hjust [numeric()] Axis labels hjust
#' @param area_names [character()] Vector of area names
#' @param class_names [character()] Vector of class names
#' @param error_alpha [numeric()] Error bar alpha
#' @param error_width [numeric()] Error bar width
#' @param legend_title [character()] Legend title
#' @param legend_position [character()] Legend position
#'
#' @return An object with class \code{gg} and \code{ggplot}
#'
#' @export
#'
#' @examples
#'
mmmLineplot <- function (x = NULL,
                      x_name = "Year",
                      y_name = "Annual Movement Rate",
                      x_breaks = NULL,
                      x_labels = NULL,
                      x_labels_angle = 0,
                      x_labels_hjust = 0,
                      area_names = NULL,
                      class_names = NULL,
                      error_alpha = 1,
                      error_width = 0.2,
                      legend_title = "Length Class",
                      legend_position = "right") {

  #---------------- Check arguments -------------------------------------------#


  #---------------- Extract movement results data -----------------------------#

  if (is.element("mmmFit", class(x))) {
    d <- x$results$movement_rate
  } else if (is.data.frame(x)) {
    # TODO: Check appropriate columns are present
    d <- x
  } else {
    warning("x not as expected")
  }

  #---------------- Convert columns to factor ---------------------------------#

  d$Area_From <- factor(
    d$Area_From,
    levels = sort(unique(d$Area_From)))
  d$Area_To <- factor(
    d$Area_To,
    levels = sort(unique(d$Area_To)))
  d$Class <- factor(
    d$Class,
    levels = sort(unique(d$Class)))

  #---------------- Create area labels ----------------------------------------#

  names(area_names) <- as.character(seq_along(area_names))

  #---------------- Construct geom object -------------------------------------#

  ggplot2::ggplot(data = d) +
    # ggplot2::geom_point(
    #   ggplot2::aes(x = Pattern_Time, y = Estimate, colour = Class)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        x = Pattern_Time,
        y = Estimate,
        colour = Class,
        ymin = Estimate - SE,
        ymax = Estimate + SE),
      alpha = error_alpha,
      width = error_width) +
    ggplot2::geom_line(
      ggplot2::aes(x = Pattern_Time, y = Estimate, colour = Class)) +
    ggplot2::scale_color_brewer(
      palette = "Blues",
      labels = class_names) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Area_From),
      cols = ggplot2::vars(Area_To),
      labeller = ggplot2::labeller(
        Area_From = area_names,
        Area_To = area_names),
      switch = "y"
    ) +
    ggplot2::scale_x_continuous(
      name = paste0("\n", x_name),
      breaks = x_breaks,
      labels = x_labels) +
    ggplot2::scale_y_continuous(
      name = paste0(y_name, "\n"),
      position = "right") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = x_labels_angle, hjust = x_labels_hjust),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(
        colour = "black",
        fill = "white",
        size = 0.75,
        linetype = "solid"),
      legend.position = legend_position) +
    ggplot2::labs(colour = legend_title)
}

#' Sensibly Round
#'
#' @param x [numeric()] Value to round
#' @param digits [numeric] Number of digits
#'
#' @return [character()]
#'
sensibly_round <- function (x, digits) {
  y <- character(length = length(x))
  ind_low <- which(x < 10^(-digits))
  ind_high <- which(x >= 10^(-digits))
  # Assign
  y[ind_low] <- paste0("<", sub(".", "", 10^(-digits)))
  y[ind_high] <- as.character(round(x[ind_high], digits))
  # Return
  y
}

#' Barplot Movement Probabilities from an 'mmmTMB' Object
#'
#' @param x An object of class \code{mmmTMB}
#' @param area_names [character()] Vector of area names
#'
#' @return An object with class \code{gg} and \code{ggplot}
#' @export
#'
#' @examples
#'
mmmBarplot <- function (x = NULL, area_names = NULL) {

  #---------------- Check arguments -------------------------------------------#


  #---------------- Extract movement results data -----------------------------#

  if (is.element("mmmFit", class(x))) {
    x <- x$results$movement_rate
    # x_inds <- which(x$Class == class)
    # x <- x[x_inds,]
  } else if (is.data.frame(x)) {
    # TODO: Check appropriate columns are present
    # x_inds <- which(x$Class == class)
    # x <- x[x_inds,]
  } else {
    warning("x not as expected")
  }

  #---------------- Convert columns to factor ---------------------------------#

  x$Area_From <- factor(
    x$Area_From,
    levels = sort(unique(x$Area_From), decreasing = FALSE))
  x$Area_To <- factor(
    x$Area_To,
    levels = sort(unique(x$Area_To)))
  x$Class <- factor(
    x$Class,
    levels = sort(unique(x$Class))
  )

  #---------------- Construct geom object -------------------------------------#

  ggplot2::ggplot(
    data = x,
    mapping = ggplot2::aes(x = Pattern_Time, y = Estimate, fill = Class)
    ) +
    ggplot2::geom_col(position = ggplot2::position_dodge()) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = Estimate - SE,
        ymax = Estimate + SE),
        width = 0.2,
        position = ggplot2::position_dodge(0.9)
    ) +
    ggplot2::scale_fill_brewer(palette = "Blues") +
    ggplot2::facet_grid(
      rows = ggplot2::vars(Area_From),
      cols = ggplot2::vars(Area_To),
      switch = "y"
    ) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(position = "right")

}
