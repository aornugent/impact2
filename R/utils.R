#' Function to write to console
#'
#' @name printf
#' @rdname printf
#' @export
#' @usage printf("Hello world")

printf <- function(...) cat(sprintf(...), "\n")

#' Load dplyr
#'
#' dplyr is used extensively in this package
#'
#' @name dplyr
#' @import dplyr
NULL

#' Basic ggplot theme
#'
#' Set as default on package load.
#'
#' @import ggplot2

.onLoad <- function(libname, pkgname) {

  theme_impact <-
    theme_minimal() +
      theme(
        panel.grid = element_blank(),
        #panel.border = element_rect(fill = NA),
        axis.line.x = element_line(
          color="black",
          size = 0.5),
        axis.line.y = element_line(
          color="black",
          size = 0.5),
        axis.ticks = element_line(),
        panel.spacing = unit(1, "lines"),
        plot.title = element_text(
          margin = margin(t = 10, b = 20),
          hjust = 0),
        legend.position = "bottom",
        aspect.ratio = 1)

  theme_set(theme_impact)
}

#' expand.grid for dataframes
#'
#' @param ... a dataframe
#' @usage df_long <- expand.grid.df(df)

expand_grid <- function(...) {
  Reduce(function(...) merge(..., by=NULL), list(...))
}
