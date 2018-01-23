#' Pipe operator
#'
#' See \code{magrittr::\link[magrittr]{\%>\%}} for details.
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL

#' Function to write to console
#'
#' @name printf
#' @rdname printf
#' @export
#' @usage printf("Hello world")

printf <- function(...) cat(sprintf(...), "\n")
