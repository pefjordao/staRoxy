#' @title Interpret Correlation Strength and Direction
#'
#' @description
#' Internal helper function that converts a correlation coefficient (\eqn{r})
#' and its P-value (\eqn{p}) into a descriptive character string.
#'
#' @details
#' The interpretation follows a standard scale for the absolute value of \eqn{r}:
#' \itemize{
#'   \item \strong{0.00 - 0.19:} Very weak
#'   \item \strong{0.20 - 0.39:} Weak
#'   \item \strong{0.40 - 0.59:} Medium
#'   \item \strong{0.60 - 0.79:} Strong
#'   \item \strong{0.80 - 1.00:} Very strong
#' }
#' If \eqn{P \ge 0.05}, it returns "No significant correlation" regardless
#' of the \eqn{R} value.
#'
#' @param r Numeric. The correlation coefficient (\eqn{R \in [-1, 1]}).
#' @param p Numeric. The p-value associated with the correlation test.
#'
#' @return A character string describing the correlation (e.g., "Strong
#' positive correlation").
#'
#' @keywords internal
get_cor_label <- function(r, p) {
  if (is.na(p) || p >= 0.05) return("No significant correlation")
  abs_r <- abs(r)
  labels <- c("None", "Very weak", "Weak", "Medium", "Strong", "Very strong")
  breaks <- c(-Inf, 0, 0.199, 0.399, 0.599, 0.799, 1)
  tag <- labels[cut(abs_r, breaks = breaks, labels = FALSE)]
  direction <- if (r > 0) "positive" else "negative"
  return(paste(tag, direction, "correlation"))
}
