#' @title Automatic Correlation Method Selection
#'
#' @description
#' Internal helper function to recommend a correlation method (Pearson or Spearman)
#' based on normality tests (Shapiro-Wilk) and outlier detection (IQR method).
#'
#' @param x Numeric vector.
#' @param y Numeric vector.
#' @param silent Logical. If \code{FALSE}, prints the decision diagnostic to the console.
#' Default is \code{TRUE}.
#'
#' @return A character string: \code{"pearson"} if data is normal and has no outliers,
#' otherwise \code{"spearman"}. The result includes \code{"is_normal"} and
#' \code{"has_out"} as logical attributes.
#'
#' @importFrom stats shapiro.test quantile
#' @importFrom cli cli_alert_info
#'
#' @keywords internal
check_cor_method <- function(x, y, silent = TRUE) {
  # 1. Normality check (Shapiro-Wilk)
  sh_p      <- c(stats::shapiro.test(x)$p.value, stats::shapiro.test(y)$p.value)
  is_normal <- all(sh_p > 0.05)

  # 2. Outlier detection (IQR method)
  is_outlier <- function(v) {
    q           <- stats::quantile(v, probs = c(0.25, 0.75), na.rm = TRUE)
    lower_bound <- q[1] - 1.5 * diff(q)
    upper_bound <- q[2] + 1.5 * diff(q)
    any(v < lower_bound | v > upper_bound, na.rm = TRUE)
  }
  has_out <- is_outlier(x) | is_outlier(y)

  # 3. Decision logic
  rec <- if (is_normal && !has_out) "pearson" else "spearman"

  # Store diagnostics as attributes
  attr(rec, "is_normal") <- is_normal
  attr(rec, "has_out")   <- has_out

  if (!silent) {
    cli::cli_alert_info("Method: {toupper(rec)} (Normal={is_normal}, Outliers={has_out})")
  }

  return(rec)
}
