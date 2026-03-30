#' @title Data Distribution and Experimental Design Advisor
#'
#' @description
#' Performs a Quality Control (QC) diagnostic to evaluate distribution
#' symmetry (skewness), group sample sizes, and experimental design balance.
#' It ensures the dataset meets the assumptions required for robust linear modeling.
#'
#' @details
#' The function assesses three main pillars:
#' \itemize{
#'   \item \strong{Skewness (\eqn{\gamma_1}):} Measures asymmetry. If \eqn{|\gamma_1| > 1.2},
#'   the data is likely in raw scale and log2 transformation is recommended.
#'   \item \strong{Sample Size:} Verifies if \eqn{n \ge 3} per group to ensure
#'   sufficient power for variance estimation.
#'   \item \strong{Experimental Balance (R):} Calculates the ratio between the
#'   largest and smallest group. If \eqn{R > 2.0}, it warns about potential
#'   heteroscedasticity risks.
#' }
#'
#' @param obj A \code{staRoxy} object containing \code{data} (abundance matrix)
#' and \code{meta} (metadata).
#'
#' @return Invisibly returns a list containing:
#' \item{skewness}{The calculated Fisher-Pearson skewness coefficient (\eqn{\gamma_1}).}
#' \item{sample_sizes}{A named numeric vector with counts per experimental group.}
#'
#' @importFrom cli cli_h1 cli_alert_info cli_alert_warning cli_alert_danger cli_alert_success cli_rule
#' @importFrom stats sd
#'
#' @examples
#' # Assess the distribution and design of the pellet example dataset
#' check_data_distribution(data_oxy_pellet)
#'
#' @export
check_data_distribution <- function(obj) {
  cli::cli_h1("Model Advisor Report")

  # 1. Skewness and Scale Analysis
  # Isolate numeric values for distribution testing (excluding NAs and zeros)
  vals <- as.vector(obj$data)
  vals <- vals[vals > 0 & !is.na(vals)]

  # Internal helper for skewness calculation (Fisher-Pearson)
  calc_skew <- function(x) {
    n <- length(x)
    mu <- mean(x)
    m3 <- sum((x - mu)^3) / n
    s3 <- stats::sd(x)^3
    return(m3 / s3)
  }

  gamma1 <- calc_skew(vals)

  cli::cli_alert_info("Testing distribution symmetry using Pearson's skewness coefficient (\u03b31).")

  if (abs(gamma1) > 1.2) {
    cli::cli_alert_warning("Scale Check: High skewness detected (\u03b31 = {round(gamma1, 2)}).")
    cli::cli_alert_danger("Advice: Data appears to be in RAW scale. Log2 transformation is highly recommended.")
  } else {
    cli::cli_alert_info("Scale Check: Distribution is reasonably symmetric (\u03b31 = {round(gamma1, 2)}).")
    cli::cli_alert_success("Data scale looks appropriate for linear modeling (likely log-transformed).")
  }

  # 2. Sample Size and Power Check
  n_counts <- tapply(obj$meta$sample, obj$meta$group, length)

  cli::cli_alert_info("Checking sample sizes per group:")

  for (g in names(n_counts)) {
    n <- n_counts[g]
    if (n < 3) {
      cli::cli_alert_danger("Group '{g}': n = {n}. Insufficient statistical power for robust variance estimation!")
    } else {
      cli::cli_alert_success("Group '{g}': n = {n}. Sample size is adequate.")
    }
  }

  # 3. Experimental Design Balance
  if (length(n_counts) > 1) {
    ratio <- max(n_counts) / min(n_counts)
    cli::cli_alert_info("Imbalance Ratio (R): {round(ratio, 2)}x.")

    if (ratio > 2) {
      cli::cli_alert_warning("Unbalanced design: R > 2.0. Risk of heteroscedasticity (unequal variances).")
      cli::cli_alert_info("Consider using 'weights' in limma or robust estimation methods.")
    } else {
      cli::cli_alert_success("Experimental design is balanced (R <= 2.0). Assumptions are likely met.")
    }
  }

  cli::cli_rule()
  return(invisible(list(skewness = gamma1, sample_sizes = n_counts)))
}
