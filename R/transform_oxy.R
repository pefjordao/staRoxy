#' @title Log2 Transformation and Variance Filtering
#'
#' @description
#' Applies a log2 transformation to the oxylipin abundance matrix and removes
#' non-informative features (those with zero variance or entirely missing values).
#'
#' @details
#' Log transformation is a critical step in oxylipidomics to stabilize variance
#' and ensure that data meet the normality assumptions required for linear
#' modeling (e.g., in \code{limma_oxy}).
#'
#' Additionally, this function performs a secondary quality control check:
#' \itemize{
#'   \item \strong{Log2 Scaling:} Converts raw abundance to a logarithmic scale.
#'   \item \strong{Variance Filter:} Identifies and removes oxylipins that
#'   show no variation across all samples (constant values) or are entirely
#'   composed of NAs, as these features provide no statistical information
#'   and could cause errors in multivariate analyses.
#' }
#'
#' @param obj A \code{staRoxy} object.
#'
#' @return Returns the updated \code{staRoxy} object with log2-transformed
#' abundances and the filtered feature set.
#'
#' @importFrom stats var
#' @importFrom cli cli_alert_success cli_alert_warning cli_alert_info
#'
#' @examples
#' # Apply log2 transformation to the filtered dataset
#' staRoxy_object <- transform_oxy(staRoxy_object)
#'
#' @export
transform_oxy <- function(obj) {

  # Log2 Transformation
  # Apply log2 transformation to the data matrix
  log_d <- log2(obj$data)

  # Variance Check
  # Calculate variance per row (oxylipin) to identify non-informative profiles
  vars <- apply(log_d, 1, var, na.rm = TRUE)

  # Keep only oxylipins with valid variance (not NA and greater than zero)
  keep <- !is.na(vars) & vars > 0

  # Track removed features for reporting
  removed_names <- rownames(log_d)[!keep]
  removed_count <- length(removed_names)

  # Object Update
  # Update the data matrix and report the status
  obj$data <- log_d[keep, , drop = FALSE]
  cli::cli_alert_success("Log2 transformation complete.")

  # Feature Removal Reporting
  # Detailed console output for oxylipins that didn't pass the variance check
  if (removed_count > 0) {
    cli::cli_alert_warning("Removed {removed_count} oxylipin{?s} (zero variance or all NAs):")
    cat(paste0("  - ", removed_names, collapse = "\n"), "\n")
  } else {
    cli::cli_alert_info("All oxylipins retained: no zero variance or full NA profiles detected.")
  }

  return(obj)
}
