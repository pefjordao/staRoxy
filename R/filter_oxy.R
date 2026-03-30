#' @title Filter Oxylipins Based on Missing Values
#'
#' @description
#' Filters oxylipins from the dataset according to three different strategies
#' for handling missing values (NAs). This step is essential to ensure
#' statistical robustness before linear modeling or multivariate analysis.
#'
#' @details
#' The function offers three filtering logic types:
#' \itemize{
#'   \item \strong{Type I (Default):} Retains oxylipins with at least \code{prop}
#'   valid values in \strong{at least one} experimental group. This strategy
#'   is ideal for preserving lipids that might be exclusive to a specific
#'   condition (e.g., biomarkers that only appear after stimulation).
#'   \item \strong{Type II:} Retains oxylipins with at least \code{prop} valid
#'   values across the \strong{entire dataset} (all samples combined).
#'   \item \strong{Type III (Strict):} Only retains oxylipins with \strong{no
#'   missing values} across all samples.
#' }
#'
#' @param obj A \code{staRoxy} object containing \code{data} (abundance matrix)
#' and \code{meta} (metadata).
#' @param type Character. The filtering strategy: \code{"I"}, \code{"II"},
#' or \code{"III"}. Default is \code{"I"}.
#' @param prop Numeric. The proportion threshold (0 to 1) for valid values.
#' Default is \code{0.5} (50%).
#'
#' @return Returns the updated \code{staRoxy} object with the filtered data
#' matrix and updated metadata information (\code{n_oxylipins}).
#'
#' @importFrom cli cli_alert_info cli_alert_success
#'
#' @examples
#' # Keep lipids present in at least 50% of samples in at least one group
#' data_filtered <- filter_oxy(data_oxy_pellet, type = "I", prop = 0.5)
#'
#' # Strict filtering: keep only lipids with 100% valid values
#' data_strict <- filter_oxy(data_oxy_pellet, type = "III")
#'
#' @export
filter_oxy <- function(obj, type = "I", prop = 0.5) {

  # Data Extraction
  # Extract data matrix and metadata for filtering logic
  d <- obj$data
  m <- obj$meta

  # Filter Logic Selection
  # Define 'keep' vector based on the selected filtering strategy
  if (type == "I") {
    # Type I: Keep if valid values are >= prop in AT LEAST one experimental group
    msg <- paste0("Type I: >", prop * 100, "% valid values in AT LEAST one group.")

    keep <- apply(d, 1, function(x) {
      any(tapply(x, m$group, function(g) sum(!is.na(g)) / length(g) >= prop))
    })

  } else if (type == "II") {
    # Type II: Keep if valid values are >= prop across the entire dataset (all samples)
    msg <- paste0("Type II: >", prop * 100, "% valid values across ALL samples.")
    keep <- rowSums(!is.na(d)) / ncol(d) >= prop

  } else {
    # Type III: Strict filtering. Only keep oxylipins with NO missing values
    msg <- "Type III: Strict mode. NO missing values allowed."
    keep <- rowSums(is.na(d)) == 0
  }

  # Object Update
  # Subset the data matrix and update the metadata count
  obj$data <- d[keep, , drop = FALSE]
  obj$info$n_oxylipins <- nrow(obj$data)

  # Console Output
  cli::cli_alert_info(msg)
  cli::cli_alert_success("Filtering complete: {obj$info$n_oxylipins} oxylipins remaining.")

  return(obj)
}
