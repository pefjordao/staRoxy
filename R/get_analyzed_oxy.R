#' @title Retrieve Names of Analyzed Oxylipins
#'
#' @description
#' A simple helper function to extract the names of all oxylipin features
#' currently retained in the dataset. This reflects the state of the data
#' after any filtering or subsetting steps.
#'
#' @param obj A \code{staRoxy} object.
#'
#' @return A character vector containing the names of the lipids (row names
#' of the data matrix).
#'
#' @examples
#' \dontshow{
#' staRoxy_object <- read_oxy(data_oxy_lps_pellet, metadata_oxy_lps_pellet)
#' staRoxy_object <- filter_oxy(staRoxy_object)
#' staRoxy_object <- transform_oxy(staRoxy_object)
#' }
#'
#' # Retrieve the names of the oxylipins currently in the object
#' get_analyzed_oxy(staRoxy_object)
#'
#' @importFrom cli cli_alert_info
#'
#' @export
get_analyzed_oxy <- function(obj) {

  # Input Validation
  if (!is.list(obj) || !all(c("data", "meta") %in% names(obj))) {
    cli::cli_abort("Input must be a valid staRoxy object (a list containing 'data' and 'meta').")
  }

  # Data Extraction
  raw_data <- obj$data
  lipid_list <- rownames(raw_data)

  # Safety Check for rownames
  if (is.null(lipid_list)) {
    cli::cli_abort("Error: The provided object has no row names (rownames).")
  }

  # Console Output & Return
  cli::cli_alert_info("Found {length(lipid_list)} oxylipins detected in the dataset.")

  return(lipid_list)
}
