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
#' @importFrom cli cli_alert_info
#'
#' @examples
#' # List all lipids currently in the pellet dataset
#' get_analyzed_oxy(data_oxy_pellet)
#'
#' @export
get_analyzed_oxy <- function(obj) {
  lipid_list <- rownames(obj$data)

  cli::cli_alert_info("Found {length(lipid_list)} oxylipins detected in the dataset.")

  return(lipid_list)
}
