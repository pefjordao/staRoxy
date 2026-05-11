#' @title Get Default Oxylipin Reference Universe
#'
#' @description
#' Retrieves the default \code{staRoxy} reference universe, mapping oxylipins 
#' to their fatty acid precursors. This object is used natively by the 
#' \code{oxy_ora} function.
#'
#' @return A data frame with two columns: \code{Oxylipin} and \code{Precursor}.
#'
#' @examples
#' get_oxy_universe()
#'
#' @export
get_oxy_universe <- function() {
  # Call the internal data object saved via usethis::use_data()
  return(staRoxy::oxy_universe_default)
}
