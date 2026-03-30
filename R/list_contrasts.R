#' @title List Experimental Contrasts
#'
#' @description
#' Helper function to display all pairwise comparisons available in a
#' \code{limma_oxy} results object.
#'
#' @param stats_obj The list object returned by \code{limma_oxy}.
#'
#' @return Invisibly returns the contrast names. Prints to console.
#'
#' @importFrom cli cli_h1 cli_alert_info cli_alert_danger cli_li cli_rule
#'
#' @export
list_contrasts <- function(stats_obj) {
  if (is.null(stats_obj) || length(stats_obj) == 0)
    return(cli::cli_alert_danger("The stats object is empty or invalid."))
  contrast_names <- names(stats_obj); n <- length(contrast_names)
  cli::cli_h1("Experimental Contrasts")
  cli::cli_alert_info("Found {n} comparison{?s} in this analysis:")
  df_list <- data.frame(ID = 1:n, Contrast = contrast_names, Label = paste0("Contrast ", 1:n))
  for (i in 1:nrow(df_list)) cli::cli_li("{.strong {df_list$Label[i]}:} {df_list$Contrast[i]}")
  cli::cli_rule()
}
