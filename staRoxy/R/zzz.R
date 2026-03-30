#' @importFrom dplyr %>%
NULL

.onAttach <- function(libname, pkgname) {

  pkg_version <- utils::packageVersion("staRoxy")

  cli::cli_rule(
    left = paste0("{.strong staRoxy} ver. ", pkg_version),
    right = "{.field [University of São Paulo - USP]}"
  )

  cli::cli_alert_info("Oxylipidomics Abundance Data Analysis Pipeline")
  cli::cli_text("Developed by: {.emph Pedro Henrique F. Jordão} & {.emph Eduardo M. Reis}")
  cli::cli_text("Official Repository: {.url https://github.com/pefjordao/staRoxy}")
  cli::cli_rule()
}
