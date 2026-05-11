#' @importFrom dplyr %>%
NULL

utils::globalVariables(c(
  ".", "where",

  "Abund", "adj.P.Val", "cv_perc", "Detected", "Enrich_Score",
  "err_max", "err_min", "Feature", "G1", "G2", "Group", "group",
  "has_na", "intensity", "L1", "L2", "logFC", "m", "mean_int",
  "n_na", "n_valid", "na_count", "na_prop", "Observed", "Oxylipin",
  "oxylipin", "p_dep", "p_enr", "PC1", "PC2", "PCoA1", "PCoA2",
  "plot_group", "plot_order", "Precursor", "score", "sem", "SEM",
  "sig_label", "Total_Group_Detected", "Total_Unique_Detected",
  "Universe_M", "Universe_N", "Value", "Weight"
))

.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion("staRoxy")

  msg <- c(
    cli::rule(
      left = cli::format_inline("{.strong {.col_cyan staRoxy v{pkg_version}}}"),
      right = cli::format_inline("{.col_grey University of S\u00e3o Paulo - USP}")
    ),
    cli::format_inline("{cli::col_blue(cli::symbol$info)} Oxylipidomics Abundance Data Analysis Pipeline"),
    cli::format_inline("Developed by: {.strong Pedro Henrique F. Jord\u00e3o & Eduardo M. Reis}"),
    cli::format_inline("Official Repository: {.url https://github.com/pefjordao/staRoxy}"),
    cli::rule()
  )

  packageStartupMessage(paste(msg, collapse = "\n"))
}
