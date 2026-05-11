#' @title Clean and Repair Oxylipin Name Encoding
#'
#' @description
#' Internal helper function to fix encoding artifacts in lipid names. It replaces
#' common broken characters (resulting from encoding mismatches) with a dot
#' and then maps those patterns back to their correct biochemical symbols
#' (e.g., \eqn{\alpha}, \eqn{\beta}, \eqn{\gamma}, \eqn{\Delta}).
#'
#' If unknown encoding artifacts are detected outside the internal database,
#' it alerts the user to manually correct them using proper Unicode.
#'
#' @param names_vec Character vector containing oxylipin names with potential
#' encoding issues.
#'
#' @return A character vector with corrected lipid names and restored special
#' characters.
#'
#' @details
#' This function is specifically tuned for oxylipin nomenclature used in
#' LC-MS/MS outputs where special characters often get corrupted
#' into symbols like broken encoding artifacts.
#'
#' @importFrom cli cli_alert_warning cli_bullets cli_alert_info
#'
#' @keywords internal
#' @export
clean_labels <- function(names_vec) {

  # Track which strings actually have encoding artifacts (excluding normal spaces)
  has_artifact <- grepl("[\u00c2\u00c3\u00b6\ufffd\\?]", names_vec)

  # Replace broken chars and spaces with a dot placeholder
  dotted_vec <- gsub("[? ]|\u00c2|\u00c3|\u00b6|\ufffd", ".", names_vec)

  # Dictionary of known corrections for the staRoxy default universe
  corrections <- c(
    "13-HOTrE-." = "13-HOTrE-\u03b3",
    "15-deoxy-.12,14-PGJ2" = "15-deoxy-\u039412,14-PGJ2",
    "2,3-dinor-8-iso-PGF2." = "2,3-dinor-8-iso-PGF2\u03b1",
    "15-deoxy-.12,14-PGD2" = "15-deoxy-\u039412,14-PGD2",
    "15-keto-PGF2." = "15-keto-PGF2\u03b1",
    "8-iso-15-keto-PGF2." = "8-iso-15-keto-PGF2\u03b2",
    "PGF3." = "PGF3\u03b1",
    "8-iso-PGF3." = "8-iso-PGF3\u03b1",
    "11\u03b2-13,14-dihydro-15-keto-PGF2." = "11\u03b2-13,14-dihydro-15-keto-PGF2\u03b1",
    "11\u03b2-PGF2." = "11\u03b2-PGF2\u03b1",
    "15-keto-PGF1." = "15-keto-PGF1\u03b1",
    "PGF2." = "PGF2\u03b1",
    "8,12-iso-iPF2.-VI" = "8,12-iso-iPF2\u03b1-VI",
    "13,14-dihydro-PGF2." = "13,14-dihydro-PGF2\u03b1",
    "PGF1." = "PGF1\u03b1",
    ".17-6-keto-PGF1." = "\u039417-6-keto-PGF1\u03b1",
    "6-keto-PGF1." = "6-keto-PGF1\u03b1",
    "6,15-diketo-13,14-dihydro-PGF1." = "6,15-diketo-13,14-dihydro-PGF1\u03b1"
  )

  # Apply corrections
  fixed <- dotted_vec
  for (err in names(corrections)) {
    fixed[fixed == err] <- corrections[err]
  }

  # Detect and warn about unresolved artifacts
  unresolved <- unique(names_vec[has_artifact & !(dotted_vec %in% names(corrections))])

  if (length(unresolved) > 0) {
    cli::cli_alert_warning("Potential encoding artifacts detected in the following oxylipins:")
    cli::cli_bullets(paste0("* ", unresolved))
    cli::cli_alert_info("If these lipids contain Greek letters, please rename them in your raw data using proper Unicode (e.g., \\u03b1 for Alpha, \\u03b2 for Beta).")
  }

  return(fixed)
}
