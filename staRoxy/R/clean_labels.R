#' @title Clean and Repair Oxylipin Name Encoding
#'
#' @description
#' Internal helper function to fix encoding artifacts in lipid names. It replaces
#' common broken characters (resulting from encoding mismatches) with a dot
#' and then maps those patterns back to their correct biochemical symbols
#' (e.g., alpha, beta, gamma, delta).
#'
#' @param names_vec Character vector containing oxylipin names with potential
#' encoding issues.
#'
#' @return A character vector with corrected lipid names and restored special
#' characters (α, β, γ, Δ).
#'
#' @details
#' This function is specifically tuned for oxylipin nomenclature used in
#' LC-MS/MS outputs where characters like 'Δ' or 'α' often get corrupted
#' into symbols like 'Â', 'Ã', '?', or '¶'.
#'
#' @keywords internal
clean_labels <- function(names_vec) {
  fixed <- gsub("[?]|Â|Ã|¶", ".", names_vec)
  corrections <- c(
    "13-HOTrE-." = "13-HOTrE-γ",
    "15-deoxy-.12,14-PGJ2" = "15-deoxy-Δ12,14-PGJ2",
    "2,3-dinor-8-iso-PGF2." = "2,3-dinor-8-iso-PGF2α",
    "15-deoxy-.12,14-PGD2" = "15-deoxy-Δ12,14-PGD2",
    "15-keto-PGF2." = "15-keto-PGF2α",
    "8-iso-15-keto-PGF2." = "8-iso-15-keto-PGF2β",
    "PGF3." = "PGF3α",
    "8-iso-PGF3." = "8-iso-PGF3α",
    "11β-13,14-dihydro-15-keto-PGF2." = "11β-13,14-dihydro-15-keto-PGF2α",
    "11β-PGF2." = "11β-PGF2α",
    "15-keto-PGF1." = "15-keto-PGF1α",
    "PGF2." = "PGF2α",
    "8,12-iso-iPF2.-VI" = "8,12-iso-iPF2α-VI",
    "13,14-dihydro-PGF2." = "13,14-dihydro-PGF2α",
    "PGF1." = "PGF1α",
    ".17-6-keto-PGF1." = "Δ17-6-keto-PGF1α",
    "6-keto-PGF1." = "6-keto-PGF1α",
    "6,15-diketo-13,14-dihydro-PGF1." = "6,15-diketo-13,14-dihydro-PGF1α")
  for (err in names(corrections)) fixed[fixed == err] <- corrections[err]
  return(fixed)
}
