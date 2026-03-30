#' @title Detect Group-Exclusive Oxylipins
#'
#' @description
#' Identifies oxylipins that are detected in only one experimental group
#' and are completely absent (all NAs) in all other groups. This is a
#' crucial exploratory tool for discovering "all-or-nothing" biomarkers.
#'
#' @details
#' Exclusivity is determined strictly: a lipid must have at least one valid
#' numerical value in the target group and be 100% missing (all NAs) in every
#' other group in the dataset.
#'
#' @param obj A \code{staRoxy} object.
#'
#' @return Invisibly returns a named list where each element is a character
#' vector of exclusive oxylipins for a specific group. A summary report
#' is printed to the console.
#'
#' @importFrom cli cli_h1 cli_alert_info cli_alert_success cli_alert_danger cli_rule
#'
#' @examples
#' # Identify oxylipins unique to Scraper or Trypsin groups
#' get_exclusives_oxy(data_oxy_pellet)
#'
#' @export
get_exclusives_oxy <- function(obj) {
  # 1. Setup and Validation
  d <- obj$data
  m <- obj$meta
  grps <- levels(as.factor(m$group))

  cli::cli_h1("Exclusive Oxylipins Detection")

  if (length(grps) < 2) {
    return(cli::cli_alert_danger("Exclusive analysis requires at least 2 experimental groups."))
  }

  # 2. NA Mapping
  # Identify which lipids are completely missing (all NA) per group
  get_na_status <- function(grp) {
    apply(d[, m$group == grp, drop = FALSE], 1, function(x) all(is.na(x)))
  }

  na_map <- sapply(grps, get_na_status)
  exclusive_list <- list()

  # 3. Exclusivity Logic
  for (g in grps) {
    # A lipid is exclusive to group 'g' if:
    # It is NOT all NA in group 'g' AND it IS all NA in every other group
    others_na <- rowSums(na_map[, colnames(na_map) != g, drop = FALSE]) == (length(grps) - 1)
    is_present_in_g <- !na_map[, g]

    excl <- rownames(d)[is_present_in_g & others_na]
    exclusive_list[[g]] <- excl

    # 4. Console Reporting
    if (length(excl) > 0) {
      cli::cli_alert_success("Group '{g}': {length(excl)} exclusive oxylipin{?s} found.")
      cat(paste0("  - ", excl, collapse = "\n"), "\n")
    } else {
      cli::cli_alert_info("Group '{g}': No exclusive oxylipins found.")
    }
  }

  cli::cli_rule()
  return(invisible(exclusive_list))
}
