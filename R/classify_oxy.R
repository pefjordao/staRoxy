#' @title Classify Oxylipins by Fatty Acid Precursors
#'
#' @description
#' Maps detected oxylipins to their biological fatty acid precursors (e.g., ARA,
#' EPA, DHA, LA). It provides either a detailed presence/absence matrix per
#' group or a summarized count of unique species per precursor.
#'
#' @details
#' The classification relies on a reference universe mapping table. By default,
#' it uses the internal \code{staRoxy} database. A oxylipin is considered "Present"
#' in a group if it has at least one finite (non-NA, non-Inf) value across the
#' samples of that group.
#'
#' @param obj A \code{staRoxy} object.
#' @param universe A data frame mapping oxylipins to precursors. Defaults to
#' the output of \code{get_oxy_universe()}.
#' @param return_mode Character. Use \code{"summary"} (default) for counts of unique
#' species per precursor/group, or \code{"species"} for a binary presence (1)
#' vs. absence (0) matrix for each Oxylipin.
#'
#' @return A data frame containing either the classification summary or the
#' species-level binary matrix.
#'
#' @examples
#' \dontshow{
#' staRoxy_object <- read_oxy(data_oxy_lps_pellet, metadata_oxy_lps_pellet)
#' staRoxy_object <- filter_oxy(staRoxy_object)
#' staRoxy_object <- transform_oxy(staRoxy_object)
#' }
#'
#' # Standard summary using the default universe
#' classification_summary <- classify_oxy(staRoxy_object, return_mode = "summary")
#'
#' # Classification using a custom edited universe
#' my_custom_univ <- edit_oxy_universe(remove = c("PGE2", "PGD2"))
#' custom_classification <- classify_oxy(
#'   obj = staRoxy_object,
#'   universe = my_custom_univ,
#'   return_mode = "species"
#' )
#'
#' @importFrom dplyr distinct mutate left_join select across group_by summarise all_of arrange relocate n
#' @importFrom tidyr replace_na
#' @importFrom cli cli_h1 cli_alert_info
#'
#' @export
classify_oxy <- function(obj, universe = get_oxy_universe(), return_mode = "summary") {

  # Input Validation (Strict staRoxy object enforcement)
  if (!is.list(obj) || !all(c("data", "meta") %in% names(obj))) {
    stop("Input must be a valid staRoxy object (a list containing 'data' and 'meta').")
  }

  if (is.null(obj$meta$group)) {
    stop("The 'meta' element of the staRoxy object must contain a 'group' column.")
  }

  # Load Reference Universe from argument
  univ_map <- universe
  univ_map$Oxylipin <- trimws(univ_map$Oxylipin)

  # Data Initialization (Direct extraction)
  data_raw <- as.matrix(obj$data)
  groups <- obj$meta$group

  oxylipin_names <- clean_labels(trimws(rownames(data_raw)))
  rownames(data_raw) <- oxylipin_names

  unique_groups <- unique(groups)

  # Presence Matrix Construction
  # Initialize binary matrix for Oxylipin presence per group
  presence_mat <- data.frame(Oxylipin = oxylipin_names, stringsAsFactors = FALSE)

  for (grp in unique_groups) {
    cols_in_group <- which(groups == grp)

    # Present if at least one sample in the group has a finite value (detected)
    is_present <- apply(data_raw[, cols_in_group, drop = FALSE], 1, function(row) {
      any(is.finite(row))
    })

    presence_mat[[grp]] <- as.numeric(is_present)
  }

  # Annotation Join
  # Link presence data with precursors
  species_df <- presence_mat %>%
    dplyr::left_join(univ_map, by = "Oxylipin") %>%
    dplyr::mutate(Precursor = tidyr::replace_na(Precursor, "Other/Not found")) %>%
    # Ensure Precursor is the second column
    dplyr::relocate(Precursor, .after = Oxylipin) %>%
    dplyr::arrange(Precursor, Oxylipin)

  # Return Mode: Species Matrix
  if (return_mode == "species") {
    cli::cli_h1("Species Presence/Absence Matrix")
    cli::cli_alert_info("1 = Presence (detected in group) | 0 = Absence")

    print(as.data.frame(species_df), row.names = FALSE)
    return(invisible(species_df))
  }

  # Return Mode: Summary
  if (return_mode == "summary") {
    # Count unique species per precursor for each group
    summary_df <- species_df %>%
      dplyr::group_by(Precursor) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(unique_groups), sum), .groups = "drop")

    # Calculate global uniqueness (detected in at least one group)
    total_unique <- species_df %>%
      dplyr::mutate(Detected = rowSums(dplyr::select(., dplyr::all_of(unique_groups))) > 0) %>%
      dplyr::filter(Detected == TRUE) %>%
      dplyr::group_by(Precursor) %>%
      dplyr::summarise(Total_Unique_Detected = dplyr::n(), .groups = "drop")

    final_stats <- summary_df %>%
      dplyr::left_join(total_unique, by = "Precursor") %>%
      dplyr::mutate(Total_Unique_Detected = tidyr::replace_na(Total_Unique_Detected, 0)) %>%
      as.data.frame()

    cli::cli_h1("Oxylipin Classification Summary")
    cli::cli_alert_info("Counts represent unique species detected per precursor/group.")

    print(final_stats, row.names = FALSE)
    return(invisible(final_stats))
  }

  stop("Invalid return_mode. Choose 'summary' or 'species'.")
}
