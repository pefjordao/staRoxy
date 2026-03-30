#' @title Classify Oxylipins by Fatty Acid Precursors
#'
#' @description
#' Maps detected oxylipins to their biological fatty acid precursors (e.g., ARA,
#' EPA, DHA, LA) using an internal biochemical dictionary. It provides either
#' a detailed presence/absence matrix per group or a summarized count of
#' unique species per precursor.
#'
#' @details
#' The classification relies on a built-in mapping table of common oxylipins.
#' A lipid is considered "Present" in a group if it has at least one finite
#' (non-NA, non-Inf) value across the samples of that group.
#'
#' \strong{Note on Nomenclature:}
#' Lipid identifiers must match the package database (see vignette). Use \code{show_database_oxy()} to view the supported names and
#' their respective precursors.
#'
#' @param obj A \code{staRoxy} object containing \code{data} (abundance matrix)
#' and \code{meta} (metadata).
#' @param return_mode Character. Use \code{"summary"} (default) for counts of unique
#' species per precursor/group, or \code{"species"} for a binary presence (1)
#' vs. absence (0) matrix for each lipid.
#'
#' @return A data frame containing either the classification summary or the
#' species-level binary matrix.
#'
#' @importFrom dplyr distinct mutate left_join select across group_by summarise all_of arrange relocate n
#' @importFrom tidyr replace_na
#' @importFrom cli cli_h1 cli_alert_info
#'
#' @examples
#' # Get a summary of unique species per precursor
#' classify_oxy(data_oxy_pellet, return_mode = "summary")
#'
#' # Get the binary presence/absence matrix
#' classify_oxy(data_oxy_pellet, return_mode = "species")
#'
#' @export
classify_oxy <- function(obj, return_mode = "summary") {

  # 1. Internal Mapping Table
  univ_map <- data.frame(
    Lipid = c("tetranor-12-HETE", "12-HHTrE", "13-oxo-ODE", "9-HOTrE", "9-oxo-ODE", "13-HOTrE-γ", "13-HOTrE", "9-HODE", "10-HODE", "12-HODE", "13-HODE", "9(10)-EpOME", "12(13)-EpOME", "9-HOA", "10-HOA", "9,10-DiHOME", "12,13-DiHOME", "15-deoxy-Δ12,14-PGJ2", "15-oxo-ETE", "5-HEPE", "9-HEPE", "8-HEPE", "11-HEPE", "12-HEPE", "15-HEPE", "18-HEPE", "14(15)-EpETE", "17(18)-EpETE", "5-HETE", "8-HETE", "9-HETE", "11-HETE", "12-HETE", "14 / 15-HETE", "16-HETE", "17-HETE", "18-HETE", "19-HETE", "20-HETE", "5(6)-EpETrE", "8(9)-EpETrE", "11(12)-EpETrE", "14(15)-EpETrE", "15-oxo-EDE", "5-HETrE", "15-HETrE", "2,3-dinor-8-iso-PGF2α", "10-Nitrooleate", "9-Nitrooleate", "tetranor-PGDM", "PGA2", "PGB2", "PGJ2", "15-deoxy-Δ12,14-PGD2", "12-oxo-LTB4", "LTB4", "6-trans-LTB4", "6-trans-12-epi-LTB4", "12-epi-LTB4", "8,15-DiHETE", "14,15-DiHETE", "8,9-DiHETrE", "17-oxo-DHA", "2,3-dinor-TXB2", "4-HDoHE", "7-HDoHE", "8-HDoHE", "10-HDoHE", "11-HDoHE", "13-HDoHE", "14-HDoHE", "16 / 17-HDoHE", "19-HDoHE", "20-HDoHE", "7(8)-EpDPE", "10(11)-EpDPE", "13(14)-EpDPE", "16(17)-EpDPE", "19(20)-EpDPE", "17-oxo-DPA", "PGE3", "15-keto-PGE2", "LXA5", "RvE1", "epi-LXA4", "15-epi-LXA4", "LXB4", "20-hydroxy-LTB4", "PGD2", "13,14-dihydro-15-keto-PGD2", "PGE2", "13,14-dihydro-15-keto-PGE2", "15-keto-PGF2α", "8-iso-15-keto-PGF2β", "PGF3α", "8-iso-PGF3α", "15-F2t-Isoprostane", "11β-13,14-dihydro-15-keto-PGF2α", "11β-PGF2α", "15-keto-PGF1α", "8-iso-PGE1", "PGD1", "PGF2α", "8,12-iso-iPF2α-VI", "13,14-dihydro-PGF2α", "PGF1α", "PDX", "Maresin 1", "Maresin 2", "RvD5", "19,20-DiHDPA", "20-carboxy-LTB4", "TXB3", "20-hydroxy-PGE2", "6-keto-PGE1", "Δ17-6-keto-PGF1α", "TXB2", "6-keto-PGF1α", "6,15-diketo-13,14-dihydro-PGF1α", "TXB1", "RvD1", "RvD2", "RvD3", "16,16-dimethyl-PGD2", "1a,1b-dihomo-PGE2"),
    Precursor = c("ARA (20:4 n-6)", "ARA (20:4 n-6)", "LA (18:2 n-6)", "ALA (18:3 n-3)", "LA (18:2 n-6)", "GLA (18:3 n-6)", "ALA (18:3 n-3)", "LA (18:2 n-6)", "LA (18:2 n-6)", "LA (18:2 n-6)", "LA (18:2 n-6)", "LA (18:2 n-6)", "LA (18:2 n-6)", "OA (18:1 n-9)", "OA (18:1 n-9)", "LA (18:2 n-6)", "LA (18:2 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "EDA (20:2 n-6)", "MA (20:3 n-9)", "DGLA (20:3 n-6)", "ARA (20:4 n-6)", "OA (18:1 n-9)", "OA (18:1 n-9)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "DHA (22:6 n-3)", "ARA (20:4 n-6)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DPA (22:5 n-3)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "ARA (20:4 n-6)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "DGLA (20:3 n-6)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "ARA (20:4 n-6)", "ARA (20:4 n-6)")
  ) %>% dplyr::distinct(Lipid, .keep_all = TRUE)

  # Clean lipid names in mapping table
  univ_map$Lipid <- trimws(univ_map$Lipid)

  # 2. Data Initialization
  data_raw      <- obj$data
  lipid_names   <- trimws(rownames(data_raw))
  groups        <- obj$meta$group
  unique_groups <- unique(groups)

  # 3. Presence Matrix Construction
  # Initialize binary matrix for lipid presence per group
  presence_mat <- data.frame(Lipid = lipid_names, stringsAsFactors = FALSE)

  for (grp in unique_groups) {
    cols_in_group <- which(groups == grp)

    # Present if at least one sample in the group has a finite value (detected)
    is_present <- apply(data_raw[, cols_in_group, drop = FALSE], 1, function(row) {
      any(is.finite(row))
    })

    presence_mat[[grp]] <- as.numeric(is_present)
  }

  # 4. Annotation Join
  # Link presence data with precursors
  species_df <- presence_mat %>%
    dplyr::left_join(univ_map, by = "Lipid") %>%
    dplyr::mutate(Precursor = tidyr::replace_na(Precursor, "Other/Not found")) %>%
    # Ensure Precursor is the second column
    dplyr::relocate(Precursor, .after = Lipid) %>%
    dplyr::arrange(Precursor, Lipid)

  # 5. Return Mode: Species Matrix
  if (return_mode == "species") {
    cli::cli_h1("Species Presence/Absence Matrix")
    cli::cli_alert_info("1 = Presence (detected in group) | 0 = Absence")

    print(as.data.frame(species_df), row.names = FALSE)
    return(invisible(species_df))
  }

  # 6. Return Mode: Summary
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
