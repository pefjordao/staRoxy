#' @title Compute Descriptive Statistics for Oxylipins
#'
#' @description
#' Calculates comprehensive summary statistics for each lipid feature across
#' experimental groups, including counts (valid/NA), central tendency (mean, median),
#' and dispersion (SD, SEM, CV%).
#'
#' @details
#' The function transforms the abundance matrix into a long format and merges it
#' with metadata to compute statistics per group. The Coefficient of Variation (CV%)
#' is calculated as \eqn{(SD / Mean) * 100}.
#'
#' Sorting and filtering options allow for quick identification of top-abundance
#' features or those with high analytical variability.
#'
#' @param obj A \code{staRoxy} object containing \code{data} (abundance matrix)
#' and \code{meta} (metadata).
#' @param group Character vector. Specific experimental groups to include.
#' If \code{NULL} (default), all groups are processed.
#' @param n_top Integer. If provided, returns only the top N features based on
#' the highest maximum mean value across groups.
#' @param sort_by Character. Metric to sort the results in descending order.
#' Options: \code{"cv"}, \code{"mean"}, \code{"median"}, \code{"sd"},
#' \code{"sem"}, \code{"na"}, \code{"n"}, or \code{"feature"} (alphabetical).
#'
#' @return Invisibly returns a data frame with the computed statistics.
#' Results are also printed to the console as a formatted table.
#'
#' @examples
#' \dontshow{
#' staRoxy_object <- read_oxy(data_oxy_lps_pellet, metadata_oxy_lps_pellet)
#' staRoxy_object <- filter_oxy(staRoxy_object)
#' staRoxy_object <- transform_oxy(staRoxy_object)
#' }
#'
#' # Get full descriptive statistics for all lipids and groups
#' desc_all <- descriptives_oxy(staRoxy_object)
#'
#' # Find the top 5 most abundant lipids, sorted by Coefficient of Variation (CV%)
#' desc_top <- descriptives_oxy(staRoxy_object, n_top = 5, sort_by = "cv")
#'
#' # Analyze specific groups and sort by the number of missing values (NAs)
#' desc_group <- descriptives_oxy(
#'   staRoxy_object,
#'   group = c("LPS", "Control"),
#'   sort_by = "na"
#' )
#'
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr left_join group_by summarise mutate across filter slice_max pull arrange desc all_of
#' @importFrom cli cli_h1 cli_alert_info cli_rule
#'
#' @export
descriptives_oxy <- function(obj, group = NULL, n_top = NULL, sort_by = "feature") {
  suppressWarnings({
    # Data Wrangling & Calculation
    # Convert matrix to long format, join metadata, and compute summary stats
    stats <- as.data.frame(obj$data) %>%
      tibble::rownames_to_column("Feature") %>%
      tidyr::pivot_longer(cols = -Feature, names_to = "sample", values_to = "Value") %>%
      dplyr::left_join(obj$meta, by = "sample") %>%
      dplyr::group_by(Feature, group) %>%
      dplyr::summarise(
        n_valid  = sum(!is.na(Value)),
        n_na     = sum(is.na(Value)),
        mean     = mean(Value, na.rm = TRUE),
        median   = stats::median(Value, na.rm = TRUE),
        sd       = stats::sd(Value, na.rm = TRUE),
        sem      = sd / sqrt(n_valid),
        cv_perc  = (sd / mean) * 100,
        .groups  = "drop"
      ) %>%
      dplyr::mutate(dplyr::across(tidyselect::where(is.numeric), ~ round(., 2)))

    # Filtering Logic
    # Filter by specific groups if provided
    if (!is.null(group)) {
      stats <- stats %>% dplyr::filter(group %in% !!group)
    }

    # Filter for top N features based on maximum mean value
    if (!is.null(n_top)) {
      top_features <- stats %>%
        dplyr::group_by(Feature) %>%
        dplyr::summarise(m = max(mean, na.rm = TRUE), .groups = "drop") %>%
        dplyr::slice_max(m, n = n_top, with_ties = FALSE) %>%
        dplyr::pull(Feature)

      stats <- stats %>% dplyr::filter(Feature %in% top_features)
    }

    # Sorting Logic
    # Reorder results based on the selected metric
    stats <- switch(sort_by,
                    "cv"      = stats %>% dplyr::arrange(dplyr::desc(cv_perc)),
                    "mean"    = stats %>% dplyr::arrange(dplyr::desc(mean)),
                    "median"  = stats %>% dplyr::arrange(dplyr::desc(median)),
                    "sd"      = stats %>% dplyr::arrange(dplyr::desc(sd)),
                    "sem"     = stats %>% dplyr::arrange(dplyr::desc(sem)),
                    "na"      = stats %>% dplyr::arrange(dplyr::desc(n_na)),
                    "n"       = stats %>% dplyr::arrange(dplyr::desc(n_valid)),
                    "feature" = stats %>% dplyr::arrange(Feature),
                    stats %>% dplyr::arrange(Feature)
    )
  })

  # Reporting & Console Output
  cli::cli_h1("Descriptive Statistics")

  available_groups <- paste(unique(stats$group), collapse = ', ')
  cli::cli_alert_info("Sorting by: {.val {sort_by}} | Groups: {.val {available_groups}}")

  print(as.data.frame(stats), row.names = FALSE)
  cli::cli_rule()

  return(invisible(stats))
}
