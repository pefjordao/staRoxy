#' @title Over-representation Analysis (ORA) of Fatty Acid Precursors
#'
#' @description
#' Performs an Over-representation Analysis (ORA) to determine if the diversity
#' of oxylipins detected in experimental groups deviates significantly from
#' a theoretical reference universe.
#'
#' @details
#' The function uses a hypergeometric test (one-vs-all) to compare the observed
#' diversity of each precursor class against its frequency in a reference universe.
#'
#' \strong{Enrichment Score (ES):}
#' Calculated as the ratio between the observed proportion in the group and
#' the expected proportion in the universe:
#' \deqn{ES = \frac{k/n}{K/N}}
#' Where \eqn{k} is the number of species from a precursor class detected,
#' \eqn{n} is the total detected in the group, \eqn{K} is the total in the
#' universe, and \eqn{N} is the universe size.
#'
#' \strong{Statistical Significance:}
#' Computed via the upper tail of the hypergeometric distribution:
#' \deqn{P(X \ge k) = \sum_{i=k}^{min(n,K)} \frac{\binom{K}{i}\binom{N-K}{n-i}}{\binom{N}{n}}}
#' An \eqn{ES > 1} indicates over-representation (enrichment), while
#' \eqn{ES < 1} indicates under-representation (depletion).
#'
#' @param obj A \code{staRoxy} object.
#' @param universe A data frame mapping oxylipins to precursors. Defaults to
#' the output of \code{get_oxy_universe()}.
#' @param compare_two_groups Logical. If \code{TRUE}, filters analysis for
#' specific group comparisons. Default is \code{FALSE}.
#' @param group1,group2 Character. Groups to compare if \code{compare_two_groups = TRUE}.
#' @param score_plot Logical. If \code{TRUE}, generates an Enrichment Score bar plot.
#' @param palette Character. ColorBrewer palette for the plot (default "Spectral").
#' @param title_size,x_axis_size,y_axis_size,legend_size,legend_title_size,sig_size Numeric. Plot aesthetics and font sizes.
#'
#' @return Invisibly returns a data frame with ORA statistics, including
#' P-values for enrichment and depletion, and Enrichment Scores.
#'
#' @examples
#' \dontshow{
#' staRoxy_object <- read_oxy(data_oxy_lps_pellet, metadata_oxy_lps_pellet)
#' staRoxy_object <- filter_oxy(staRoxy_object)
#' staRoxy_object <- transform_oxy(staRoxy_object)
#' }
#'
#'\dontrun{
#' # Run ORA using the default universe
#' oxy_ora(staRoxy_object)
#'
#' # Run ORA using a custom universe
#' custom_univ <- edit_oxy_universe(remove = c("PGE2", "PGD2"))
#' oxy_ora(staRoxy_object, universe = custom_univ)
#'}
#'
#' @importFrom dplyr distinct count mutate left_join filter case_when
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom stats phyper
#' @importFrom ggplot2 ggplot aes geom_bar geom_text geom_hline scale_fill_brewer labs theme element_text position_dodge
#' @importFrom cowplot theme_half_open
#' @importFrom cli cli_h1 cli_alert_success cli_bullets cli_alert_info cli_alert_warning cli_rule
#'
#' @export
oxy_ora <- function(obj,
                    universe = get_oxy_universe(),
                    compare_two_groups = FALSE,
                    group1 = NULL,
                    group2 = NULL,
                    score_plot = TRUE,
                    palette = "Spectral",
                    title_size = 14,
                    x_axis_size = 12,
                    y_axis_size = 12,
                    legend_title_size = 12,
                    legend_size = 12,
                    sig_size = 6) {

  # Reference Universe Mapping
  univ_map <- universe
  univ_size <- nrow(univ_map)
  univ_counts_df <- univ_map %>% dplyr::count(Precursor, name = "Universe_M")

  # Data Wrangling
  # Pivot and join data with biochemical mapping; check.names=FALSE preserves IDs
  suppressMessages({
    all_data <- as.data.frame(t(as.matrix(obj$data)), check.names = FALSE) %>%
      tibble::rownames_to_column("sample") %>%
      tidyr::pivot_longer(-sample, names_to = "Oxylipin", values_to = "Abund") %>%
      dplyr::mutate(Oxylipin = trimws(Oxylipin)) %>%
      dplyr::left_join(obj$meta, by = "sample") %>%
      dplyr::left_join(univ_map, by = "Oxylipin") %>%
      dplyr::filter(!is.na(Precursor))
  })

  if (compare_two_groups) {
    all_data <- all_data %>% dplyr::filter(group %in% c(group1, group2))
    all_data$group <- factor(all_data$group, levels = c(group1, group2))
  }

  # Identify detected Oxylipins (is.finite handles log-scale data/negatives)
  detected_Oxylipins <- all_data %>%
    dplyr::filter(is.finite(Abund)) %>%
    dplyr::distinct(group, Oxylipin, Precursor)

  summary <- detected_Oxylipins %>%
    dplyr::count(group, Precursor, name = "Observed")

  group_totals <- detected_Oxylipins %>%
    dplyr::count(group, name = "Total_Group_Detected")

  # Enrichment Statistics
  # Hypergeometric test for Over-representation and Depletion
  summary <- summary %>%
    dplyr::left_join(group_totals, by = "group") %>%
    dplyr::left_join(univ_counts_df, by = "Precursor") %>%
    dplyr::mutate(
      Universe_N = univ_size - Universe_M,
      p_enr = stats::phyper(Observed - 1, Universe_M, Universe_N, Total_Group_Detected, lower.tail = FALSE),
      p_dep = stats::phyper(Observed, Universe_M, Universe_N, Total_Group_Detected, lower.tail = TRUE),
      p_val = pmin(p_enr, p_dep),
      Enrich_Score = (Observed / Total_Group_Detected) / (Universe_M / univ_size),
      sig_label = dplyr::case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        TRUE          ~ ""
      )
    )

  # CLI Reporting
  cli::cli_h1("Detailed Over-representation analysis (ORA)")

  # Report Significant Enrichment
  enr <- summary %>% dplyr::filter(p_enr < 0.05)
  if (nrow(enr) > 0) {
    cli::cli_alert_success("Significant ENRICHMENT found:")
    for (i in 1:nrow(enr)) {
      cli::cli_bullets(c("+" = "{enr$Precursor[i]} ({enr$group[i]}): ES = {round(enr$Enrich_Score[i], 2)}x | P = {format.pval(enr$p_enr[i], digits = 3)}"))
    }
  } else {
    cli::cli_alert_info("No significant precursor ENRICHMENT detected.")
  }

  # Report Significant Depletion
  dep <- summary %>% dplyr::filter(p_dep < 0.05)
  if (nrow(dep) > 0) {
    cli::cli_alert_warning("Significant DEPLETION found:")
    for (i in 1:nrow(dep)) {
      cli::cli_bullets(c("-" = "{dep$Precursor[i]} ({dep$group[i]}): ES = {round(dep$Enrich_Score[i], 2)}x | P = {format.pval(dep$p_dep[i], digits = 3)}"))
    }
  }
  cli::cli_rule()

  # Visualization
  if (score_plot) {
    p2 <- ggplot2::ggplot(summary, ggplot2::aes(x = group, y = Enrich_Score, fill = Precursor)) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(width = 0.9), color = "black", linewidth = 0.4) +
      ggplot2::geom_text(ggplot2::aes(label = sig_label, y = Enrich_Score + 0.05),
                         position = ggplot2::position_dodge(width = 0.9), vjust = 0, size = sig_size) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
      ggplot2::scale_fill_brewer(palette = palette) +
      ggplot2::labs(
        title = "Precursor Enrichment",
        subtitle = "ES > 1 (Enriched) or < 1 (Depleted) relative to the reference universe",
        x = "", y = "Enrichment Score (ES)"
      ) +
      cowplot::theme_half_open() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5),
        plot.subtitle = ggplot2::element_text(size = title_size - 2, hjust = 0.5),
        axis.text = ggplot2::element_text(size = x_axis_size),
        axis.title = ggplot2::element_text(size = y_axis_size),
        legend.text = ggplot2::element_text(size = legend_size),
        legend.title = ggplot2::element_text(size = legend_title_size)
      )

    suppressWarnings(print(p2))
  }

  return(invisible(summary))
}
