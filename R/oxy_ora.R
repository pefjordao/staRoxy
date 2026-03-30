#' @title Over-representation Analysis (ORA) of Fatty Acid Precursors
#'
#' @description
#' Performs an Over-representation Analysis (ORA) to determine if the diversity
#' of oxylipins detected in experimental groups deviates significantly from
#' a theoretical reference universe.
#'
#' @details
#' The function uses a hypergeometric test (one-vs-all) to compare the observed
#' diversity of each precursor class against its frequency in a reference
#' universe of 125 oxylipins.
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
#' @param compare_two_groups Logical. If \code{TRUE}, filters analysis for
#' specific group comparisons. Default is \code{FALSE}.
#' @param group1,group2 Character. Groups to compare if \code{compare_two_groups = TRUE}.
#' @param score_plot Logical. If \code{TRUE}, generates an Enrichment Score bar plot.
#' @param palette Character. ColorBrewer palette for the plot (default "Spectral").
#' @param title_size,x_axis_size,y_axis_size,legend_size,legend_title_size,sig_size Numeric.
#' Plot aesthetics and font sizes.
#'
#' @return Invisibly returns a data frame with ORA statistics, including
#' P-values for enrichment and depletion, and Enrichment Scores.
#'
#' @importFrom dplyr distinct count mutate left_join filter case_when
#' @importFrom tidyr pivot_longer
#' @importFrom tibble rownames_to_column
#' @importFrom stats phyper
#' @importFrom ggplot2 ggplot aes geom_bar geom_text geom_hline scale_fill_brewer labs theme element_text position_dodge
#' @importFrom cowplot theme_half_open
#' @importFrom cli cli_h1 cli_alert_success cli_bullets cli_alert_info cli_alert_warning cli_rule
#'
#' @examples
#' # Perform ORA on the default pellet dataset
#' ora_results <- oxy_ora(data_oxy_pellet)
#'
#' # Compare ORA specifically between two groups
#' oxy_ora(data_oxy_pellet, compare_two_groups = TRUE,
#'         group1 = "Scraper", group2 = "Trypsin")
#'
#' @export
oxy_ora <- function(obj,
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

  # 1. Reference Universe Mapping
  # Define the biochemical relationship between oxylipins and fatty acid precursors
  univ_map <- data.frame(
    Lipid = c("tetranor-12-HETE", "12-HHTrE", "13-oxo-ODE", "9-HOTrE", "9-oxo-ODE", "13-HOTrE-γ", "13-HOTrE", "9-HODE", "10-HODE", "12-HODE", "13-HODE", "9(10)-EpOME", "12(13)-EpOME", "9-HOA", "10-HOA", "9,10-DiHOME", "12,13-DiHOME", "15-deoxy-Δ12,14-PGJ2", "15-oxo-ETE", "5-HEPE", "9-HEPE", "8-HEPE", "11-HEPE", "12-HEPE", "15-HEPE", "18-HEPE", "14(15)-EpETE", "17(18)-EpETE", "5-HETE", "8-HETE", "9-HETE", "11-HETE", "12-HETE", "14 / 15-HETE", "16-HETE", "17-HETE", "18-HETE", "19-HETE", "20-HETE", "5(6)-EpETrE", "8(9)-EpETrE", "11(12)-EpETrE", "14(15)-EpETrE", "15-oxo-EDE", "5-HETrE", "15-HETrE", "2,3-dinor-8-iso-PGF2α", "10-Nitrooleate", "9-Nitrooleate", "tetranor-PGDM", "PGA2", "PGB2", "PGJ2", "15-deoxy-Δ12,14-PGD2", "12-oxo-LTB4", "LTB4", "6-trans-LTB4", "6-trans-12-epi-LTB4", "12-epi-LTB4", "8,15-DiHETE", "14,15-DiHETE", "8,9-DiHETrE", "17-oxo-DHA", "2,3-dinor-TXB2", "4-HDoHE", "7-HDoHE", "8-HDoHE", "10-HDoHE", "11-HDoHE", "13-HDoHE", "14-HDoHE", "16 / 17-HDoHE", "19-HDoHE", "20-HDoHE", "7(8)-EpDPE", "10(11)-EpDPE", "13(14)-EpDPE", "16(17)-EpDPE", "19(20)-EpDPE", "17-oxo-DPA", "PGE3", "15-keto-PGE2", "LXA5", "RvE1", "epi-LXA4", "15-epi-LXA4", "LXB4", "20-hydroxy-LTB4", "PGD2", "13,14-dihydro-15-keto-PGD2", "PGE2", "13,14-dihydro-15-keto-PGE2", "15-keto-PGF2α", "8-iso-15-keto-PGF2β", "PGF3α", "8-iso-PGF3α", "15-F2t-Isoprostane", "11β-13,14-dihydro-15-keto-PGF2α", "11β-PGF2α", "15-keto-PGF1α", "8-iso-PGE1", "PGD1", "PGF2α", "8,12-iso-iPF2α-VI", "13,14-dihydro-PGF2α", "PGF1α", "PDX", "Maresin 1", "Maresin 2", "RvD5", "19,20-DiHDPA", "20-carboxy-LTB4", "TXB3", "20-hydroxy-PGE2", "6-keto-PGE1", "Δ17-6-keto-PGF1α", "TXB2", "6-keto-PGF1α", "6,15-diketo-13,14-dihydro-PGF1α", "TXB1", "RvD1", "RvD2", "RvD3", "16,16-dimethyl-PGD2", "1a,1b-dihomo-PGE2"),
    Precursor = c("ARA (20:4 n-6)", "ARA (20:4 n-6)", "LA (18:2 n-6)", "ALA (18:3 n-3)", "LA (18:2 n-6)", "GLA (18:3 n-6)", "ALA (18:3 n-3)", "LA (18:2 n-6)", "LA (18:2 n-6)", "LA (18:2 n-6)", "LA (18:2 n-6)", "LA (18:2 n-6)", "LA (18:2 n-6)", "OA (18:1 n-9)", "OA (18:1 n-9)", "LA (18:2 n-6)", "LA (18:2 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "EDA (20:2 n-6)", "MA (20:3 n-9)", "DGLA (20:3 n-6)", "ARA (20:4 n-6)", "OA (18:1 n-9)", "OA (18:1 n-9)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "DHA (22:6 n-3)", "ARA (20:4 n-6)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DPA (22:5 n-3)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "EPA (20:5 n-3)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "ARA (20:4 n-6)", "EPA (20:5 n-3)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "ARA (20:4 n-6)", "DGLA (20:3 n-6)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "DHA (22:6 n-3)", "ARA (20:4 n-6)", "ARA (20:4 n-6)")
  ) %>% dplyr::distinct(Lipid, .keep_all = TRUE)

  univ_size <- nrow(univ_map)
  univ_counts_df <- univ_map %>% dplyr::count(Precursor, name = "Universe_M")

  # 2. Data Wrangling
  # Pivot and join data with biochemical mapping; check.names=FALSE preserves IDs
  suppressMessages({
    all_data <- as.data.frame(t(obj$data), check.names = FALSE) %>%
      tibble::rownames_to_column("sample") %>%
      tidyr::pivot_longer(-sample, names_to = "Lipid", values_to = "Abund") %>%
      dplyr::mutate(Lipid = trimws(Lipid)) %>%
      dplyr::left_join(obj$meta, by = "sample") %>%
      dplyr::left_join(univ_map, by = "Lipid") %>%
      dplyr::filter(!is.na(Precursor))
  })

  if (compare_two_groups) {
    all_data <- all_data %>% dplyr::filter(group %in% c(group1, group2))
    all_data$group <- factor(all_data$group, levels = c(group1, group2))
  }

  # Identify detected lipids (is.finite handles log-scale data/negatives)
  detected_lipids <- all_data %>%
    dplyr::filter(is.finite(Abund)) %>%
    dplyr::distinct(group, Lipid, Precursor)

  resumo <- detected_lipids %>%
    dplyr::count(group, Precursor, name = "Observed")

  group_totals <- detected_lipids %>%
    dplyr::count(group, name = "Total_Group_Detected")

  # 3. Enrichment Statistics
  # Hypergeometric test for Over-representation and Depletion
  resumo <- resumo %>%
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

  # 4. CLI Reporting
  cli::cli_h1("Detailed Over-representation analysis (ORA)")

  # Report Significant Enrichment
  enr <- resumo %>% dplyr::filter(p_enr < 0.05)
  if (nrow(enr) > 0) {
    cli::cli_alert_success("Significant ENRICHMENT found:")
    for (i in 1:nrow(enr)) {
      cli::cli_bullets(c("+" = "{enr$Precursor[i]} ({enr$group[i]}): ES = {round(enr$Enrich_Score[i], 2)}x | P = {format.pval(enr$p_enr[i], digits = 3)}"))
    }
  } else {
    cli::cli_alert_info("No significant precursor ENRICHMENT detected.")
  }

  # Report Significant Depletion
  dep <- resumo %>% dplyr::filter(p_dep < 0.05)
  if (nrow(dep) > 0) {
    cli::cli_alert_warning("Significant DEPLETION found:")
    for (i in 1:nrow(dep)) {
      cli::cli_bullets(c("-" = "{dep$Precursor[i]} ({dep$group[i]}): ES = {round(dep$Enrich_Score[i], 2)}x | P = {format.pval(dep$p_dep[i], digits = 3)}"))
    }
  }
  cli::cli_rule()

  # 5. Visualization
  if (score_plot) {
    p2 <- ggplot2::ggplot(resumo, ggplot2::aes(x = group, y = Enrich_Score, fill = Precursor)) +
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

  return(invisible(resumo))
}
