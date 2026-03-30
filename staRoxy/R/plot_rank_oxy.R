#' @title Plot Ranked Differential Abundance
#'
#' @description
#' Visualizes the top differentially abundant oxylipins using a horizontal
#' bar plot. Features are ranked by either Fold-Change or Significance.
#'
#' @details
#' Standard Error of the Mean (SEM) is estimated as \code{abs(logFC / t)}.
#' Significance levels are indicated by asterisks based on adjusted P-values.
#'
#' @param stats List of dataframes from \code{limma_oxy} or a single dataframe.
#' @param contrast Numeric (index) or Character (name) of the contrast to plot.
#' @param top_n Integer. Number of oxylipins to display.
#' @param rank_by Character. Rank by \code{"logFC"} (absolute) or \code{"p"} (significance).
#' @param palette Character. Viridis palette option (e.g., "magma", "viridis").
#' @param direction Integer. Direction of the color scale (1 or -1).
#' @param title_size,x_axis_size,y_axis_size,legend_title_size,legend_size Font sizes.
#' @param sig_size Numeric. Size of the significance asterisks.
#'
#' @return A ggplot2 object.
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr mutate slice_max case_when
#' @importFrom forcats fct_reorder
#' @importFrom ggplot2 ggplot aes geom_errorbar geom_bar geom_text scale_fill_viridis_c labs theme element_text
#' @importFrom cowplot theme_half_open background_grid
#' @importFrom cli cli_alert_danger
#'
#' @export
plot_rank_oxy <- function(stats,
                          contrast = 1,
                          top_n = 10,
                          rank_by = "logFC",
                          palette = "magma",
                          direction = -1,
                          title_size = 14,
                          x_axis_size = 12,
                          y_axis_size = 12,
                          legend_title_size = 12,
                          legend_size = 12,
                          sig_size = 6) {

  # 1. Contrast Selection and Data Extraction
  # Handle both list of dataframes (from limma_oxy) or a single dataframe
  if (is.list(stats) && !is.data.frame(stats)) {
    if (is.character(contrast)) {
      if (!(contrast %in% names(stats))) {
        return(cli::cli_alert_danger("Contrast '{contrast}' not found in the stats list."))
      }
      comp_name <- contrast
      stats_df <- stats[[contrast]]
    } else if (is.numeric(contrast)) {
      if (contrast > length(stats)) {
        return(cli::cli_alert_danger("Index {contrast} exceeds available contrasts."))
      }
      comp_name <- names(stats)[contrast]
      stats_df <- stats[[contrast]]
    }
  } else {
    comp_name <- attr(stats, "comp")
    stats_df <- stats
  }

  # 2. Title and Scoring Logic
  full_title <- if (!is.null(comp_name)) {
    paste("Oxylipin Differential Abundance -", gsub("-", " vs. ", comp_name))
  } else {
    "Oxylipin Differential Abundance"
  }

  # Determine ranking criteria (Absolute LogFC or Significance)
  select_score <- if (rank_by == "p") (1 - stats_df$adj.P.Val) else abs(stats_df$logFC)

  # 3. Data Wrangling (dplyr & forcats)
  plot_data <- stats_df %>%
    tibble::rownames_to_column("Oxylipin") %>%
    dplyr::mutate(score = select_score) %>%
    dplyr::slice_max(order_by = score, n = top_n, with_ties = FALSE) %>%
    dplyr::mutate(
      SEM = abs(logFC / t),
      # Define plot order based on user preference
      plot_order = if (rank_by == "p") (1 - adj.P.Val) else logFC,
      Oxylipin = forcats::fct_reorder(Oxylipin, plot_order),
      # Significant labels (asterisks)
      sig_label = dplyr::case_when(
        adj.P.Val < 0.001 ~ "***",
        adj.P.Val < 0.01  ~ "**",
        adj.P.Val < 0.05  ~ "*",
        TRUE              ~ ""
      ),
      # Calculate error bar boundaries
      err_min = ifelse(logFC > 0, logFC, logFC - SEM),
      err_max = ifelse(logFC > 0, logFC + SEM, logFC)
    )

  # 4. Visualization (ggplot2)
  ggplot2::ggplot(plot_data, ggplot2::aes(x = logFC, y = Oxylipin, fill = adj.P.Val)) +
    # Error bars for precision
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = err_min, xmax = err_max),
      width = 0.4,
      linewidth = 0.5
    ) +
    # Rank bars
    ggplot2::geom_bar(stat = "identity", color = "black", linewidth = 0.6) +
    # Significance labels positioned dynamically
    ggplot2::geom_text(
      ggplot2::aes(x = ifelse(logFC > 0, err_max, err_min), label = sig_label),
      hjust = ifelse(plot_data$logFC > 0, -0.4, 1.4),
      vjust = 0.7,
      size = sig_size
    ) +
    # Styling and scales
    ggplot2::scale_fill_viridis_c(
      option = palette,
      direction = direction,
      name = "Adj. P-value"
    ) +
    ggplot2::labs(title = full_title, x = "Log2FC ± SEM", y = "") +
    cowplot::theme_half_open() +
    cowplot::background_grid(major = "x") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5),
      axis.title = ggplot2::element_text(size = x_axis_size),
      axis.text = ggplot2::element_text(size = y_axis_size),
      legend.text = ggplot2::element_text(size = legend_size),
      legend.title = ggplot2::element_text(size = legend_title_size)
    )
}
