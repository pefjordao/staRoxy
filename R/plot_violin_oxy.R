#' @title Visualize Oxylipin Abundance with Violin and Boxplots
#'
#' @description
#' Generates a multi-layered visualization for a specific oxylipin, combining
#' violin plots (distribution), boxplots (summary statistics), and jittered points
#' (individual samples). It automatically handles scale transformations and
#' significance annotations from differential analysis results.
#'
#' @details
#' \strong{Scale Transformation:}
#' By default, data is plotted on the log2 scale. If \code{log_scale = FALSE},
#' values are back-transformed using \eqn{2^x - 1} for linear visualization
#' (e.g., in concentrations).
#'
#' \strong{Significance Bars:}
#' If the \code{stats} object (from \code{limma_oxy}) is provided, the function
#' parses the comparison names (e.g., "GroupA vs. GroupB") and adds
#' significance bars (\eqn{*}, \eqn{**}, \eqn{***}) for comparisons where the
#' adjusted P-value < 0.05.
#'
#' @param obj A \code{staRoxy} object.
#' @param stats List. The results list returned by \code{limma_oxy}.
#' If provided, significant comparisons are automatically annotated.
#' @param oxylipin Character. Name of the oxylipin feature to plot.
#' @param show_sig Logical. Whether to display significance bars. Default is \code{TRUE}.
#' @param colors A named character vector for group colors.
#' @param sig_y_nudge Numeric. Vertical adjustment for the significance bars.
#' @param log_scale Logical. If \code{TRUE} (default), uses the log2 scale.
#' If \code{FALSE}, transforms data to the linear scale.
#' @param unit Character. Measurement unit for the Y-axis label when \code{log_scale = FALSE}.
#' @param show_violin Logical. Whether to include the violin layer. Default is \code{TRUE}.
#' @param title_size,x_axis_size,y_axis_size Numeric. Font sizes for plot elements.
#' @param violin_width,boxplot_width Numeric. Aesthetic widths for the plot layers.
#' @param sig_bar_size,sig_text_size Numeric. Sizes for the significance bars and asterisks.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_jitter scale_fill_manual labs theme element_text
#' @importFrom cowplot theme_half_open background_grid
#' @importFrom ggsignif geom_signif
#' @importFrom scales hue_pal
#' @importFrom dplyr case_when
#' @importFrom cli cli_alert_danger cli_alert_warning
#'
#' @examples
#' # Plot 9,10-DiHOME abundance with significance bars
#' plot_violin_oxy(data_oxy_pellet,
#'                 stats = results_limma,
#'                 oxylipin = "9,10-DiHOME")
#'
#' # Plot in linear scale with custom units
#' plot_violin_oxy(data_oxy_pellet,
#'                 oxylipin = "PGE2",
#'                 log_scale = FALSE,
#'                 unit = "nmol/mg")
#'
#' @export
plot_violin_oxy <- function(obj, stats = NULL, oxylipin, show_sig = TRUE, colors = NULL,
                            sig_y_nudge = 0, log_scale = TRUE, unit = "unit",
                            title_size = 14, x_axis_size = 12, y_axis_size = 12,
                            show_violin = TRUE, violin_width = 0.8, boxplot_width = 0.15,
                            sig_bar_size = 0.6, sig_text_size = 6) {

  # 1. Data Validation and Extraction
  # Check if the requested oxylipin exists in the dataset
  if (!(oxylipin %in% rownames(obj$data))) {
    return(cli::cli_alert_danger("Oxylipin '{oxylipin}' not found."))
  }

  # Extract values and handle scale transformation (Log2 to Linear)
  vals <- as.numeric(obj$data[oxylipin, ])
  if (!log_scale) vals <- 2^vals - 1

  plot_df <- data.frame(Value = vals, Group = obj$meta$group)

  # Remove non-finite values (NAs, Inf) that might break the plot
  is_valid <- is.finite(plot_df$Value)
  na_count <- sum(!is_valid)
  if (na_count > 0) {
    cli::cli_alert_warning("{na_count} missing/infinite value(s) removed from '{oxylipin}' plot.")
    plot_df <- plot_df[is_valid, ]
  }

  # Setup color palette
  if (is.null(colors)) {
    colors <- scales::hue_pal()(length(levels(plot_df$Group)))
  }

  # 2. Base Plot Construction
  # Initialize ggplot with violin and boxplot layers
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Group, y = Value, fill = Group))

  if (show_violin) {
    p <- p + ggplot2::geom_violin(
      trim = FALSE, alpha = 0.4, linewidth = 0.3, width = violin_width
    )
  }

  p <- p +
    ggplot2::geom_boxplot(width = boxplot_width, outlier.shape = NA, alpha = 0.5) +
    ggplot2::geom_jitter(shape = 21, width = 0.1, size = 2, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = oxylipin,
      y = if (log_scale) "Abundance (log2)" else paste0("Abundance (", unit, ")"),
      x = ""
    ) +
    cowplot::theme_half_open() +
    cowplot::background_grid() +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_size),
      axis.title = ggplot2::element_text(size = x_axis_size),
      axis.text = ggplot2::element_text(size = y_axis_size)
    )

  # 3. Significance Annotation
  # If stats (limma) are provided, calculate p-value positions and labels
  if (show_sig && !is.null(stats)) {
    rng <- diff(range(plot_df$Value, na.rm = TRUE))
    sig_df <- data.frame(
      group1 = character(),
      group2 = character(),
      label = character(),
      y.position = numeric()
    )

    for (comp in names(stats)) {
      # Split comparison string (e.g., "GroupA vs. GroupB") to identify groups
      grps <- unlist(strsplit(comp, " vs\\. "))

      if (length(grps) < 2) next

      if (oxylipin %in% rownames(stats[[comp]])) {
        p_v <- stats[[comp]][oxylipin, "adj.P.Val"]
        lab <- dplyr::case_when(
          p_v < 0.001 ~ "***",
          p_v < 0.01  ~ "**",
          p_v < 0.05  ~ "*",
          TRUE        ~ "ns"
        )

        # Only add to plot if significant
        if (lab != "ns") {
          y_base <- max(plot_df$Value, na.rm = TRUE)
          y_pos <- y_base + (rng * (0.1 + (nrow(sig_df) * 0.15) + sig_y_nudge))

          sig_df <- rbind(sig_df, data.frame(
            group1 = grps[1],
            group2 = grps[2],
            label = lab,
            y.position = y_pos
          ))
        }
      }
    }

    # Apply annotations using ggsignif
    if (nrow(sig_df) > 0) {
      comps_list <- lapply(1:nrow(sig_df), function(i) {
        c(as.character(sig_df$group1[i]), as.character(sig_df$group2[i]))
      })

      p <- p + ggsignif::geom_signif(
        comparisons = comps_list,
        annotations = sig_df$label,
        y_position = sig_df$y.position,
        tip_length = 0.02,
        size = sig_bar_size,
        textsize = sig_text_size,
        vjust = 0.5
      )
    }
  }

  return(p)
}
