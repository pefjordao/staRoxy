#' @title Visualize Oxylipin Abundance Profiles
#'
#' @description
#' Generates a multi-layered visualization combining violin plots, boxplots,
#' and jittered points. It dynamically adapts its layout: plotting experimental
#' groups on the X-axis for a single oxylipin, or plotting multiple oxylipins
#' on the X-axis dodged by group for broader profiling.
#'
#' @details
#' \strong{Scale Transformation:}
#' By default, data is plotted on the log2 scale. If \code{log_scale = FALSE},
#' values are back-transformed using \eqn{2^x - 1} for linear visualization.
#'
#' \strong{Dynamic Layout & Significance Bars:}
#' \itemize{
#'   \item \strong{Single Oxylipin:} If \code{stats} from \code{limma_oxy} are
#'   provided, significance bars (\eqn{*}, \eqn{**}, \eqn{***}) are automatically
#'   added for comparisons where the adjusted P-value < 0.05.
#'   \item \strong{Multiple Oxylipins:} The X-axis switches to display oxylipin
#'   names, and groups are dodged side-by-side. Significance bars are disabled
#'   in this mode to prevent visual clutter.
#' }
#'
#' @param obj A \code{staRoxy} object.
#' @param oxylipins Character vector. Name(s) of the oxylipin(s) to plot.
#' @param stats List. The results list returned by \code{limma_oxy}.
#' Only applied when plotting a single oxylipin.
#' @param target_group Character. Optional filter for a specific experimental group
#' (useful when plotting multiple oxylipins).
#' @param show_sig Logical. Whether to display significance bars (single lipid only).
#' @param colors A named character vector for group colors.
#' @param sig_y_nudge Numeric. Vertical adjustment for the significance bars.
#' @param log_scale Logical. If \code{TRUE} (default), uses the log2 scale.
#' @param unit Character. Y-axis unit label when \code{log_scale = FALSE}.
#' @param show_violin Logical. Whether to include the violin layer. Default is \code{TRUE}.
#' @param title_size,x_axis_size,y_axis_size,legend_title_size,legend_size Numeric. Font sizes.
#' @param violin_width,boxplot_width,point_size,sig_bar_size,sig_text_size Numeric. Sizes for plot layers.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' \dontshow{
#' staRoxy_object <- read_oxy(data_oxy_lps_pellet, metadata_oxy_lps_pellet)
#' staRoxy_object <- filter_oxy(staRoxy_object)
#' staRoxy_object <- transform_oxy(staRoxy_object)
#' diff_res <- limma_oxy(staRoxy_object, save_results = FALSE)
#' }
#'
#'\dontrun{
#' # Single Oxylipin with Significance Bars (Log2 Scale)
#' plot_violin_oxy(
#'   obj = staRoxy_object,
#'   stats = diff_res,
#'   oxylipins = "PGE2"
#' )
#'
#' # Multiple Oxylipins Profiling (Significance bars auto-disabled)
#' target_lipids <- c("PGE2", "PGD2", "11-HETE", "15-HETrE")
#' plot_violin_oxy(
#'   obj = staRoxy_object,
#'   oxylipins = target_lipids
#' )
#'
#' # Multiple Oxylipins isolated to a single target group (Linear Scale)
#' plot_violin_oxy(
#'   obj = staRoxy_object,
#'   oxylipins = target_lipids,
#'   target_group = "LPS",
#'   log_scale = FALSE,
#'   unit = "fg/mL"
#' )
#'}
#'
#' @importFrom ggplot2 ggplot aes geom_violin geom_boxplot geom_point scale_fill_manual labs theme element_text position_dodge position_jitterdodge position_jitter
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter all_of case_when
#' @importFrom cowplot theme_half_open background_grid
#' @importFrom ggsignif geom_signif
#' @importFrom scales hue_pal
#' @importFrom cli cli_alert_danger cli_alert_warning cli_alert_info
#'
#' @export
plot_violin_oxy <- function(obj,
                            oxylipins,
                            stats = NULL,
                            target_group = NULL,
                            show_sig = TRUE,
                            colors = NULL,
                            sig_y_nudge = 0.6,
                            log_scale = TRUE,
                            unit = "unit",
                            title_size = 14,
                            x_axis_size = 12,
                            y_axis_size = 12,
                            legend_title_size = 12,
                            legend_size = 12,
                            show_violin = TRUE,
                            violin_width = 0.8,
                            boxplot_width = 0.2,
                            point_size = 1.5,
                            sig_bar_size = 0.6,
                            sig_text_size = 6) {

  # Validation and Extraction
  valid_oxys <- oxylipins[oxylipins %in% rownames(obj$data)]
  if (length(valid_oxys) == 0) {
    cli::cli_alert_danger("None of the requested oxylipins were found.")
    return(invisible(NULL))
  }

  is_single <- length(valid_oxys) == 1

  # Extract values and handle scale transformation
  ext_data <- as.data.frame(t(obj$data[valid_oxys, , drop = FALSE]))
  if (!log_scale) ext_data <- 2^ext_data - 1
  ext_data$Group <- obj$meta$group

  plot_df <- tidyr::pivot_longer(ext_data, cols = dplyr::all_of(valid_oxys),
                                 names_to = "Oxylipin", values_to = "Value")

  # Remove non-finite values (NAs, Inf)
  is_valid <- is.finite(plot_df$Value)
  na_count <- sum(!is_valid)
  if (na_count > 0) {
    cli::cli_alert_warning("{na_count} missing/infinite value(s) removed.")
    plot_df <- plot_df[is_valid, ]
  }

  # Target Group Filter
  if (!is.null(target_group)) {
    plot_df <- dplyr::filter(plot_df, Group == target_group)
  }

  # Setup color palette
  if (is.null(colors)) {
    colors <- scales::hue_pal()(length(unique(plot_df$Group)))
  }

  # Dynamic Mapping Logic
  if (is_single && is.null(target_group)) {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Group, y = Value, fill = Group))
    dodge_w <- 0
    p_title <- valid_oxys[1]
    x_angle <- 0
    x_hjust <- 0.5
    show_leg <- FALSE
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Oxylipin, y = Value, fill = Group))
    dodge_w <- if (is.null(target_group)) 0.8 else 0
    p_title <- if (is.null(target_group)) "Oxylipin Abundance Profile" else target_group
    x_angle <- 90
    x_hjust <- 1
    show_leg <- is.null(target_group)
  }

  # Position for Jitter
  if (dodge_w > 0) {
    pos_jitter <- ggplot2::position_jitterdodge(jitter.width = 0.1, dodge.width = dodge_w)
  } else {
    pos_jitter <- ggplot2::position_jitter(width = 0.1)
  }

  # Base Plot Construction
  if (show_violin) {
    p <- p + ggplot2::geom_violin(
      position = ggplot2::position_dodge(width = dodge_w),
      trim = FALSE, alpha = 0.4, linewidth = 0.3, width = violin_width,
      show.legend = show_leg
    )
  }

  p <- p +
    ggplot2::geom_boxplot(
      position = ggplot2::position_dodge(width = dodge_w),
      width = boxplot_width, outlier.shape = NA, alpha = 0.5,
      color = "black", linewidth = 0.5, key_glyph = "rect", show.legend = FALSE
    ) +
    ggplot2::geom_point(
      position = pos_jitter,
      shape = 21, size = point_size, alpha = 0.8, show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = p_title,
      y = if (log_scale) "Abundance (log2)" else paste0("Abundance (", unit, ")"),
      x = "", fill = "Group"
    ) +
    cowplot::theme_half_open() +
    cowplot::background_grid(major = "y") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_size),
      axis.title = ggplot2::element_text(size = x_axis_size),
      axis.text.y = ggplot2::element_text(size = y_axis_size),
      axis.text.x = ggplot2::element_text(size = x_axis_size, angle = x_angle, vjust = 0.5, hjust = x_hjust),
      legend.position = if (show_leg) "right" else "none",
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.text = ggplot2::element_text(size = legend_size)
    )

  # Significance Annotation (Single Oxylipin Only)
  if (show_sig && !is.null(stats)) {
    if (is_single && is.null(target_group)) {
      rng <- diff(range(plot_df$Value, na.rm = TRUE))
      sig_df <- data.frame(group1 = character(), group2 = character(), label = character(), y.position = numeric())
      oxylipin <- valid_oxys[1]

      for (comp in names(stats)) {
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

          if (lab != "ns") {
            y_base <- max(plot_df$Value, na.rm = TRUE)
            y_pos <- y_base + (rng * (0.1 + (nrow(sig_df) * 0.15) + sig_y_nudge))
            sig_df <- rbind(sig_df, data.frame(group1 = grps[1], group2 = grps[2], label = lab, y.position = y_pos))
          }
        }
      }

      if (nrow(sig_df) > 0) {
        comps_list <- lapply(1:nrow(sig_df), function(i) c(as.character(sig_df$group1[i]), as.character(sig_df$group2[i])))
        p <- p + ggsignif::geom_signif(
          comparisons = comps_list, annotations = sig_df$label, y_position = sig_df$y.position,
          tip_length = 0.02, size = sig_bar_size, textsize = sig_text_size, vjust = 0.5
        )
      }
    } else {
      cli::cli_alert_info("Significance bars are disabled when plotting multiple oxylipins simultaneously.")
    }
  }

  return(p)
}
