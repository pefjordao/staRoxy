#' @title Principal Component Analysis (PCA) Score Plot
#'
#' @description
#' Generates a publication-quality PCA score plot to visualize sample clustering,
#' group separation, and potential outliers. It includes automatic confidence
#' ellipses and variance explanation labels.
#'
#' @details
#' \strong{Data Pre-processing:}
#' The function utilizes the internal \code{prep_data_for_stats} pipeline,
#' which handles missing values (MNAR/MAR) and filters out group-exclusive
#' lipids. Features with near-zero variance (\eqn{SD < 1e-4}) are automatically
#' removed to ensure numerical stability during the \code{prcomp} execution.
#' Data is centered and scaled by default.
#'
#' \strong{Visualization:}
#' Confidence ellipses are drawn using \code{ggplot2::stat_ellipse} to
#' represent the distribution of each experimental group.
#'
#' @param obj A \code{staRoxy} object.
#' @param na_method Character. Imputation method: \code{"minprob"} (for MNAR)
#' or \code{"rf"} (Random Forest for MAR). Default is \code{"minprob"}.
#' @param na_threshold Numeric. Proportion threshold (0-1) for retaining features.
#' @param remove_exclusive Logical. If \code{TRUE} (default), excludes lipids
#' detected in only one group to prevent artificial clustering.
#' @param colors A named character vector for group colors. If \code{NULL},
#' a default palette is used.
#' @param title_size,x_axis_size,y_axis_size,legend_title_size,legend_size Numeric.
#' Font sizes for various plot elements.
#'
#' @return A \code{ggplot2} object.
#'
#' @importFrom stats sd prcomp
#' @importFrom ggplot2 ggplot aes geom_vline geom_hline stat_ellipse geom_point scale_color_manual scale_fill_manual labs theme element_text
#' @importFrom cowplot theme_half_open background_grid
#' @importFrom scales hue_pal
#'
#' @examples
#' # PCA with default colors
#' plot_pca_oxy(data_oxy_pellet)
#'
#' # PCA with custom group colors
#' my_colors <- c("Trypsin" = "#F65F72", "Scraper" = "#492376")
#' plot_pca_oxy(data_oxy_pellet, colors = my_colors, title_size = 16)
#'
#' @export
plot_pca_oxy <- function(obj,
                         na_method = "minprob",
                         na_threshold = 0.5,
                         remove_exclusive = TRUE,
                         colors = NULL,
                         title_size = 14,
                         x_axis_size = 12,
                         y_axis_size = 12,
                         legend_title_size = 12,
                         legend_size = 12) {

  # 1. Data Imputation and Pre-processing
  # Prepare the data matrix handling missing values and exclusive features
  m_imp <- prep_data_for_stats(
    obj,
    na_threshold = na_threshold,
    method = na_method,
    remove_exclusive = remove_exclusive
  )

  # Transpose for PCA (samples as rows, oxylipins as columns)
  pca_d <- t(m_imp)

  # 2. Feature Selection by Variation
  # Remove features with near-zero variance to avoid errors in prcomp
  sds <- apply(pca_d, 2, sd, na.rm = TRUE)
  pca_d <- pca_d[, sds > 1e-4, drop = FALSE]

  # 3. Principal Component Computation
  # Calculate PCA with centering and scaling
  pca_res <- prcomp(pca_d, center = TRUE, scale. = TRUE)

  # Calculate explained variance per component for axis labeling
  var_p <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

  # Prepare plotting dataframe with metadata
  df <- data.frame(pca_res$x, group = obj$meta$group)

  # Handle color palette if not provided
  if (is.null(colors)) {
    colors <- scales::hue_pal()(length(unique(df$group)))
  }

  # 4. Visualization (ggplot2)
  # Generate the PCA score plot with confidence ellipses
  ggplot2::ggplot(df, ggplot2::aes(PC1, PC2, color = group, fill = group)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    ggplot2::stat_ellipse(geom = "polygon", alpha = 0.15, show.legend = FALSE) +
    ggplot2::geom_point(size = 3, alpha = 0.8) +
    ggplot2::scale_color_manual(values = colors, name = "Group") +
    ggplot2::scale_fill_manual(values = colors, name = "Group") +
    ggplot2::labs(
      title = "Principal Component Analysis",
      x = paste0("PC1 (", var_p[1], "%)"),
      y = paste0("PC2 (", var_p[2], "%)")
    ) +
    cowplot::theme_half_open() +
    cowplot::background_grid() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(size = x_axis_size),
      axis.text = ggplot2::element_text(size = y_axis_size),
      legend.text = ggplot2::element_text(size = legend_size),
      legend.title = ggplot2::element_text(size = legend_title_size)
    )
}
