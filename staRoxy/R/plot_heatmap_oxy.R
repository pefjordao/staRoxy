#' @title Z-score Heatmap for Oxylipidomics Profiles
#'
#' @description
#' Generates a publication-ready hierarchical clustering heatmap. Abundances are
#' standardized via Z-score transformation to allow comparison across oxylipins
#' with different magnitudes.
#'
#' @details
#' \strong{Data Processing:}
#' \itemize{
#'   \item \strong{Z-score:} Calculated as \eqn{(x - \mu) / \sigma} per oxylipin.
#'   \item \strong{Winsorization (z_cap):} Limits extreme Z-scores to a specified
#'   range (default \eqn{\pm 4}) to prevent outliers from compressing the color scale.
#'   \item \strong{Dual-Matrix Logic:} Hierarchical clustering is performed on a
#'   fully imputed matrix to ensure stable dendrograms. For visualization,
#'   original \code{NA} positions are restored and displayed in a distinct color
#'   (default \code{"grey90"}), unless \code{hide_na = TRUE}.
#' }
#'
#' \strong{Clustering:} Uses Euclidean distance and the 'complete' linkage
#' method for both rows (oxylipins) and columns (samples).
#'
#' @param obj A \code{staRoxy} object.
#' @param na_method Imputation method for the clustering pipeline (default \code{"minprob"}).
#' @param na_threshold Missingness threshold for feature retention (default \code{0.5}).
#' @param remove_exclusive Logical. If \code{TRUE}, excludes lipids unique to one group.
#' @param colors A named character vector for group annotations.
#' @param show_values Logical. If \code{TRUE}, prints Z-score values inside cells.
#' @param hide_na Logical. If \code{TRUE}, filters out oxylipins with any missing values
#' before plotting.
#' @param na_color Color for missing values (default \code{"grey90"}).
#' @param z_cap Numeric. Maximum absolute Z-score to display (default \code{4}).
#' @param palette Character. Viridis palette option (e.g., "magma", "viridis").
#' @param x_axis_size,y_axis_size,legend_title_size,legend_size Numeric. Font sizes
#' for plot labels and legends.
#'
#' @return Draws a \code{ComplexHeatmap} and invisibly returns the \code{HeatmapList} object.
#'
#' @importFrom stats sd hclust dist
#' @importFrom circlize colorRamp2
#' @importFrom viridis viridis
#' @importFrom grid gpar grid.text
#' @importFrom ComplexHeatmap HeatmapAnnotation Heatmap Legend draw
#' @importFrom scales hue_pal
#' @importFrom cli cli_alert_danger
#'
#' @examples
#' # Generate heatmap with default settings
#' plot_heatmap_oxy(data_oxy_pellet)
#'
#' # Heatmap hiding NAs and capping Z-score at 2
#' plot_heatmap_oxy(data_oxy_pellet, hide_na = TRUE, z_cap = 2)
#'
#' @export
plot_heatmap_oxy <- function(obj,
                             na_method = "minprob",
                             na_threshold = 0.5,
                             remove_exclusive = TRUE,
                             colors = NULL,
                             show_values = FALSE,
                             hide_na = FALSE,
                             na_color = "grey90",
                             z_cap = 4,
                             palette = "magma",
                             x_axis_size = 12,
                             y_axis_size = 12,
                             legend_title_size = 12,
                             legend_size = 12) {

  # 1. Data Preparation and Filtering
  # Impute missing values for statistical consistency (clustering/z-score)
  m_imp <- prep_data_for_stats(
    obj,
    na_threshold = na_threshold,
    method = na_method,
    remove_exclusive = remove_exclusive
  )

  # Remove features with near-zero variance
  valid_rows <- apply(m_imp, 1, function(x) sd(x, na.rm = TRUE) > 1e-4)
  m_imp <- m_imp[valid_rows, , drop = FALSE]

  # Match the raw data matrix for visual NA preservation
  m_raw <- as.matrix(obj$data)[rownames(m_imp), , drop = FALSE]

  # 2. Z-score Calculation
  # Standardize rows (oxylipins) based on the imputed matrix
  z <- t(apply(m_imp, 1, function(x) (x - mean(x)) / sd(x)))

  # Apply capping (winsorization) to prevent outliers from dominating the scale
  z[z > z_cap] <- z_cap
  z[z < -z_cap] <- -z_cap

  # 3. Visual Logic: Re-introducing Missing Values
  # Map NAs from the original data back to the Z-score matrix for visualization
  z_vis <- z
  z_vis[is.na(m_raw)] <- NA

  # Optional: Remove any oxylipin that contains a missing value
  if (hide_na) {
    keep_vis <- rowSums(is.na(z_vis)) == 0
    z_vis <- z_vis[keep_vis, , drop = FALSE]
    z <- z[keep_vis, , drop = FALSE]

    if (nrow(z_vis) == 0) {
      return(cli::cli_alert_danger("No oxylipins left after visual NA filtering!"))
    }
  }

  # 4. Hierarchical Clustering
  # Cluster based on imputed data to ensure a continuous dendrogram
  hr <- hclust(dist(z), method = "complete")
  hc <- hclust(dist(t(z)), method = "complete")

  # 5. Aesthetics and Color Mapping
  lim <- max(abs(z), na.rm = TRUE)
  col_fun <- circlize::colorRamp2(
    seq(-lim, lim, length.out = 100),
    viridis::viridis(100, option = palette)
  )

  # Define text styles using grid graphical parameters
  style_t <- grid::gpar(fontsize = legend_title_size, fontface = "plain")
  style_l <- grid::gpar(fontsize = legend_size)

  # Setup group colors for annotation
  grps <- levels(as.factor(obj$meta$group))
  if (is.null(colors)) {
    colors <- scales::hue_pal()(length(grps))
    names(colors) <- grps
  } else {
    if (is.null(names(colors))) names(colors) <- grps
  }

  # 6. Heatmap Components
  # Create top annotation for experimental groups
  anno <- ComplexHeatmap::HeatmapAnnotation(
    Group = obj$meta$group,
    col = list(Group = colors),
    show_annotation_name = FALSE,
    annotation_legend_param = list(Group = list(title_gp = style_t, labels_gp = style_l))
  )

  # Build the main heatmap object
  ht <- ComplexHeatmap::Heatmap(
    z_vis,
    name = "Z-score",
    col = col_fun,
    na_col = na_color,
    top_annotation = anno,
    cluster_rows = hr,
    cluster_columns = hc,
    row_names_gp = grid::gpar(fontsize = y_axis_size),
    column_names_gp = grid::gpar(fontsize = x_axis_size),
    heatmap_legend_param = list(title_gp = style_t, labels_gp = style_l),
    cell_fun = function(j, i, x, y, w, h, f) {
      if (show_values && !is.na(z_vis[i, j])) {
        grid::grid.text(
          sprintf("%.1f", z_vis[i, j]), x, y,
          gp = grid::gpar(
            fontsize = 7,
            col = ifelse(abs(z_vis[i, j]) > lim / 2, "white", "black")
          )
        )
      }
    }
  )

  # 7. Final Rendering
  # Draw heatmap with or without the 'Missing Value' legend item
  if (!hide_na && any(is.na(z_vis))) {
    ComplexHeatmap::draw(
      ht,
      annotation_legend_list = list(
        ComplexHeatmap::Legend(
          labels = "Not Detected",
          title = "Missing",
          legend_gp = grid::gpar(fill = na_color),
          title_gp = style_t,
          labels_gp = style_l
        )
      )
    )
  } else {
    ComplexHeatmap::draw(ht)
  }
}
