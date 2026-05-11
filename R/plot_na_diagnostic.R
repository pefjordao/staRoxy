#' @title Missing Data Diagnostic for Oxylipins
#'
#' @description
#' Visualizes missing value patterns to distinguish between Missing Not At Random
#' (MNAR), typically associated with the Limit of Detection (LOD), and
#' Missing At Random (MAR), arising from experimental or technical variability.
#'
#' @details
#' The function provides two diagnostic modes:
#' \itemize{
#'   \item \strong{Correlation mode ("cor"):} Generates a scatter plot of Mean
#'   Log2 Abundance vs. Proportion of Missing Values. A negative correlation
#'   suggests an MNAR pattern, where lipids are more likely to be missing
#'   at lower abundances.
#'   \item \strong{Distribution mode ("dist"):} Compares the abundance
#'   distributions of features with and without missing values using Density
#'   and ECDF plots. A leftward shift for features with NAs (\code{TRUE})
#'   is a strong indicator of MNAR.
#' }
#'
#' \strong{Implication for Imputation:}
#' \itemize{
#'   \item If MNAR is detected, use \code{na_method = "minprob"} in downstream functions.
#'   \item If MAR is detected, use \code{na_method = "rf"} (Random Forest).
#' }
#'
#' @param obj A \code{staRoxy} object.
#' @param mode Character. Either \code{"cor"} (default) or \code{"dist"}.
#' @param point_color,true_color,false_color,line_color Character. Hex codes or names for plot colors.
#' @param title_size,x_axis_size,y_axis_size,legend_title_size,legend_size Numeric. Font sizes for plot elements.
#'
#' @return A \code{ggplot2} object (mode \code{"cor"}) or a \code{cowplot} grid object (mode \code{"dist"}).
#'
#' @examples
#' \dontshow{
#' staRoxy_object <- read_oxy(data_oxy_lps_pellet, metadata_oxy_lps_pellet)
#' staRoxy_object <- filter_oxy(staRoxy_object)
#' staRoxy_object <- transform_oxy(staRoxy_object)
#' }
#'
#'\dontrun{
#' # Diagnostic using correlation mode (default)
#' # Ideal to see if lower abundance correlates with more missing values
#' plot_na_diagnostic(staRoxy_object, mode = "cor")
#'
#' # Diagnostic using distribution mode
#' # Ideal to compare density and cumulative distribution shifts
#' plot_na_diagnostic(staRoxy_object, mode = "dist")
#'}
#'
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by summarise n mutate ungroup filter distinct
#' @importFrom ggplot2 ggplot aes theme element_text geom_point geom_smooth labs geom_density scale_color_manual stat_ecdf
#' @importFrom cowplot theme_half_open background_grid plot_grid
#' @importFrom stats sd median
#' @importFrom cli cli_alert_info
#'
#' @export
plot_na_diagnostic <- function(obj,
                               mode = "cor",
                               point_color = "#2C56A0",
                               true_color = "#C9273E",
                               false_color = "#57B894",
                               line_color = "black",
                               title_size = 14,
                               x_axis_size = 12,
                               y_axis_size = 12,
                               legend_title_size = 12,
                               legend_size = 12) {

  # Data Preparation
  # Reshape data to long format for easier ggplot mapping
  df <- as.data.frame(t(obj$data))
  df$group <- obj$meta$group

  df_long <- tidyr::pivot_longer(
    df,
    cols = -group,
    names_to = "oxylipin",
    values_to = "intensity"
  )

  # Helper to maintain staRoxy visual identity
  apply_staRoxy_theme <- function(plt) {
    plt + ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = 10, hjust = 0.5, color = "grey40"),
      axis.title = ggplot2::element_text(size = x_axis_size),
      axis.text = ggplot2::element_text(size = y_axis_size),
      legend.text = ggplot2::element_text(size = legend_size),
      legend.title = ggplot2::element_text(size = legend_title_size),
      legend.position = "right"
    )
  }

  # Mode: Correlation (Abundance vs. Missingness)
  if (mode == "cor") {
    na_stats <- df_long %>%
      dplyr::group_by(oxylipin) %>%
      dplyr::summarise(
        mean_int = mean(intensity, na.rm = TRUE),
        na_prop = sum(is.na(intensity)) / dplyr::n(),
        .groups = "drop"
      )

    p <- ggplot2::ggplot(na_stats, ggplot2::aes(x = mean_int, y = na_prop)) +
      cowplot::theme_half_open() +
      cowplot::background_grid() +
      ggplot2::geom_point(alpha = 0.8, color = point_color, size = 2)

    # TREND LINE CHECK: Only plot if there's variation in NAs across lipids
    if (stats::sd(na_stats$na_prop, na.rm = TRUE) > 0) {
      p <- p + ggplot2::geom_smooth(
        method = "gam",
        formula = y ~ s(x, bs = "cs"),
        color = line_color,
        se = TRUE,
        linetype = "dashed",
        linewidth = 1
      )
    } else {
      cli::cli_alert_info("No variation in missing data detected. Skipping trend line.")
    }

    p <- p + ggplot2::labs(
      title = "Missing Data Diagnostic - Correlation",
      x = "Mean Log2 Abundance",
      y = "Proportion of Missing Values"
    )

    return(apply_staRoxy_theme(p))

  }

  # Mode: Distribution Comparison
  else if (mode == "dist") {

    dist_data <- df_long %>%
      dplyr::group_by(oxylipin) %>%
      dplyr::mutate(
        na_count = sum(is.na(intensity)),
        has_na = na_count > 0
      ) %>%
      dplyr::ungroup()

    if (length(unique(dist_data$has_na)) == 1 && dist_data$has_na[1] == TRUE) {

      ox_na_summary <- dist_data %>% dplyr::distinct(oxylipin, na_count)
      median_na <- stats::median(ox_na_summary$na_count)

      cli::cli_alert_info("No oxylipins with 100% detection. Splitting distribution by median missingness (> {median_na} NAs).")

      dist_data <- dist_data %>%
        dplyr::mutate(
          plot_group = ifelse(na_count > median_na, "High Missingness", "Low/Avg Missingness")
        )
      legend_title <- "Missingness Level"

    } else {
      dist_data <- dist_data %>%
        dplyr::mutate(
          plot_group = ifelse(has_na, "Has NAs", "100% Detected")
        )
      legend_title <- "Missing Values"
    }

    dist_data <- dist_data %>% dplyr::filter(!is.na(intensity))

    color_mapping <- c(
      "High Missingness" = true_color,
      "Low/Avg Missingness" = false_color,
      "Has NAs" = true_color,
      "100% Detected" = false_color
    )

    # Density Plot
    p1 <- ggplot2::ggplot(dist_data, ggplot2::aes(x = intensity, color = plot_group)) +
      ggplot2::geom_density(linewidth = 1) +
      ggplot2::scale_color_manual(
        values = color_mapping,
        name = legend_title
      ) +
      ggplot2::labs(
        title = "Missing Data Diagnostic - Distribution",
        x = "Log2 Abundance",
        y = "Density"
      ) +
      cowplot::theme_half_open() +
      cowplot::background_grid()

    # ECDF Plot (Cumulative Distribution)
    p2 <- ggplot2::ggplot(dist_data, ggplot2::aes(x = intensity, color = plot_group)) +
      ggplot2::stat_ecdf(linewidth = 1) +
      ggplot2::scale_color_manual(
        values = color_mapping,
        name = legend_title
      ) +
      ggplot2::labs(
        x = "Log2 Abundance",
        y = "Cumulative Fraction"
      ) +
      cowplot::theme_half_open() +
      cowplot::background_grid()

    p1 <- apply_staRoxy_theme(p1)
    p2 <- apply_staRoxy_theme(p2)

    # Stack plots vertically
    return(cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 0.8)))

  } else {
    stop("Invalid mode. Choose 'cor' or 'dist'.")
  }
}
