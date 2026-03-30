#' @title Loading Pattern Dissimilarity Analysis
#'
#' @description
#' Evaluates the structural stability of oxylipin profiles between two groups
#' by comparing their contributions (loadings) to a specific Principal Component.
#' It identifies whether the biological "architecture" of lipid relationships
#' is preserved or reorganized between conditions.
#'
#' @details
#' The function performs independent PCAs for each group and extracts the
#' loading vectors (\eqn{L_{g,i}}). Similarity is quantified by the squared
#' correlation (\eqn{R^2}), which represents global congruence.
#'
#' \strong{Mathematical Framework:}
#' \itemize{
#'   \item \strong{Sign Reflection:} Since Eigenvector signs are arbitrary, the
#'   function automatically reflects axes if \eqn{cor(L_1, L_2) < 0} to maximize
#'   congruence.
#'   \item \strong{Divergence:} Calculated as \eqn{1 - R^2}.
#'   \item \strong{Significance:} Assessed via Fisher Z-transformation:
#'   \deqn{Z = \frac{atanh(r_{obs}) - atanh(r_{stability})}{1/\sqrt{n-3}}}
#'   The P-value represents the probability of observing the profile correlation
#'   equal to or lower than \eqn{r_{obs}} under the null hypothesis of stability
#'   (\eqn{r \ge r_{stability}}).
#' }
#'
#' \strong{Visualization:} Generates a quadrant plot where oxylipins in
#' Quadrants II and IV represent features with the strongest divergence in
#' contribution between groups.
#'
#' @param obj A \code{staRoxy} object.
#' @param group1,group2 Character. The experimental groups to compare.
#' @param pc_index Integer. The Principal Component index to analyze (default is 1).
#' @param na_method Character. Imputation method for the QC pipeline.
#' @param na_threshold Numeric. Missingness threshold (default 0.5).
#' @param remove_exclusive Logical. If \code{TRUE}, excludes lipids unique to one group.
#' @param method Character. Correlation method: \code{"auto"} (default), \code{"pearson"}, or \code{"spearman"}.
#' @param threshold Numeric. \eqn{\Delta}Loading threshold for labeling divergent oxylipins.
#' @param r_stability Numeric. Stability threshold for the Z-test (default 0.6).
#' @param plot Logical. If \code{TRUE}, generates a scatter plot.
#' @param palette Character. Viridis palette for the \eqn{\Delta}Loading color scale.
#' @param title_size,subtitle_size,axis_title_size,axis_text_size,legend_title_size,legend_size,base_size Numeric.
#' Plot aesthetics and font sizes.
#'
#' @return Invisibly returns a data frame containing the loadings for both groups,
#' the calculated \eqn{\Delta}Loading, and metadata for plotting.
#'
#' @importFrom stats sd prcomp shapiro.test quantile IQR cor pnorm
#' @importFrom dplyr inner_join rename
#' @importFrom tidyr drop_na
#' @importFrom ggplot2 ggplot aes annotate geom_vline geom_hline geom_point geom_smooth scale_color_viridis_c labs theme element_text .pt
#' @importFrom cowplot theme_half_open background_grid
#' @importFrom ggrepel geom_text_repel
#' @importFrom cli cli_h1 cli_alert_info cli_alert_warning cli_rule
#'
#' @export
loading_dissimilarity <- function(obj,
                                  group1,
                                  group2,
                                  pc_index = 1,
                                  na_method = "minprob",
                                  na_threshold = 0.5,
                                  remove_exclusive = TRUE,
                                  method = "auto",
                                  threshold = 0.3,
                                  r_stability = 0.6,
                                  plot = TRUE,
                                  palette = "magma",
                                  title_size = 14,
                                  subtitle_size = 12,
                                  axis_title_size = 12,
                                  axis_text_size = 12,
                                  legend_title_size = 12,
                                  legend_size = 12,
                                  base_size = 12) {

  # 1. Global Data Preparation
  # Impute missing values globally to ensure comparable feature spaces
  m_imp_global <- prep_data_for_stats(
    obj,
    na_threshold = na_threshold,
    method = na_method,
    remove_exclusive = remove_exclusive
  )

  # Helper function to process PCA per specific group
  process_pca <- function(grp) {
    s <- obj$meta$sample[obj$meta$group == grp]
    m_c <- t(m_imp_global[, s, drop = FALSE])

    # Variance filter: essential for group-specific PCA stability
    sds <- apply(m_c, 2, stats::sd)
    m_c <- m_c[, sds > 1e-4, drop = FALSE]

    res <- stats::prcomp(m_c, center = TRUE, scale. = TRUE)
    list(pca = res, var = (res$sdev^2) / sum(res$sdev^2) * 100)
  }

  # 2. Comparative PCA Execution
  r1 <- process_pca(group1)
  r2 <- process_pca(group2)

  # Helper to extract and normalize loadings (force positive orientation)
  extract_l <- function(p, i) {
    ld <- data.frame(Lipid = rownames(p$rotation), Weight = p$rotation[, i])
    if (sum(ld$Weight) < 0) ld$Weight <- -ld$Weight
    return(ld)
  }

  # Join loadings from both groups for comparison
  df_c <- dplyr::inner_join(
    extract_l(r1$pca, pc_index) %>% dplyr::rename(!!group1 := Weight),
    extract_l(r2$pca, pc_index) %>% dplyr::rename(!!group2 := Weight),
    by = "Lipid"
  ) %>% tidyr::drop_na()

  # 3. Congruence Statistics
  # Normality check (Shapiro-Wilk)
  is_n <- stats::shapiro.test(df_c[[group1]])$p.value > 0.05 &&
    stats::shapiro.test(df_c[[group2]])$p.value > 0.05

  # Outlier detection logic (IQR method)
  det_out <- function(x) {
    q <- stats::quantile(x, probs = c(0.25, 0.75))
    h <- 1.5 * stats::IQR(x)
    any(x < (q[1] - h) | x > (q[2] + h))
  }
  has_out <- det_out(df_c[[group1]]) || det_out(df_c[[group2]])

  # Select correlation method based on data distribution
  if (method == "auto") {
    method <- if (is_n && !has_out) "pearson" else "spearman"
  }

  r_init <- stats::cor(df_c[[group1]], df_c[[group2]], method = method)

  # Handle sign reflection (eigenvector orientation)
  flipped <- FALSE
  if (r_init < 0) {
    df_c[[group2]] <- df_c[[group2]] * -1
    flipped <- TRUE
    r_obs <- stats::cor(df_c[[group1]], df_c[[group2]], method = method)
  } else {
    r_obs <- r_init
  }

  # Divergence metrics and Significance (Z-test vs stability threshold)
  n <- nrow(df_c)
  r_sq <- r_obs^2
  divergence <- 1 - r_sq
  z_score <- (atanh(r_obs) - atanh(r_stability)) / (1 / sqrt(n - 3))
  p_div <- stats::pnorm(z_score)
  p_div_f <- if (p_div < 0.001) "< 0.001" else round(p_div, 3)

  # 4. CLI Reporting
  cli::cli_h1("Loading Pattern Dissimilarity Analysis")
  cli::cli_alert_info("Comparison: '{group1}' vs '{group2}' (PC{pc_index})")
  cli::cli_alert_info("Method: {toupper(method)} (Normal={is_n}, Outliers={has_out})")

  if (flipped) cli::cli_alert_warning("Reflection: Eigenvector sign was inverted to maximize congruence.")

  cli::cli_alert_info("Correlation (R): {round(r_obs, 3)} | Global Congruence (R²): {round(r_sq * 100, 1)}%")
  cli::cli_alert_info("Divergence (1-R²): {round(divergence * 100, 1)}% | Significance: P(div) = {p_div_f}")
  cli::cli_rule()

  if (!plot) return(invisible(df_c))

  # 5. Visualization (ggplot2)
  ggplot2::ggplot(df_c, ggplot2::aes(x = .data[[group1]], y = .data[[group2]])) +
    # Background rects to highlight sign shifts
    ggplot2::annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = "grey90", alpha = 0.4) +
    ggplot2::annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = "grey90", alpha = 0.4) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.2) +
    ggplot2::geom_point(ggplot2::aes(color = abs(.data[[group1]] - .data[[group2]])), size = 3, alpha = 0.7) +
    ggplot2::geom_smooth(method = "lm", color = "black", linetype = "dashed", se = FALSE, formula = y ~ x, linewidth = 1) +
    # Labeling top divergent oxylipins
    ggrepel::geom_text_repel(
      ggplot2::aes(label = ifelse(abs(.data[[group1]] - .data[[group2]]) >= threshold, Lipid, "")),
      size = (axis_text_size * 0.9 / ggplot2::.pt),
      max.overlaps = 20
    ) +
    ggplot2::scale_color_viridis_c(option = palette, name = expression(Delta * Loading)) +
    ggplot2::labs(
      title = "Loading Pattern Dissimilarity",
      subtitle = paste0("PC", pc_index, " | Divergence: ", round(divergence * 100, 1), "% | P(div) = ", p_div_f),
      x = paste0(group1, " Loading (", round(r1$var[pc_index], 1), "%)"),
      y = paste0(group2, " Loading (", round(r2$var[pc_index], 1), "%)")
    ) +
    cowplot::theme_half_open() +
    cowplot::background_grid() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(size = subtitle_size, hjust = 0.5),
      axis.title = ggplot2::element_text(size = axis_title_size),
      axis.text = ggplot2::element_text(size = axis_text_size),
      legend.text = ggplot2::element_text(size = legend_size),
      legend.title = ggplot2::element_text(size = legend_title_size)
    )
}
