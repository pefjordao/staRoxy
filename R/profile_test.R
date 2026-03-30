#' @title Multivariate Oxylipin Profile Analysis (PERMANOVA/PCoA)
#'
#' @description
#' Performs a multivariate comparison of global oxylipin profiles across
#' experimental groups. It uses PERMANOVA to test for statistical significance
#' and PCoA (Principal Coordinates Analysis) for visualization of beta-diversity.
#'
#' @details
#' \strong{Statistical Framework:}
#' \itemize{
#'   \item \strong{Standardization:} Data is Z-score scaled (mean=0, var=1) to
#'   ensure that high-abundance oxylipins do not overshadow biologically
#'   relevant lower-abundance mediators.
#'   \item \strong{PERMANOVA:} Executes a non-parametric Permutational
#'   Multivariate Analysis of Variance using Euclidean distances. The \eqn{R^2}
#'   indicates the proportion of total variance explained by the group structure.
#'   \item \strong{Post-hoc Testing:} If more than two groups are present and
#'   the global test is significant (\eqn{P-value < 0.05}), it performs pairwise
#'   PERMANOVA with Bonferroni correction.
#' }
#'
#' \strong{Visualization:}
#' Generates a PCoA plot (equivalent to PCA when using Euclidean distance)
#' showcasing the coordinates that explain the most variation in the
#' multidimensional lipid space.
#'
#' @param obj A \code{staRoxy} object.
#' @param na_method Character. Imputation strategy (MNAR/MAR pipeline).
#' Default is \code{"minprob"}.
#' @param na_threshold Numeric. Proportion of missing values allowed.
#' Default is \code{0.5}.
#' @param remove_exclusive Logical. If \code{TRUE}, removes lipids unique
#' to a single group to ensure robust multivariate modeling.
#' @param colors A named character vector for group colors.
#' @param permutations Integer. Number of permutations for the PERMANOVA
#' test. Default is \code{10000}.
#' @param plot Logical. If \code{TRUE}, generates a PCoA plot with
#' confidence ellipses.
#' @param title_size,x_axis_size,y_axis_size,legend_title_size,legend_size Numeric.
#' Font sizes for plot elements.
#'
#' @return Invisibly returns the global PERMANOVA results (\code{adonis2} object).
#'
#' @importFrom stats sd dist cmdscale
#' @importFrom ggplot2 ggplot aes geom_vline geom_hline stat_ellipse geom_point scale_color_manual scale_fill_manual labs theme element_text
#' @importFrom cowplot theme_half_open background_grid
#' @importFrom scales hue_pal
#' @importFrom cli cli_h1 cli_h2 cli_alert_info cli_rule
#'
#' @export
profile_test <- function(obj,
                         na_method = "minprob",
                         na_threshold = 0.5,
                         remove_exclusive = TRUE,
                         colors = NULL,
                         permutations = 10000,
                         plot = FALSE,
                         title_size = 14,
                         x_axis_size = 12,
                         y_axis_size = 12,
                         legend_title_size = 12,
                         legend_size = 12) {

  # 1. Dependencies & Data Preparation
  if (!requireNamespace("vegan", quietly = TRUE)) stop("Package 'vegan' required.")

  # Impute and clean data based on provided arguments
  m_imp <- prep_data_for_stats(
    obj,
    na_threshold = na_threshold,
    method = na_method,
    remove_exclusive = remove_exclusive
  )

  data_mat <- t(m_imp)
  meta_df  <- obj$meta

  # Filter out features with near-zero variance
  sds      <- apply(data_mat, 2, sd)
  data_mat <- data_mat[, sds > 1e-4, drop = FALSE]

  # Scale data for Euclidean-based PERMANOVA
  data_mat_scaled <- scale(data_mat, center = TRUE, scale = TRUE)

  # 2. Global PERMANOVA Execution
  res_glob <- vegan::adonis2(
    data_mat_scaled ~ group,
    data = meta_df,
    method = "euclidean",
    permutations = permutations
  )

  p_val <- round(res_glob$`Pr(>F)`[1], 3)
  r2    <- round(res_glob$R2[1], 3)

  # 3. CLI Reporting
  cli::cli_h1("Multivariate Profile Analysis (PERMANOVA)")
  cli::cli_alert_info("Groups: {paste(unique(meta_df$group), collapse = ', ')} | Oxylipins: {ncol(data_mat)} | Permutations: {permutations}")
  cli::cli_alert_info("Total Variance Explained (R²): {r2} | P-value: {p_val}")

  # 4. Post-hoc Pairwise Comparisons
  if (p_val < 0.05 && length(unique(meta_df$group)) > 2) {
    if (requireNamespace("pairwiseAdonis", quietly = TRUE)) {
      cli::cli_h2("Pairwise Comparisons (Bonferroni corrected)")

      pw <- suppressWarnings(pairwiseAdonis::pairwise.adonis(
        data_mat_scaled,
        meta_df$group,
        sim.method = "euclidean",
        p.adjust.m = "bonferroni"
      ))

      pw$R2         <- round(pw$R2, 3)
      pw$p.value    <- round(pw$p.value, 3)
      pw$p.adjusted <- round(pw$p.adjusted, 3)

      print(pw[, c("pairs", "R2", "p.value", "p.adjusted")])
    }
  }
  cli::cli_rule()

  if (!plot) return(invisible(res_glob))

  # 5. PCoA Calculation
  dist_obj <- dist(data_mat_scaled)
  pcoa_res <- cmdscale(dist_obj, k = 2, eig = TRUE)

  # Calculate variance explained per axis
  var_exp <- round(pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)

  pcoa_df <- data.frame(
    PCoA1 = pcoa_res$points[, 1],
    PCoA2 = pcoa_res$points[, 2],
    group = meta_df$group
  )

  if (is.null(colors)) colors <- scales::hue_pal()(length(unique(pcoa_df$group)))

  # 6. Data Visualization
  ggplot(pcoa_df, aes(x = PCoA1, y = PCoA2, color = group, fill = group)) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.3) +
    stat_ellipse(geom = "polygon", alpha = 0.15, show.legend = FALSE) +
    geom_point(size = 3, alpha = 0.8) +
    scale_color_manual(values = colors, name = "Group") +
    scale_fill_manual(values = colors, name = "Group") +
    labs(
      title = "PCoA: Oxylipin Profile Beta-Diversity",
      subtitle = bquote("PERMANOVA P =" ~ .(p_val) ~ "|" ~ "R2" == .(r2)),
      x = paste0("PCo 1 (", var_exp[1], "%)"),
      y = paste0("PCo 2 (", var_exp[2], "%)")
    ) +
    theme_half_open() +
    background_grid() +
    theme(
      plot.title    = element_text(size = title_size, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = element_text(size = title_size - 2, hjust = 0.5, color = "black"),
      axis.title    = element_text(size = x_axis_size),
      axis.text     = element_text(size = y_axis_size),
      legend.text   = element_text(size = legend_size),
      legend.title  = element_text(size = legend_title_size)
    )
}
