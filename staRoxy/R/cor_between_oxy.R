#' @title Inter-group Correlation Analysis for a Single Lipid
#'
#' @description
#' Calculates the correlation for a specific oxylipin across two different
#' experimental groups. This is typically used for paired data or to assess
#' consistency of a single feature across conditions.
#'
#' @details
#' The function automatically selects between Pearson and Spearman methods
#' based on normality and outlier presence unless a method is explicitly
#' specified. It requires that both groups have the same number of samples
#' (paired design).
#'
#' @param obj A \code{staRoxy} object.
#' @param oxylipin Character. The name of the lipid to analyze.
#' @param group1 Character. The first experimental group.
#' @param group2 Character. The second experimental group.
#' @param method Character. "auto" (default), "pearson", or "spearman".
#' @param plot Logical. If \code{TRUE}, returns a scatter plot with a linear
#' regression line.
#' @param title_size,x_axis_size,y_axis_size Numeric. Font sizes for plot elements.
#'
#' @return Invisibly returns \code{NULL}. Prints statistics to the console
#' and generates a plot if \code{plot = TRUE}.
#'
#' @importFrom stats cor.test
#' @importFrom ggplot2 ggplot aes geom_smooth geom_point labs theme element_text
#' @importFrom cowplot theme_half_open
#' @importFrom stringr str_to_title
#' @importFrom cli cli_h1 cli_alert_info cli_alert_success cli_alert_danger cli_rule
#'
#' @examples
#' # Compare PGE2 levels between Scraper and Trypsin groups
#' cor_between_oxy(data_oxy_pellet,
#'                 oxylipin = "PGE2",
#'                 group1 = "Scraper",
#'                 group2 = "Trypsin")
#'
#' @export
cor_between_oxy <- function(obj,
                            oxylipin,
                            group1,
                            group2,
                            method = "auto",
                            plot = TRUE,
                            title_size = 14,
                            x_axis_size = 12,
                            y_axis_size = 12) {

  # 1. Data Selection and Paired Validation
  vals1 <- as.numeric(obj$data[oxylipin, obj$meta$group == group1])
  vals2 <- as.numeric(obj$data[oxylipin, obj$meta$group == group2])

  if (length(vals1) != length(vals2)) {
    return(cli::cli_alert_danger("Unequal sizes (check if samples are paired)."))
  }

  # Mathematical validity filter
  valid <- is.finite(vals1) & is.finite(vals2)
  vals1 <- vals1[valid]
  vals2 <- vals2[valid]

  # 2. Statistical Correlation
  # Run advisor to get flags and method
  method_res <- check_cor_method(vals1, vals2, silent = TRUE)
  is_normal  <- attr(method_res, "is_normal")
  has_out    <- attr(method_res, "has_out")

  if (method == "auto") method <- as.character(method_res)

  cor_res <- tryCatch(
    stats::cor.test(vals1, vals2, method = method),
    error = function(e) NULL
  )

  if (is.null(cor_res)) {
    return(cli::cli_alert_danger("Insufficient data for correlation analysis."))
  }

  # Extract metrics
  r <- as.numeric(cor_res$estimate)
  p <- as.numeric(cor_res$p.value)
  interpretation <- get_cor_label(r, p)

  # 3. Console Reporting
  cli::cli_h1("Inter-group Correlation Analysis")
  cli::cli_alert_info("Lipid: '{oxylipin}' | Comparison: '{group1}' vs. '{group2}'")
  cli::cli_alert_info("Method: {toupper(method)} (Normal={is_normal}, Outliers={has_out})")

  if (p < 0.05) {
    cli::cli_alert_success("R = {round(r, 3)} | P = {format.pval(p, digits = 3)} | Result: {interpretation}")
  } else {
    cli::cli_alert_info("R = {round(r, 3)} | P = {format.pval(p, digits = 3)} | Result: {interpretation}")
  }
  cli::cli_rule()

  if (!plot) return(invisible(NULL))

  # 4. Visualization
  plot_df <- data.frame(G1 = vals1, G2 = vals2)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = G1, y = G2)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = "grey30", fill = "grey30", alpha = 0.15, linewidth = 0.5) +
    ggplot2::geom_point(shape = 21, fill = "grey70", color = "black", size = 3, alpha = 0.8) +
    ggplot2::labs(
      title    = paste0("Between ", group1, " vs. ", group2),
      subtitle = paste0(interpretation, "\n", stringr::str_to_title(method), " | R = ", round(r, 3), " | P = ", round(p, 3)),
      x        = paste(oxylipin, group1), y = paste(oxylipin, group2)
    ) +
    cowplot::theme_half_open() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_size),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = title_size - 2),
      axis.title    = ggplot2::element_text(size = x_axis_size),
      axis.text     = ggplot2::element_text(size = y_axis_size)
    )
}
