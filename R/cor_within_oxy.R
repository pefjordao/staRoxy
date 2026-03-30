#' @title Intra-group Correlation Analysis Between Two Lipids
#'
#' @description
#' Computes the correlation between two different oxylipins within a single
#' experimental group. Useful for identifying metabolic co-regulation or
#' pathway associations.
#'
#' @details
#' Like its counterpart, this function uses an automatic method selection
#' (Pearson/Spearman) based on data distribution.
#'
#' @param obj A \code{staRoxy} object.
#' @param oxylipin1 Character. The first lipid name.
#' @param oxylipin2 Character. The second lipid name.
#' @param group Character. The experimental group to analyze.
#' @param method Character. "auto" (default), "pearson", or "spearman".
#' @param plot Logical. If \code{TRUE}, returns a scatter plot.
#' @param colors A named character vector for group coloring.
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
#' # Correlate PGE2 and 9-HETE within the Scraper group
#' cor_within_oxy(data_oxy_pellet,
#'                oxylipin1 = "PGE2",
#'                oxylipin2 = "9-HETE",
#'                group = "Scraper")
#'
#' @export
cor_within_oxy <- function(obj,
                           oxylipin1,
                           oxylipin2,
                           group,
                           method = "auto",
                           plot = TRUE,
                           colors = NULL,
                           title_size = 14,
                           x_axis_size = 12,
                           y_axis_size = 12) {

  # 1. Data Validation
  if (!(oxylipin1 %in% rownames(obj$data)) || !(oxylipin2 %in% rownames(obj$data))) {
    return(cli::cli_alert_danger("Lipids not found in the dataset."))
  }

  # Extract group-specific values
  idx     <- which(obj$meta$group == group)
  x_vals  <- as.numeric(obj$data[oxylipin1, idx])
  y_vals  <- as.numeric(obj$data[oxylipin2, idx])

  # Remove non-finite values
  valid   <- is.finite(x_vals) & is.finite(y_vals)
  x_vals  <- x_vals[valid]
  y_vals  <- y_vals[valid]

  # 2. Statistical Correlation
  # Run advisor to get flags and method
  method_res <- check_cor_method(x_vals, y_vals, silent = TRUE)
  is_normal  <- attr(method_res, "is_normal")
  has_out    <- attr(method_res, "has_out")

  if (method == "auto") method <- as.character(method_res)

  cor_res <- tryCatch(
    stats::cor.test(x_vals, y_vals, method = method),
    error = function(e) NULL
  )

  if (is.null(cor_res)) {
    return(cli::cli_alert_danger("Insufficient data for correlation in group '{group}'."))
  }

  # Extract metrics
  r <- as.numeric(cor_res$estimate)
  p <- as.numeric(cor_res$p.value)
  interpretation <- get_cor_label(r, p)

  # 3. Console Reporting
  cli::cli_h1("Intra-group Correlation Analysis")
  cli::cli_alert_info("Comparison: '{oxylipin1}' vs. '{oxylipin2}' | Group: '{group}'")
  cli::cli_alert_info("Method: {toupper(method)} (Normal={is_normal}, Outliers={has_out})")

  if (p < 0.05) {
    cli::cli_alert_success("R = {round(r, 3)} | P = {format.pval(p, digits = 3)} | Result: {interpretation}")
  } else {
    cli::cli_alert_info("R = {round(r, 3)} | P = {format.pval(p, digits = 3)} | Result: {interpretation}")
  }
  cli::cli_rule()

  if (!plot) return(invisible(NULL))

  # 4. Visualization
  line_color <- if (!is.null(colors) && group %in% names(colors)) colors[group] else "grey30"
  plot_df    <- data.frame(L1 = x_vals, L2 = y_vals)

  ggplot2::ggplot(plot_df, ggplot2::aes(x = L1, y = L2)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, color = line_color, fill = line_color, alpha = 0.15, linewidth = 0.5) +
    ggplot2::geom_point(shape = 21, fill = line_color, color = "black", size = 3, alpha = 0.8) +
    ggplot2::labs(
      title    = paste("Within:", group),
      subtitle = paste0(interpretation, "\n", stringr::str_to_title(method), " | R = ", round(r, 3), " | P = ", round(p, 3)),
      x        = oxylipin1, y = oxylipin2
    ) +
    cowplot::theme_half_open() +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold", size = title_size),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = title_size - 2),
      axis.title    = ggplot2::element_text(size = x_axis_size),
      axis.text     = ggplot2::element_text(size = y_axis_size)
    )
}
