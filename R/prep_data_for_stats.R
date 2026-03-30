#' @title Internal Pre-processing and Imputation Pipeline
#'
#' @description
#' A robust internal pipeline for filtering and imputing missing oxylipidomics
#' data. It distinguishes between technical missingness (MAR) and biological
#' non-detection (MNAR) to ensure statistical stability.
#'
#' @details
#' \strong{1. Filtering Phase:}
#' \itemize{
#'   \item \strong{Threshold:} Retains oxylipins with valid measurements in at
#'   least \eqn{(1 - na\_threshold)} proportion of samples in at least one group.
#'   \item \strong{Exclusivity:} If \code{remove_exclusive = TRUE}, oxylipins
#'   missing 100% of values in any group are removed to avoid artificial
#'   variance stabilization issues in linear models.
#' }
#'
#' \strong{2. Imputation Phase:}
#' \itemize{
#'   \item \strong{MinProb (MNAR):} For Missing Not At Random values. Imputes
#'   group-wise missingness using the group mean. If an entire group is missing,
#'   it uses the global minimum value for that lipid, adding stochastic noise
#'   (10% of global SD) to mimic the Limit of Detection (LOD).
#'   \item \strong{Random Forest (MAR):} For Missing At Random values. Leverages
#'   the global metabolic profile via the \code{missForest} package to
#'   estimate stochastic missingness.
#' }
#'
#' @param obj A \code{staRoxy} object.
#' @param na_threshold Numeric. Maximum proportion of missing values allowed
#' (0 to 1). Default is \code{0.5} (requires 50% valid values).
#' @param method Character. Imputation strategy: \code{"minprob"} (MNAR-focused)
#' or \code{"rf"} (Random Forest for MAR-focused).
#' @param remove_exclusive Logical. Whether to filter out lipids detected in
#' only one group. Default is \code{TRUE}.
#' @param seed Numeric. Seed for reproducibility of stochastic imputation.
#'
#' @return A numeric matrix of filtered and imputed log2-transformed data.
#'
#' @importFrom stats sd rnorm ave
#' @importFrom cli cli_alert_warning cli_alert_info cli_alert_success
#'
#' @keywords internal
prep_data_for_stats <- function(obj,
                                na_threshold = 0.5,
                                method = "minprob",
                                remove_exclusive = TRUE,
                                seed = 5742) {

  # 1. Initialization
  set.seed(seed)
  m <- as.matrix(obj$data)
  grps <- obj$meta$group
  n_start <- nrow(m)

  # 2. Filtering Phase

  # STEP 2.1: Threshold Filter (Minimum Survival)
  keep_thr <- apply(m, 1, function(x) {
    any(tapply(x, grps, function(g) sum(!is.na(g)) / length(g)) >= (1 - na_threshold))
  })

  n_filtered_thr <- sum(!keep_thr)
  if (n_filtered_thr > 0) {
    cli::cli_alert_warning("Filtered out {n_filtered_thr} oxylipin{?s} failing the {na_threshold * 100}% threshold in all groups.")
  }
  m <- m[keep_thr, , drop = FALSE]

  # STEP 2.2: Exclusivity Filter
  if (remove_exclusive) {
    keep_exc <- apply(m, 1, function(x) {
      group_counts <- tapply(x, grps, function(g) sum(!is.na(g)))
      return(all(group_counts > 0))
    })

    n_filtered_exc <- sum(!keep_exc)
    if (n_filtered_exc > 0) {
      cli::cli_alert_warning("Removed {n_filtered_exc} oxylipin{?s} with zero variance/data in at least one group.")
    }
    m <- m[keep_exc, , drop = FALSE]
  } else {
    cli::cli_alert_info("Exclusivity filter disabled. Groups with 100% NAs will be fully imputed.")
  }

  if (nrow(m) == 0) {
    stop("All oxylipins were filtered out. Relax your 'na_threshold' or check data quality.")
  }

  n_final <- nrow(m)
  cli::cli_alert_success("Filtering complete: {n_final}/{n_start} oxylipins retained.")

  # 3. Imputation Phase

  if (method == "minprob") {
    cli::cli_alert_info("Imputing via 'minprob': assuming missingness is due to low concentration (MNAR).")

    m_imp <- t(apply(m, 1, function(x) {
      global_min <- min(x, na.rm = TRUE)
      global_sd  <- stats::sd(x, na.rm = TRUE)

      if (is.na(global_sd) || global_sd == 0) {
        global_sd <- abs(global_min) * 0.05
      }

      stats::ave(x, grps, FUN = function(g) {
        if (all(is.na(g))) {
          if (is.infinite(global_min)) return(rep(0, length(g)))
          return(stats::rnorm(length(g), mean = global_min, sd = global_sd * 0.1))
        }
        g[is.na(g)] <- mean(g, na.rm = TRUE)
        return(g)
      })
    }))

  } else if (method == "rf") {
    cli::cli_alert_info("Imputing via Random Forest: assuming missingness is stochastic (MAR).")

    if (!requireNamespace("missForest", quietly = TRUE)) {
      stop("Package 'missForest' is required for RF imputation.")
    }

    rf_res <- suppressWarnings(missForest::missForest(t(m)))
    m_imp <- t(rf_res$ximp)

  } else {
    stop("Invalid imputation method. Use 'minprob' or 'rf'.")
  }

  # 4. Finalization
  colnames(m_imp) <- colnames(m)
  rownames(m_imp) <- rownames(m)

  return(m_imp)
}
