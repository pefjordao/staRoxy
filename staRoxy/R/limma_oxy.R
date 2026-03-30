#' @title Differential Abundance Analysis of Oxylipins
#'
#' @description
#' Fits linear models to oxylipin abundance data and applies Empirical Bayes
#' moderation to identify differentially abundant features across experimental
#' groups.
#'
#' @details
#' The function utilizes the \code{limma} framework, ideal for datasets with
#' small sample sizes. It automatically generates all possible pairwise
#' comparisons (contrasts). Missing values are handled via a multi-step
#' pipeline (MNAR vs. MAR) integrated into the internal data preparation.
#'
#' @param obj A \code{staRoxy} object.
#' @param covar Character. Optional covariates for the model (e.g., "Age").
#' Default is \code{NULL}.
#' @param na_method Character. Imputation method: \code{"minprob"} (default)
#' or others supported by \code{prep_data_for_stats}.
#' @param na_threshold Numeric. Missingness threshold (0-1). Default is \code{0.5}.
#' @param remove_exclusive Logical. If \code{TRUE}, removes group-exclusive lipids.
#' @param save_results Logical. If \code{TRUE}, exports CSV files for each contrast.
#'
#' @return A named list of data frames containing Log2FC, P-values, and FDR.
#'
#' @importFrom limma lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom stats setNames model.matrix as.formula
#' @importFrom utils combn write.csv2
#' @importFrom cli cli_alert_success
#'
#' @export
limma_oxy <- function(obj,
                      covar = NULL,
                      na_method = "minprob",
                      na_threshold = 0.5,
                      remove_exclusive = TRUE,
                      save_results = TRUE) {

  # 1. Data Preparation
  # Prepare data using the standardized imputation and filtering pipeline
  m_prep <- prep_data_for_stats(
    obj,
    na_threshold = na_threshold,
    method = na_method,
    remove_exclusive = remove_exclusive
  )

  # 2. Syntactic Protection for Group Names
  # Ensure group levels are safe for model formulas (e.g., removing spaces)
  orig_levels <- levels(as.factor(obj$meta$group))
  safe_levels <- make.names(orig_levels)
  name_map <- stats::setNames(orig_levels, safe_levels)

  meta_temp <- obj$meta
  meta_temp$group <- factor(make.names(meta_temp$group), levels = safe_levels)

  # 3. Model Configuration
  # Build formula: including covariates if provided
  f_string <- if (is.null(covar)) "~ 0 + group" else paste0("~ 0 + group + ", covar)
  des <- stats::model.matrix(as.formula(f_string), data = meta_temp)

  grps <- levels(meta_temp$group)
  colnames(des)[1:length(grps)] <- grps

  # 4. Limma Differential Analysis
  # Linear modeling and Empirical Bayes moderation
  suppressWarnings({
    fit <- limma::lmFit(m_prep, des)

    # Generate all pairwise combinations for contrasts
    cont_matrix <- apply(utils::combn(grps, 2), 2, function(x) paste0(x[1], "-", x[2]))

    cont <- limma::makeContrasts(contrasts = cont_matrix, levels = des)
    fit_e <- limma::eBayes(limma::contrasts.fit(fit, cont))
  })

  # 5. Result Processing and Label Translation
  # Map internal safe names back to original user-friendly group names
  res <- lapply(1:ncol(fit_e$contrasts), function(i) {
    df <- limma::topTable(fit_e, coef = i, number = Inf)

    comp_safe <- colnames(fit_e$contrasts)[i]
    parts <- unlist(strsplit(comp_safe, "-"))

    # Formatted comparison name: "GroupA vs. GroupB"
    comp_pretty <- paste0(name_map[parts[1]], " vs. ", name_map[parts[2]])

    attr(df, "comp") <- comp_pretty
    return(df)
  })

  names(res) <- sapply(res, function(x) attr(x, "comp"))

  # 6. Results Export (CSV)
  if (save_results) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    for (i in seq_along(res)) {
      pretty_name <- names(res)[i]

      # File naming: replace spaces with underscores for compatibility
      safe_fname <- gsub(" ", "_", pretty_name)
      fname <- paste0("limma_", safe_fname, "_", ts, ".csv")

      df_to_save <- cbind(Oxylipin = rownames(res[[i]]), res[[i]])
      utils::write.csv2(df_to_save, file = fname, row.names = FALSE, fileEncoding = "UTF-8")
    }
  }

  cli::cli_alert_success("limma differential analysis complete.")

  return(res)
}
