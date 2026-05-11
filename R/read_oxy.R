#' @title Import and Harmonize Oxylipidomics Data
#'
#' @description
#' The primary entry point for the \code{staRoxy} pipeline. It imports abundance
#' data and metadata from various file formats (Excel, CSV, or TSV), standardizes
#' lipid nomenclature, and ensures perfect alignment between samples and their
#' experimental annotations.
#'
#' @details
#' \strong{Data Import and Cleaning:}
#' \itemize{
#'   \item \strong{File Formats:} Supports \code{.xlsx}, \code{.xls}, \code{.csv},
#'   and \code{.txt} (via \code{readxl} and \code{data.table}). It also directly accepts
#'   \code{data.frame} or \code{matrix} objects already loaded in R environment.
#'   \item \strong{Label Standardization:} Automatically applies \code{clean_labels()}
#'   to fix encoding artifacts and biochemical symbols in oxylipin names.
#'   \item \strong{Matrix Sanitization:} Automatically handles common data entry
#'   issues, such as converting commas to dots and treating zeros as missing values
#'   (\code{NA}).
#' }
#'
#' \strong{Sample Alignment:}
#' The function identifies the intersection of sample IDs between the abundance
#' matrix and the metadata. Any samples present in only one of the files are
#' automatically removed, and a summary report is printed to the console.
#'
#' @param data_path Character path to the oxylipin abundance file OR a loaded data.frame/matrix. The first column or rownames must contain oxylipin identifiers.
#' @param meta_path Character path to the metadata file OR a loaded data.frame. It must contain at least a \code{sample} column (matching the data headers) and a \code{group} column.
#' @param auto_clean Logical. If \code{TRUE} (default), applies internal nomenclature
#' cleaning to fix encoding artifacts common in oxylipidomics. Set to \code{FALSE}
#' if using staRoxy for other omics data (e.g., proteomics) to preserve original names.
#'
#' @return An object of class \code{staRoxy}, which is a list containing:
#' \itemize{
#'   \item \strong{data:} A numeric matrix of oxylipin abundances.
#'   \item \strong{meta:} A data frame of experimental metadata.
#'   \item \strong{info:} A list of summary statistics and creation date.
#' }
#'
#' @examples
#' \dontshow{
#' df_abundance <- data_oxy_lps_pellet
#' df_metadata <- metadata_oxy_lps_pellet
#' }
#'
#' # Loading from objects already in the R environment
#' staRoxy_object <- read_oxy(data_path = df_abundance, meta_path = df_metadata)
#'
#' # Loading directly from files (example syntax)
#' # staRoxy_object <- read_oxy(
#' #   data_path = "path/to/abundance_data.xlsx",
#' #   meta_path = "path/to/metadata.csv"
#' # )
#'
#' @importFrom tools file_ext
#' @importFrom readxl read_excel
#' @importFrom data.table fread
#' @importFrom cli cli_h1 cli_alert_success cli_alert_warning cli_alert_info cli_h2 cli_li cli_rule cli_text
#'
#' @export
read_oxy <- function(data_path, meta_path, auto_clean = TRUE) {

  # Internal Helper Functions
  load_any <- function(input) {
    if (is.data.frame(input) || is.matrix(input)) {
      df <- as.data.frame(input)
      if (!identical(rownames(df), as.character(seq_len(nrow(df))))) {
        df <- cbind(ID = rownames(df), df)
        rownames(df) <- NULL
      }
      return(df)
    }

    if (is.character(input)) {
      ext <- tolower(tools::file_ext(input))
      if (ext %in% c("xls", "xlsx")) {
        return(as.data.frame(readxl::read_excel(input)))
      }
      return(data.table::fread(input, data.table = FALSE, check.names = FALSE, encoding = "UTF-8"))
    }

    stop("Input must be a data.frame, matrix, or a file path (csv, xlsx).")
  }

  # Data Loading and ID Standardization
  raw <- load_any(data_path)
  raw_ids <- as.character(raw[[1]])
  raw_ids <- stringi::stri_unescape_unicode(raw_ids)

  if (auto_clean) {
    ids <- clean_labels(raw_ids)

    # Handle cases where standardization fails
    if (any(is.na(ids))) {
      na_count <- sum(is.na(ids))
      cli::cli_alert_info("Note: {na_count} names could not be standardized. Keeping original labels.")
      ids[is.na(ids)] <- raw_ids[is.na(ids)]
    }
  } else {
    # Omics bypass: keeps exact names from the file
    ids <- raw_ids
  }

  # Ensure all feature labels are unique
  if (any(duplicated(ids))) {
    cli::cli_alert_warning("Duplicated labels detected! Appending suffixes to keep them unique.")
    ids <- make.unique(ids)
  }

  # Data Matrix Cleaning
  clean_mat <- apply(raw[, -1, drop = FALSE], 2, function(x) {
    if (is.character(x)) {
      x <- gsub("[ \t\r\n]", "", gsub(",", ".", x))
    }
    val <- as.numeric(x)
    val[val == 0] <- NA
    return(val)
  })

  rownames(clean_mat) <- ids

  # Metadata Loading and Sample Alignment
  meta <- load_any(meta_path)

  meta_ids <- as.character(meta$sample)
  data_ids <- colnames(clean_mat)

  # Compare sample IDs between data and metadata
  ex_meta <- setdiff(meta_ids, data_ids)
  ex_data <- setdiff(data_ids, meta_ids)
  common <- intersect(meta_ids, data_ids)

  # Report alignment discrepancies
  if (length(ex_meta) > 0) {
    cli::cli_alert_warning("Removed from Metadata: {length(ex_meta)} samples (No data match):")
    cli::cli_text("{.field {ex_meta}}")
  }

  if (length(ex_data) > 0) {
    cli::cli_alert_warning("Removed from Data: {length(ex_data)} samples (No metadata match):")
    cli::cli_text("{.field {ex_data}}")
  }

  # Final Data Harmonization
  meta <- meta[match(common, meta$sample), ]
  meta$group <- factor(meta$group)
  data <- clean_mat[, meta$sample, drop = FALSE]

  # Global Quality Control - 100% NA samples removing
  empty_samples <- colSums(!is.na(data)) == 0

  if (any(empty_samples)) {
    n_empty <- sum(empty_samples)
    names_empty <- colnames(data)[empty_samples]

    cli::cli_alert_danger("Removed {n_empty} sample{?s} with 100% missing values:")
    cli::cli_text("{.field {names_empty}}")

    data <- data[, !empty_samples, drop = FALSE]
    meta <- meta[!empty_samples, , drop = FALSE]
  }

  n_oxy <- nrow(data)
  n_samples <- ncol(data)
  group_counts <- table(meta$group)

  # Console Output Summary
  cli::cli_h1("staRoxy Data Summary")
  cli::cli_alert_success("Successful upload of {n_oxy} features across {n_samples} samples.")

  cli::cli_h2("Experimental Groups:")
  for (g_name in names(group_counts)) {
    cli::cli_li("{g_name}: n = {group_counts[g_name]}")
  }
  cli::cli_rule()

  # Object Construction
  obj <- list(
    data = data,
    meta = meta,
    info = list(
      n_features = n_oxy,
      n_samples = n_samples,
      groups = names(group_counts),
      date_created = Sys.time()
    )
  )

  class(obj) <- "staRoxy"
  return(invisible(obj))
}
