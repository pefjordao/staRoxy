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
#'   and \code{.txt} (via \code{readxl} and \code{data.table}).
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
#' @param data_path Character. Path to the oxylipin abundance file. The first
#' column must contain oxylipin identifiers.
#' @param meta_path Character. Path to the metadata file. It must contain
#' at least a \code{sample} column (matching the data headers) and a
#' \code{group} column.
#'
#' @return An object of class \code{staRoxy}, which is a list containing:
#' \itemize{
#'   \item \strong{data:} A numeric matrix of oxylipin abundances.
#'   \item \strong{meta:} A data frame of experimental metadata.
#'   \item \strong{info:} A list of summary statistics and creation date.
#' }
#'
#' @importFrom tools file_ext
#' @importFrom readxl read_excel
#' @importFrom data.table fread
#' @importFrom cli cli_h1 cli_alert_success cli_alert_warning cli_alert_info cli_h2 cli_li cli_rule
#'
#' @examples
#' # Example using paths to internal or external files
#' # staRoxy_object <- read_oxy(data_path = "data_oxylipins.csv",
#' #                           meta_path = "metadata.csv")
#'
#' @export
read_oxy <- function(data_path, meta_path) {

  # Internal Helper Functions
  # Define a helper to load different file formats (Excel vs. Text/CSV)
  load_any <- function(path) {
    if (is.data.frame(path)) return(as.data.frame(path))

    if (is.character(path)) {
      ext <- tolower(tools::file_ext(path))
      if (ext %in% c("xls", "xlsx")) {
        return(as.data.frame(readxl::read_excel(path)))
      }
      return(data.table::fread(path, data.table = FALSE, check.names = FALSE, encoding = "UTF-8"))
    }

    stop("Input must be a data.frame or a file path (csv, xlsx).")
  }

  # Data Loading and ID Standardization
  raw <- load_any(data_path)
  raw_ids <- as.character(raw[[1]])
  ids <- clean_labels(raw_ids)

  # Handle cases where standardization fails
  if (any(is.na(ids))) {
    na_count <- sum(is.na(ids))
    cli::cli_alert_info("Note: {na_count} names could not be standardized. Keeping original labels.")
    ids[is.na(ids)] <- raw_ids[is.na(ids)]
  }

  # Ensure all feature labels are unique
  if (any(duplicated(ids))) {
    cli::cli_alert_warning("Duplicated labels detected! Appending suffixes to keep them unique.")
    ids <- make.unique(ids)
  }

  # Data Matrix Cleaning
  # Convert character strings to numeric, handle separators, and treat zeros as NAs
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
  # Reorder and filter both datasets to match common samples
  meta <- meta[match(common, meta$sample), ]
  meta$group <- factor(meta$group)
  data <- clean_mat[, meta$sample, drop = FALSE]

  n_oxy <- nrow(data)
  n_samples <- ncol(data)
  group_counts <- table(meta$group)

  # Console Output Summary
  cli::cli_h1("staRoxy Data Summary")
  cli::cli_alert_success("Successful upload of {n_oxy} oxylipins across {n_samples} samples.")

  cli::cli_h2("Experimental Groups:")
  for (g_name in names(group_counts)) {
    cli::cli_li("{g_name}: n = {group_counts[g_name]}")
  }
  cli::cli_rule()

  # Object Construction
  # Assemble the final staRoxy object
  obj <- list(
    data = data,
    meta = meta,
    info = list(
      n_oxylipins = n_oxy,
      n_samples = n_samples,
      groups = names(group_counts),
      date_created = Sys.time()
    )
  )

  class(obj) <- "staRoxy"
  return(invisible(obj))
}
