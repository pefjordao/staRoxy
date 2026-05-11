#' @title Modify the Oxylipin Reference Universe
#'
#' @description
#' Allows the user to easily customize the reference universe for Over-representation
#' Analysis (ORA) by removing specific oxylipins or adding new ones alongside
#' their biological precursor class.
#'
#' @details
#' \strong{Nomenclature & Unicode Compatibility:}
#' To ensure perfect string matching with the internal \code{staRoxy} database,
#' always use Unicode escape sequences for Greek letters when adding new features
#' or precursors. Common codes include:
#' \itemize{
#'   \item \strong{Alpha (\eqn{\alpha}):} \code{\\u03b1} (e.g., \code{"PGF2\\u03b1"})
#'   \item \strong{Beta (\eqn{\beta}):} \code{\\u03b2}
#'   \item \strong{Gamma (\eqn{\gamma}):} \code{\\u03b3}
#'   \item \strong{Delta (\eqn{\Delta}):} \code{\\u0394}
#' }
#'
#' @param univ A data frame containing the base universe. Defaults to the
#' output of \code{get_oxy_universe()}.
#' @param remove Character vector. Names of the oxylipins to completely remove
#' from the universe (deletes the entire row). Default is \code{NULL}.
#' @param add_name Character vector. Names of the new oxylipins to add.
#' Default is \code{NULL}.
#' @param add_precursor Character vector. The precursor classes corresponding
#' to \code{add_name}. Must be the same length as \code{add_name}. See Details
#' for Unicode formatting. Default is \code{NULL}.
#'
#' @return A customized data frame representing the new oxylipin universe.
#'
#' @examples
#' # Get the default universe
#' my_univ <- get_oxy_universe()
#'
#' # Remove specific lipids (e.g., removing PGE2 and PGD2)
#' my_univ <- edit_oxy_universe(my_univ, remove = c("PGE2", "PGD2"))
#'
#' # Add a new lipid and precursor using Standard ASCII
#' my_univ <- edit_oxy_universe(
#'   univ = my_univ,
#'   add_name = "New-HETE",
#'   add_precursor = "ARA (20:4 n-6)"
#' )
#'
#' # Add a new lipid requiring a Greek letter (Alpha)
#' my_univ <- edit_oxy_universe(
#'   univ = my_univ,
#'   add_name = "New-PGF2\u03b1",
#'   add_precursor = "ARA (20:4 n-6)"
#' )
#'
#' @importFrom dplyr filter bind_rows distinct
#' @importFrom cli cli_alert_success cli_alert_danger cli_alert_warning
#'
#' @export
edit_oxy_universe <- function(univ = get_oxy_universe(),
                              remove = NULL,
                              add_name = NULL,
                              add_precursor = NULL) {

  # Removal Phase
  if (!is.null(remove)) {
    initial_n <- nrow(univ)
    univ <- dplyr::filter(univ, !(Oxylipin %in% remove))
    removed_n <- initial_n - nrow(univ)

    if (removed_n > 0) {
      cli::cli_alert_success("Removed {removed_n} oxylipin{?s} from the universe.")
    } else {
      cli::cli_alert_warning("None of the specified oxylipins to remove were found in the universe.")
    }
  }

  # Addition Phase
  if (!is.null(add_name)) {

    # Safety Check: Must have matching precursors
    if (is.null(add_precursor) || length(add_name) != length(add_precursor)) {
      cli::cli_alert_danger("Error: 'add_name' and 'add_precursor' must be provided together and have the exact same length.")
      return(invisible(univ))
    }

    # Create the new dataframe block
    new_data <- data.frame(
      Oxylipin = as.character(add_name),
      Precursor = as.character(add_precursor),
      stringsAsFactors = FALSE
    )

    # Bind and ensure uniqueness (in case user tries to add something already there)
    univ <- dplyr::bind_rows(univ, new_data)
    univ <- dplyr::distinct(univ, Oxylipin, .keep_all = TRUE)

    cli::cli_alert_success("Added {length(add_name)} new oxylipin{?s}. Current universe size: {nrow(univ)}.")
  }

  return(univ)
}
