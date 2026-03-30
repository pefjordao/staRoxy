#' Oxylipidomics Abundance Data (Cellular Pellets)
#'
#' A dataset containing oxylipidomics abundance levels from RAW 264.7 murine
#' macrophage pellets. This dataset compares two harvesting methods:
#' mechanical scraping (Scraper, n = 6) and enzymatic treatment (Trypsin, n = 6).
#'
#' @format A data frame with lipids as rows and samples as columns.
#' @source Internal experimental data from RAW 264.7 cell line studies.
"data_oxy_pellet"

#' Metadata for Cellular Pellet Samples
#'
#' Mapping of sample identifiers to experimental groups (Scraper vs. Trypsin)
#' for the cellular pellet dataset.
#'
#' @format A data frame with 12 rows and 2 columns:
#' \describe{
#'   \item{sample}{Unique sample identifier}
#'   \item{group}{Harvesting method used (Scraper or Trypsin)}
#' }
"metadata_oxy_pellet"

#' Oxylipidomics Abundance Data (Cell Culture Supernatant)
#'
#' A dataset of secreted oxylipins in RAW 264.7 culture media. It evaluates
#' the impact of harvesting procedures across three conditions: Scraper (n = 6),
#' Trypsin (n = 6), and Control (no manipulation, n = 6).
#'
#' @format A data frame with lipids as rows and samples as columns.
#' @source Secreted inflammatory mediator profiles in macrophage supernatants.
"data_oxy_supernatant"

#' Metadata for Cell Culture Supernatant Samples
#'
#' Mapping of sample identifiers to experimental groups (Scraper, Trypsin,
#' and Control) for the culture supernatant dataset.
#'
#' @format A data frame with 18 rows and 2 columns:
#' \describe{
#'   \item{sample}{Unique sample identifier}
#'   \item{group}{Experimental condition (Scraper, Trypsin, or Control)}
#' }
"metadata_oxy_supernatant"

#' Oxylipidomics Abundance Data (Blood Plasma - ALS Mouse Model)
#'
#' A comprehensive oxylipidomics dataset from blood plasma of SOD1-G93A mice
#' (a model for Amyotrophic Lateral Sclerosis) and age-matched wild-type (WT)
#' controls. Samples were collected at two biological stages: asymptomatic
#' (70 days old) and symptomatic (120 days old).
#'
#' @format A data frame with lipids as rows and 24 samples as columns.
#' @source Chaves-Filho, A. B., et al. (2023). Free Radical Biology and Medicine.
#' \doi{10.1016/j.freeradbiomed.2023.08.019}
"data_oxy_plasma"

#' Metadata for Blood Plasma Samples (ALS Model)
#'
#' Mapping of 24 sample identifiers to their respective experimental groups
#' within the ALS mouse model study.
#'
#' @format A data frame with 24 rows and 2 columns:
#' \describe{
#'   \item{sample}{Unique sample identifier matching data headers}
#'   \item{group}{Experimental group: ALS_70d, WT_70d, ALS_120d, or WT_120d}
#' }
#' @source Chaves-Filho, A. B., et al. (2023). \doi{10.1016/j.freeradbiomed.2023.08.019}
"metadata_oxy_plasma"
