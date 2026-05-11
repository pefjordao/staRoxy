#' Blood Plasma Oxylipin Data (ALS Mouse Model)
#'
#' A comprehensive oxylipidomics dataset from blood plasma of SOD1-G93A mice
#' (a model for Amyotrophic Lateral Sclerosis) and age-matched wild-type (WT)
#' controls. Samples were collected at two biological stages: asymptomatic
#' (70 days old) and symptomatic (120 days old).
#'
#' @format A data frame with oxylipins as rows and 24 samples as columns.
#' @source Chaves-Filho, A. B., et al. (2023). Free Radical Biology and Medicine.
#' \doi{10.1016/j.freeradbiomed.2023.08.019}
"data_oxy_plasma"

#' Metadata for Blood Plasma Samples
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

#' RAW 264.7 Macrophage Pellet Oxylipin Data (LPS Stimulated)
#'
#' Intracellular oxylipin abundance profiles from RAW 264.7 murine
#' macrophages. This dataset compares basal inflammatory states (Control)
#' against activation with Lipopolysaccharide (LPS, 100 ng/mL).
#'
#' @format A data frame with oxylipin species as rows and 8 samples as columns.
#' Abundances are normalized to total protein content (fmol/mg protein).
#' @source Internal experimental dataset following the standardized
#' oxylipidomics protocol.
"data_oxy_lps_pellet"

#' Metadata for RAW 264.7 Pellet Samples
#'
#' Experimental metadata mapping the 8 cellular pellet samples to their
#' respective groups: "Control" or "LPS".
#'
#' @format A data frame with 8 rows and 2 columns:
#' \describe{
#'   \item{sample}{Unique sample identifier matching data headers}
#'   \item{group}{Experimental group (Control or LPS)}
#' }
"metadata_oxy_lps_pellet"

#' RAW 264.7 Macrophage Supernatant Oxylipin Data (LPS Stimulated)
#'
#' Secreted oxylipin signatures (exometabolome) from RAW 264.7 murine
#' macrophages. It captures the flux of lipid mediators released into the
#' extracellular environment following TLR4 activation.
#'
#' @format A data frame with oxylipin species as rows and 7 samples as columns.
#' Abundances are normalized by sample volume (fmol/mL).
#' @source Internal experimental dataset following the standardized
#' oxylipidomics protocol.
"data_oxy_lps_supernatant"

#' Metadata for RAW 264.7 Supernatant Samples
#'
#' Experimental metadata mapping the 8 culture supernatant samples to their
#' respective groups: "Control" or "LPS".
#'
#' @format A data frame with 7 rows and 2 columns:
#' \describe{
#'   \item{sample}{Unique sample identifier matching data headers}
#'   \item{group}{Experimental group (Control or LPS)}
#' }
"metadata_oxy_lps_supernatant"

#' Default Oxylipin Reference Universe
#'
#' A pre-defined analytical universe mapping 125 oxylipins to their
#' respective fatty acid precursors (e.g., ARA, LA, DHA). This dataset
#' serves as the default background for Over-representation Analysis (ORA).
#'
#' @format A data frame with 125 rows and 2 columns:
#' \describe{
#'   \item{Oxylipin}{The standardized identifier used within the \code{staRoxy} pipeline.}
#'   \item{Precursor}{The biological parent fatty acid precursor class.}
#' }
#' @source Internal biochemical mapping based on literature and
#' commercial standards (Cayman Chemical/LIPID MAPS).
"oxy_universe_default"

# Criar data do LPS
