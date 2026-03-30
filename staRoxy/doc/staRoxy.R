## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  message = TRUE,
  warning = FALSE,
  comment = "#>",
  fig.width = 7,
  fig.height = 6
)

## -----------------------------------------------------------------------------
library(staRoxy)

## -----------------------------------------------------------------------------
data(data_oxy_pellet)
data(metadata_oxy_pellet)

## -----------------------------------------------------------------------------
staRoxy_object <- read_oxy(data_path = data_oxy_pellet,
                           meta_path = metadata_oxy_pellet)

## -----------------------------------------------------------------------------
staRoxy_object <- filter_oxy(staRoxy_object,
                             type = "I",
                             prop = 0.5)

## -----------------------------------------------------------------------------
staRoxy_object <- transform_oxy(staRoxy_object)

## -----------------------------------------------------------------------------
check_data_distribution(staRoxy_object)

## -----------------------------------------------------------------------------
plot_na_diagnostic(staRoxy_object,
                   mode = "cor")

## -----------------------------------------------------------------------------
plot_na_diagnostic(staRoxy_object,
                   mode = "dist")

## -----------------------------------------------------------------------------
get_exclusives_oxy(staRoxy_object)

## -----------------------------------------------------------------------------
get_analyzed_oxy(staRoxy_object)

## -----------------------------------------------------------------------------
descriptives_oxy(staRoxy_object,
                 sort_by = "mean")

## -----------------------------------------------------------------------------
descriptives_oxy(staRoxy_object,
                 group = "Trypsin")

## -----------------------------------------------------------------------------
# Define group-specific colors for visualization
my_colors <- c("Trypsin" = "#F65F72", "Scraper" = "#492376")

plot_pca_oxy(staRoxy_object,
             colors = my_colors)

## -----------------------------------------------------------------------------
plot_heatmap_oxy(staRoxy_object,
                 colors = my_colors)

## -----------------------------------------------------------------------------
plot_heatmap_oxy(staRoxy_object,
                 colors = my_colors,
                 hide_na = TRUE)

## -----------------------------------------------------------------------------
plot_heatmap_oxy(obj = staRoxy_object,
                 colors = my_colors,
                 z_cap = 2)

## -----------------------------------------------------------------------------
profile_test(staRoxy_object,
             colors = my_colors,
             plot = TRUE)

## -----------------------------------------------------------------------------
stats <- limma_oxy(staRoxy_object)

## -----------------------------------------------------------------------------
list_contrasts(stats)

## -----------------------------------------------------------------------------
plot_rank_oxy(stats,
              contrast = 1,
              top_n = 15,
              rank_by = "logFC")

## -----------------------------------------------------------------------------
plot_rank_oxy(stats = stats,
              contrast = 1,
              top_n = 15,
              rank_by = "p")

## -----------------------------------------------------------------------------
plot_violin_oxy(staRoxy_object,
                stats,
                oxylipin = "TXB3",
                colors = my_colors,
                show_sig = TRUE,
                sig_y_nudge = 0.4)

## -----------------------------------------------------------------------------
plot_violin_oxy(staRoxy_object,
                stats,
                oxylipin = "TXB3",
                colors = my_colors,
                show_sig = TRUE,
                sig_y_nudge = 0.6,
                log_scale = FALSE,
                unit = "nmol/mg")

## -----------------------------------------------------------------------------
cor_within_oxy(staRoxy_object,
               oxylipin1 = "TXB3",
               oxylipin2 = "12-HHTrE",
               group = "Trypsin",
               colors = my_colors)

## -----------------------------------------------------------------------------
cor_between_oxy(staRoxy_object,
                oxylipin = "TXB3",
                group1 = "Scraper",
                group2 = "Trypsin")

## -----------------------------------------------------------------------------
loading_dissimilarity(staRoxy_object,
                      group1 = "Scraper",
                      group2 = "Trypsin",
                      pc_index = 1)

## -----------------------------------------------------------------------------
classify_oxy(staRoxy_object,
             return_mode = "summary")

## -----------------------------------------------------------------------------
classify_oxy(staRoxy_object,
             return_mode = "species")

## -----------------------------------------------------------------------------
oxy_ora(staRoxy_object)

## ----echo=FALSE---------------------------------------------------------------
sessionInfo()

