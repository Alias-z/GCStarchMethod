# Plot Export Script
#
# This script gathers plots from the analysis modules and saves them as PDFs
# under assets/plots.

# nolint start
library('tidyverse')

source(file.path('code', 'stomata_aperture.R'))
source(file.path('code', 'blended_peels_starch_staining.R'))
source(file.path('code', 'live_dead_staining.R'))

assets_dir <- file.path('assets')
plots_dir <- file.path(assets_dir, 'plots')
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
run_plots_dir <- file.path(plots_dir, 'runs')
dir.create(run_plots_dir, showWarnings = FALSE, recursive = TRUE)
species_plots_dir <- file.path(plots_dir, 'species')
dir.create(species_plots_dir, showWarnings = FALSE, recursive = TRUE)
species_run_plots_dir <- file.path(species_plots_dir, 'runs')
dir.create(species_run_plots_dir, showWarnings = FALSE, recursive = TRUE)

golden_ratio <- 1.618

compute_pretty_y_limits <- function(values, lower = 0, step = 5, pad_frac = 0.2) {
    values <- values[is.finite(values)]
    if (length(values) == 0) return(c(lower, lower + step))

    max_val <- max(values, na.rm = TRUE)
    min_val <- min(values, na.rm = TRUE)
    range_val <- max_val - min_val
    upper_raw <- if (is.finite(range_val) && range_val > 0) {
        max_val + range_val * pad_frac
    } else {
        max_val + abs(max_val) * pad_frac + step
    }

    upper <- ceiling(upper_raw / step) * step
    if (!is.finite(upper) || upper <= lower) upper <- lower + step
    c(lower, upper)
}

# --------------------------------------------------------------------------------------------
# STOMATA APERTURE
# --------------------------------------------------------------------------------------------

arabidopsis_image_paths <- list(
    'SD' = list(
        'NaOH Buffer' = list(
            'Mock' = list(
                'EoN' = file.path('data', 'Stomata aperture', 'examples',
                                 '8h_daloso_mock_eon.png'),
                '1hBL' = file.path('data', 'Stomata aperture', 'examples',
                                  '8h_daloso_mock_1hbl.png')
            ),
            'Mannitol' = list(
                'EoN' = file.path('data', 'Stomata aperture', 'examples',
                                 '8h_daloso_mannitol_eon.png'),
                '1hBL' = file.path('data', 'Stomata aperture', 'examples',
                                  '8h_daloso_mannitol_1hbl.png')
            )
        ),
        'MES-BTP Buffer' = list(
            'Mock' = list(
                'EoN' = file.path('data', 'Stomata aperture', 'examples',
                                 '8h_santelia_mock_eon.png'),
                '1hBL' = file.path('data', 'Stomata aperture', 'examples',
                                  '8h_santelia_mock_1hbl.png')
            ),
            'Mannitol' = list(
                'EoN' = file.path('data', 'Stomata aperture', 'examples',
                                 '8h_santelia_mannitol_eon.png'),
                '1hBL' = file.path('data', 'Stomata aperture', 'examples',
                                  '8h_santelia_mannitol_1hbl.png')
            )
        )
    ),
    'ND' = list(
        'NaOH Buffer' = list(
            'Mock' = list(
                'EoN' = file.path('data', 'Stomata aperture', 'examples',
                                 '12h_daloso_mock_eon.png'),
                '1hBL' = file.path('data', 'Stomata aperture', 'examples',
                                  '12h_daloso_mock_1hbl.png')
            ),
            'Mannitol' = list(
                'EoN' = file.path('data', 'Stomata aperture', 'examples',
                                 '12h_daloso_mannitol_eon.png'),
                '1hBL' = file.path('data', 'Stomata aperture', 'examples',
                                  '12h_daloso_mannitol_1hbl.png')
            )
        ),
        'MES-BTP Buffer' = list(
            'Mock' = list(
                'EoN' = file.path('data', 'Stomata aperture', 'examples',
                                 '12h_santelia_mock_eon.png'),
                '1hBL' = file.path('data', 'Stomata aperture', 'examples',
                                  '12h_santelia_mock_1hbl.png')
            ),
            'Mannitol' = list(
                'EoN' = file.path('data', 'Stomata aperture', 'examples',
                                 '12h_santelia_mannitol_eon.png'),
                '1hBL' = file.path('data', 'Stomata aperture', 'examples',
                                  '12h_santelia_mannitol_1hbl.png')
            )
        )
    )
)

species_image_paths <- list(
    'V. unguiculata' = list(
        'EoN' = file.path('data', 'Stomata aperture', 'examples', 'cowpea_eon.png'),
        '1hBL' = file.path('data', 'Stomata aperture', 'examples', 'cowpea_1hbl.png')
    ),
    'N. tabacum' = list(
        'EoN' = file.path('data', 'Stomata aperture', 'examples', 'tobacco_eon.png'),
        '1hBL' = file.path('data', 'Stomata aperture', 'examples', 'tobacco_1hbl.png')
    )
)

arabidopsis_data_for_limits <- load_aperture_data(arabidopsis_aperture_file)
arabidopsis_ratios_for_limits <- compute_ratios(arabidopsis_data_for_limits, type = 'arabidopsis')
arabidopsis_width_for_limits <- extract_width_data(arabidopsis_data_for_limits, type = 'arabidopsis')

species_data_for_limits <- load_aperture_data_species(species_aperture_file)
species_ratios_for_limits <- compute_ratios(species_data_for_limits, type = 'species')
species_width_for_limits <- extract_width_data(species_data_for_limits, type = 'species')

stomata_area_y_limits <- compute_pretty_y_limits(
    c(arabidopsis_ratios_for_limits$area_ratio_percent,
      species_ratios_for_limits$area_ratio_percent),
    lower = 0, step = 5, pad_frac = 0.2
)
stomata_width_y_limits <- compute_pretty_y_limits(
    c(arabidopsis_width_for_limits$outer_ledge_width,
      species_width_for_limits$outer_ledge_width),
    lower = 0, step = 1, pad_frac = 0.2
)

arabidopsis_pooled <- analyze_arabidopsis_pooled(
    image_paths = arabidopsis_image_paths,
    area_y_limits = stomata_area_y_limits,
    width_y_limits = stomata_width_y_limits
)
ggsave(file.path(plots_dir, 'stomata_aperture_arabidopsis_area.pdf'),
       arabidopsis_pooled$area_plot,
       device = 'pdf', height = 15, width = 10 * golden_ratio)
ggsave(file.path(plots_dir, 'stomata_aperture_arabidopsis_width.pdf'),
       arabidopsis_pooled$width_plot,
       device = 'pdf', height = 15, width = 10 * golden_ratio)

arabidopsis_run <- analyze_arabidopsis_by_run(
    area_y_limits = stomata_area_y_limits,
    width_y_limits = stomata_width_y_limits
)
ggsave(file.path(run_plots_dir, 'stomata_aperture_arabidopsis_area_run.pdf'),
       arabidopsis_run$area_plot,
       device = 'pdf', height = 15, width = 10 * golden_ratio)
ggsave(file.path(run_plots_dir, 'stomata_aperture_arabidopsis_width_run.pdf'),
       arabidopsis_run$width_plot,
       device = 'pdf', height = 15, width = 10 * golden_ratio)

species_pooled <- analyze_species_pooled(
    image_paths = species_image_paths,
    area_y_limits = stomata_area_y_limits,
    width_y_limits = stomata_width_y_limits
)
ggsave(file.path(species_plots_dir, 'stomata_aperture_species_area.pdf'),
       species_pooled$area_plot,
       device = 'pdf', height = 8, width = 8 * golden_ratio)
ggsave(file.path(species_plots_dir, 'stomata_aperture_species_width.pdf'),
       species_pooled$width_plot,
       device = 'pdf', height = 8, width = 8 * golden_ratio)

species_run <- analyze_species_by_run(
    area_y_limits = stomata_area_y_limits,
    width_y_limits = stomata_width_y_limits
)
ggsave(file.path(species_run_plots_dir, 'stomata_aperture_species_area_run.pdf'),
       species_run$area_plot,
       device = 'pdf', height = 8, width = 8 * golden_ratio)
ggsave(file.path(species_run_plots_dir, 'stomata_aperture_species_width_run.pdf'),
       species_run$width_plot,
       device = 'pdf', height = 8, width = 8 * golden_ratio)

# --------------------------------------------------------------------------------------------
# BLENDED PEELS STARCH STAINING
# --------------------------------------------------------------------------------------------

blended_12 <- load_blended_peels_starch_data(
    file.path('data', 'BL blended peels starch staining', '12-12-summary', '1212-summary.csv'),
    photoperiod = 'ND'
)
blended_8 <- load_blended_peels_starch_data(
    file.path('data', 'BL blended peels starch staining', '8-16-summary', 'Daloso-816-summary.csv'),
    photoperiod = 'SD'
)
blended_data <- dplyr::bind_rows(blended_12, blended_8)
blended_stats <- compare_timepoints_by_group(blended_data)
blended_image_paths <- build_blended_example_images()
blended_y_limits <- compute_pretty_y_limits(blended_data$starch_guard_cell_area_ratio,
                                            lower = 0, step = 5, pad_frac = 0.2)
blended_plot <- plot_blended_peels_starch(blended_data,
                                          stats = blended_stats,
                                          y_limits = blended_y_limits,
                                          facet_photoperiod = TRUE,
                                          image_paths = blended_image_paths)
ggsave(file.path(plots_dir, 'blended_peels_starch_staining.pdf'),
       blended_plot, device = 'pdf', height = 10, width = 12)

# --------------------------------------------------------------------------------------------
# LIVE DEAD STAINING
# --------------------------------------------------------------------------------------------

live_dead_data <- load_live_dead_data(
    file.path('data', 'Live dead staining', 'DalosoRebuttal_LiveDead_Overall_Data.xls'),
    sheet = 'Daloso_Rebuttal_Overall_Data'
)
live_dead_filtered <- filter_live_dead_data(
    live_dead_data,
    buffer_levels = c('MES-BTP Buffer', 'NaOH Buffer'),
    additive_levels = c('Mock', 'Mannitol'),
    growing_levels = c('ND', 'SD'),
    blending_levels = c('30s+30s')
)
live_dead_stats <- compare_additives_by_group(live_dead_filtered)
live_dead_image_paths <- build_live_dead_example_images()
live_dead_plot <- plot_live_dead_viability(
    live_dead_filtered,
    stats = live_dead_stats,
    y_limits = c(0, 110),
    image_paths = live_dead_image_paths,
    facet_growing_condition = TRUE
)
ggsave(file.path(plots_dir, 'live_dead_staining.pdf'),
       live_dead_plot, device = 'pdf', height = 10, width = 12)

# --------------------------------------------------------------------------------------------
# LEAF STARCH
# --------------------------------------------------------------------------------------------

source(file.path('code', 'leaf_starch.R'))

leaf_starch_result <- analyze_leaf_starch(y_limits = c(0, 10))
ggsave(file.path(plots_dir, 'leaf_starch.pdf'),
       leaf_starch_result$plot, device = 'pdf', width = 6, height = 6)

# nolint end
