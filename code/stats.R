# Statistics Export Script
#
# This script gathers statistical results from analysis modules and saves them as PDFs
# under assets/stats.

# nolint start
library('tidyverse')

source(file.path('code', 'stomata_aperture.R'))
source(file.path('code', 'blended_peels_starch_staining.R'))
source(file.path('code', 'live_dead_staining.R'))

assets_dir <- file.path('assets')
stats_dir <- file.path(assets_dir, 'stats')
dir.create(stats_dir, showWarnings = FALSE, recursive = TRUE)
run_stats_dir <- file.path(stats_dir, 'runs')
dir.create(run_stats_dir, showWarnings = FALSE, recursive = TRUE)

get_significance_symbol <- function(p_value) {
    vapply(p_value, function(value) {
        if (is.na(value)) return('NA')
        if (value < 0.001) return('***')
        if (value < 0.01) return('**')
        if (value < 0.05) return('*')
        'ns'
    }, character(1))
}

format_p_values <- function(stats_table) {
    if ('p_value' %in% names(stats_table)) {
        stats_table$p_value <- format.pval(stats_table$p_value, digits = 3, eps = 0.001)
    }
    stats_table
}

flatten_pairwise_results <- function(pairwise_results) {
    rows <- list()
    if (length(pairwise_results) == 0) return(tibble::tibble())

    for (group_name in names(pairwise_results)) {
        group_results <- pairwise_results[[group_name]]
        for (comp_name in names(group_results)) {
            comparison <- group_results[[comp_name]]
            rows[[length(rows) + 1]] <- data.frame(
                group = group_name,
                comparison = comp_name,
                timepoint_1 = comparison$timepoint_1,
                timepoint_2 = comparison$timepoint_2,
                n_1 = comparison$n_1,
                n_2 = comparison$n_2,
                statistic = comparison$statistic,
                p_value = comparison$p_value,
                significance = get_significance_symbol(comparison$p_value),
                stringsAsFactors = FALSE
            )
        }
    }

    dplyr::bind_rows(rows)
}

save_stats_pdf <- function(stats_table, file_path, title = NULL) {
    stats_table <- format_p_values(stats_table)
    stats_table <- tibble::as_tibble(stats_table)
    lines <- capture.output({
        if (!is.null(title)) {
            cat(title, "\n")
            cat(strrep("-", nchar(title)), "\n\n")
        }
        if (nrow(stats_table) == 0) {
            cat("(no rows)\n")
        } else {
            print(stats_table, n = nrow(stats_table), width = Inf)
        }
    })

    pdf(file_path, width = 8.5, height = 11)
    plot.new()
    text(0, 1, paste(lines, collapse = "\n"), adj = c(0, 1), family = 'mono', cex = 0.7)
    dev.off()
}

# --------------------------------------------------------------------------------------------
# STOMATA APERTURE
# --------------------------------------------------------------------------------------------

arabidopsis_data <- load_aperture_data(arabidopsis_aperture_file)
arabidopsis_ratios <- compute_ratios(arabidopsis_data, type = 'arabidopsis')
arabidopsis_width <- extract_width_data(arabidopsis_data, type = 'arabidopsis')

arabidopsis_area_stats <- pairwise_comparisons(arabidopsis_ratios,
                                               type = 'arabidopsis', metric = 'area')
arabidopsis_width_stats <- pairwise_comparisons(arabidopsis_width,
                                                type = 'arabidopsis', metric = 'width')

arabidopsis_area_table <- flatten_pairwise_results(arabidopsis_area_stats)
arabidopsis_width_table <- flatten_pairwise_results(arabidopsis_width_stats)

save_stats_pdf(arabidopsis_area_table,
               file.path(stats_dir, 'stomata_aperture_arabidopsis_area_stats.pdf'),
               title = 'Arabidopsis pooled: area ratios')
save_stats_pdf(arabidopsis_width_table,
               file.path(stats_dir, 'stomata_aperture_arabidopsis_width_stats.pdf'),
               title = 'Arabidopsis pooled: outer ledge width')

arabidopsis_area_run_stats <- pairwise_comparisons(arabidopsis_ratios,
                                                   type = 'arabidopsis_run', metric = 'area')
arabidopsis_width_run_stats <- pairwise_comparisons(arabidopsis_width,
                                                    type = 'arabidopsis_run', metric = 'width')
arabidopsis_area_run_table <- flatten_pairwise_results(arabidopsis_area_run_stats)
arabidopsis_width_run_table <- flatten_pairwise_results(arabidopsis_width_run_stats)

save_stats_pdf(arabidopsis_area_run_table,
               file.path(run_stats_dir, 'stomata_aperture_arabidopsis_area_run_stats.pdf'),
               title = 'Arabidopsis by run: area ratios')
save_stats_pdf(arabidopsis_width_run_table,
               file.path(run_stats_dir, 'stomata_aperture_arabidopsis_width_run_stats.pdf'),
               title = 'Arabidopsis by run: outer ledge width')

species_data <- load_aperture_data_species(species_aperture_file)
species_ratios <- compute_ratios(species_data, type = 'species')
species_width <- extract_width_data(species_data, type = 'species')

species_area_stats <- pairwise_comparisons(species_ratios,
                                           type = 'species', metric = 'area')
species_width_stats <- pairwise_comparisons(species_width,
                                            type = 'species', metric = 'width')

species_area_table <- flatten_pairwise_results(species_area_stats)
species_width_table <- flatten_pairwise_results(species_width_stats)

save_stats_pdf(species_area_table,
               file.path(stats_dir, 'stomata_aperture_species_area_stats.pdf'),
               title = 'Species pooled: area ratios')
save_stats_pdf(species_width_table,
               file.path(stats_dir, 'stomata_aperture_species_width_stats.pdf'),
               title = 'Species pooled: outer ledge width')

species_area_run_stats <- pairwise_comparisons(species_ratios,
                                               type = 'species_run', metric = 'area')
species_width_run_stats <- pairwise_comparisons(species_width,
                                                type = 'species_run', metric = 'width')
species_area_run_table <- flatten_pairwise_results(species_area_run_stats)
species_width_run_table <- flatten_pairwise_results(species_width_run_stats)

save_stats_pdf(species_area_run_table,
               file.path(run_stats_dir, 'stomata_aperture_species_area_run_stats.pdf'),
               title = 'Species by run: area ratios')
save_stats_pdf(species_width_run_table,
               file.path(run_stats_dir, 'stomata_aperture_species_width_run_stats.pdf'),
               title = 'Species by run: outer ledge width')

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

save_stats_pdf(blended_stats,
               file.path(stats_dir, 'blended_peels_starch_staining_stats.pdf'),
               title = 'Blended peels starch staining')

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

save_stats_pdf(live_dead_stats,
               file.path(stats_dir, 'live_dead_staining_stats.pdf'),
               title = 'Live dead staining')

# nolint end
