# Combined Plots Script
#
# Builds composite figures from analysis modules.
# Figure 1: blended peels starch staining + stomata aperture area.

# nolint start
library('tidyverse')

source(file.path('code', 'stomata_aperture.R'))
source(file.path('code', 'blended_peels_starch_staining.R'))
source(file.path('code', 'live_dead_staining.R'))

assets_dir <- file.path('assets')
plots_dir <- file.path(assets_dir, 'plots')
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# Use a larger canvas than A4 so dense multi-panel figures keep readable spacing.
figure_page_width_in <- 14
figure_page_height_in <- 18

open_pdf_device <- function(file_path, width, height) {
    tryCatch({
        grDevices::pdf(file_path, width = width, height = height)
        file_path
    }, error = function(e) {
        stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        fallback <- file.path(
            dirname(file_path),
            paste0(tools::file_path_sans_ext(basename(file_path)), "_", stamp, ".pdf")
        )
        grDevices::pdf(fallback, width = width, height = height)
        message("Could not write ", file_path, " (maybe open). Writing to ", fallback)
        fallback
    })
}

# Compute a single, shared y-range across related plots so panels are comparable.
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
# DATA SOURCES
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
blended_image_paths <- build_blended_example_images()

aperture_data <- load_aperture_data(arabidopsis_aperture_file)
aperture_ratios <- compute_ratios(aperture_data, type = 'arabidopsis')
aperture_width <- extract_width_data(aperture_data, type = 'arabidopsis')

blended_y_limits <- compute_pretty_y_limits(blended_data$starch_guard_cell_area_ratio,
                                            lower = 0, step = 5, pad_frac = 0.2)
stomata_area_y_limits <- compute_pretty_y_limits(aperture_ratios$area_ratio_percent,
                                                 lower = 0, step = 5, pad_frac = 0.2)
stomata_width_y_limits <- compute_pretty_y_limits(aperture_width$outer_ledge_width,
                                                  lower = 0, step = 1, pad_frac = 0.2)

build_arabidopsis_image_paths <- function(base_dir) {
    list(
        'SD' = list(
            'NaOH Buffer' = list(
                'Mock' = list(
                    'EoN' = file.path(base_dir, '8h_daloso_mock_eon.png'),
                    '1hBL' = file.path(base_dir, '8h_daloso_mock_1hbl.png')
                ),
                'Mannitol' = list(
                    'EoN' = file.path(base_dir, '8h_daloso_mannitol_eon.png'),
                    '1hBL' = file.path(base_dir, '8h_daloso_mannitol_1hbl.png')
                )
            ),
            'MES-BTP Buffer' = list(
                'Mock' = list(
                    'EoN' = file.path(base_dir, '8h_santelia_mock_eon.png'),
                    '1hBL' = file.path(base_dir, '8h_santelia_mock_1hbl.png')
                ),
                'Mannitol' = list(
                    'EoN' = file.path(base_dir, '8h_santelia_mannitol_eon.png'),
                    '1hBL' = file.path(base_dir, '8h_santelia_mannitol_1hbl.png')
                )
            )
        ),
        'ND' = list(
            'NaOH Buffer' = list(
                'Mock' = list(
                    'EoN' = file.path(base_dir, '12h_daloso_mock_eon.png'),
                    '1hBL' = file.path(base_dir, '12h_daloso_mock_1hbl.png')
                ),
                'Mannitol' = list(
                    'EoN' = file.path(base_dir, '12h_daloso_mannitol_eon.png'),
                    '1hBL' = file.path(base_dir, '12h_daloso_mannitol_1hbl.png')
                )
            ),
            'MES-BTP Buffer' = list(
                'Mock' = list(
                    'EoN' = file.path(base_dir, '12h_santelia_mock_eon.png'),
                    '1hBL' = file.path(base_dir, '12h_santelia_mock_1hbl.png')
                ),
                'Mannitol' = list(
                    'EoN' = file.path(base_dir, '12h_santelia_mannitol_eon.png'),
                    '1hBL' = file.path(base_dir, '12h_santelia_mannitol_1hbl.png')
                )
            )
        )
    )
}

arabidopsis_image_paths <- build_arabidopsis_image_paths(
    file.path('data', 'Stomata aperture', 'examples')
)

arabidopsis_width_image_paths <- build_arabidopsis_image_paths(
    file.path('data', 'Stomata aperture', 'examples', 'width_sreen_shots', 'cropped')
)

# --------------------------------------------------------------------------------------------
# HELPERS
# --------------------------------------------------------------------------------------------

build_image_paths_subset <- function(image_paths, photoperiod_value, buffer_value) {
    buffer_paths <- image_paths[[photoperiod_value]][[buffer_value]]
    if (is.null(buffer_paths)) return(NULL)
    setNames(list(setNames(list(buffer_paths), buffer_value)), photoperiod_value)
}

build_image_paths_subset_by_treatment <- function(image_paths, treatment_value, buffer_value) {
    subset_paths <- list()
    for (photoperiod_value in c('ND', 'SD')) {
        treatment_paths <- image_paths[[photoperiod_value]][[buffer_value]][[treatment_value]]
        if (is.null(treatment_paths)) next
        subset_paths[[photoperiod_value]] <- setNames(
            list(setNames(list(treatment_paths), treatment_value)),
            buffer_value
        )
    }
    if (length(subset_paths) == 0) return(NULL)
    subset_paths
}

build_blended_section_plot <- function(data, buffer_value, treatment_value,
                                       image_paths = NULL, y_limits = NULL,
                                       photoperiod_letters = NULL) {
    section_data <- data %>%
        dplyr::filter(buffer == buffer_value, condition == treatment_value) %>%
        dplyr::mutate(
            buffer = factor(buffer, levels = buffer_value),
            photoperiod = factor(photoperiod, levels = c('ND', 'SD')),
            condition = factor(condition, levels = treatment_value)
        )

    section_stats <- compare_timepoints_by_group(section_data)

    section_images <- NULL
    if (!is.null(image_paths)) {
        section_images <- build_image_paths_subset_by_treatment(image_paths, treatment_value, buffer_value)
    }

    plot_obj <- plot_blended_peels_starch(section_data,
                                          stats = section_stats,
                                          y_limits = y_limits,
                                          facet_photoperiod = TRUE,
                                          facet_condition = FALSE,
                                          facet_buffer = FALSE,
                                          image_paths = section_images,
                                          panel_letters = photoperiod_letters,
                                          strip_text_face = 'plain')
    plot_obj
}

build_stomata_section_plot <- function(data, buffer_value, treatment_value,
                                       image_paths, y_limits = NULL,
                                       photoperiod_letters = NULL) {
    section_data <- data %>%
        dplyr::filter(Buffer == buffer_value, Treatment == treatment_value) %>%
        dplyr::mutate(
            Buffer = factor(Buffer, levels = buffer_value),
            Photoperiod = factor(Photoperiod, levels = c('ND', 'SD')),
            Treatment = factor(Treatment, levels = treatment_value)
        )

    section_stats <- pairwise_comparisons(section_data, type = 'arabidopsis', metric = 'area')
    section_images <- build_image_paths_subset_by_treatment(image_paths, treatment_value, buffer_value)

    plot_obj <- plot_violin(section_data, section_stats, type = 'arabidopsis', metric = 'area',
                            image_paths = section_images,
                            facet_buffer = FALSE,
                            facet_photoperiod = TRUE,
                            facet_treatment = FALSE,
                            y_limits = y_limits,
                            panel_letters = photoperiod_letters,
                            strip_text_face = 'plain')
    plot_obj
}

build_stomata_width_section_plot <- function(data, buffer_value, treatment_value,
                                             image_paths, y_limits = NULL,
                                             photoperiod_letters = NULL) {
    section_data <- data %>%
        dplyr::filter(Buffer == buffer_value, Treatment == treatment_value) %>%
        dplyr::mutate(
            Buffer = factor(Buffer, levels = buffer_value),
            Photoperiod = factor(Photoperiod, levels = c('ND', 'SD')),
            Treatment = factor(Treatment, levels = treatment_value)
        )

    section_stats <- pairwise_comparisons(section_data, type = 'arabidopsis', metric = 'width')
    section_images <- build_image_paths_subset_by_treatment(image_paths, treatment_value, buffer_value)

    plot_obj <- plot_violin(section_data, section_stats, type = 'arabidopsis', metric = 'width',
                            image_paths = section_images,
                            facet_buffer = FALSE,
                            facet_photoperiod = TRUE,
                            facet_treatment = FALSE,
                            y_limits = y_limits,
                            panel_letters = photoperiod_letters,
                            strip_text_face = 'plain')
    plot_obj
}

draw_section <- function(section, plots) {
    grid::pushViewport(grid::viewport(
        layout.pos.row = section$row,
        layout.pos.col = section$col
    ))

    section_layout <- grid::grid.layout(
        nrow = 2, ncol = 1,
        heights = grid::unit(c(0.04, 0.96), 'npc')
    )
    grid::pushViewport(grid::viewport(layout = section_layout))

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::grid.text(
        label = section$label,
        x = grid::unit(0.5, 'npc'),
        y = grid::unit(0.35, 'npc'),
        just = c('center', 'center'),
        gp = grid::gpar(fontface = 'bold', fontsize = 20)
    )
    grid::popViewport()

    panel_layout <- grid::grid.layout(
        nrow = 3, ncol = 1,
        heights = grid::unit(c(0.49, 0.02, 0.49), 'npc')
    )

    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 0.96, height = 0.97))
    grid::pushViewport(grid::viewport(layout = panel_layout))

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 0.98, height = 0.98))
    grid::grid.draw(plots$blended)
    grid::popViewport(2)

    grid::pushViewport(grid::viewport(layout.pos.row = 3, layout.pos.col = 1))
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 0.98, height = 0.98))
    grid::grid.draw(plots$stomata)
    grid::popViewport(2)

    grid::popViewport(2)
    grid::popViewport()

    grid::grid.rect(gp = grid::gpar(col = 'black', fill = NA, lwd = 1))
    grid::popViewport()
    grid::popViewport()
}

# --------------------------------------------------------------------------------------------
# FIGURE 1
# --------------------------------------------------------------------------------------------

sections <- list(
    list(buffer = 'MES-BTP Buffer', treatment = 'Mock',
         row = 1, col = 1, label = 'MES-BTP Buffer - Mock'),
    list(buffer = 'NaOH Buffer', treatment = 'Mock',
         row = 1, col = 2, label = 'NaOH Buffer - Mock'),
    list(buffer = 'MES-BTP Buffer', treatment = 'Mannitol',
         row = 2, col = 1, label = 'MES-BTP Buffer - Mannitol'),
    list(buffer = 'NaOH Buffer', treatment = 'Mannitol',
         row = 2, col = 2, label = 'NaOH Buffer - Mannitol')
)

panel_letters <- letters[1:16]  # 4 sections × 4 plots each = 16 letters
width_panel_letters <- letters[1:8]  # 4 sections x 2 plots each = 8 letters
for (i in seq_along(sections)) {
    section_row <- sections[[i]]$row
    section_col <- sections[[i]]$col

    # Figure 1 letter order by display matrix:
    # row 1: a c e g
    # row 2: b d f h
    # row 3: i k m o
    # row 4: j l n p
    base_index <- (section_row - 1) * 8 + (section_col - 1) * 4
    sections[[i]]$photoperiod_letters <- c(
        'ND' = panel_letters[base_index + 1],
        'SD' = panel_letters[base_index + 3]
    )
    sections[[i]]$stomata_photoperiod_letters <- c(
        'ND' = panel_letters[base_index + 2],
        'SD' = panel_letters[base_index + 4]
    )

    # Supplementary width letter order:
    # row 1: a c e g
    # row 2: b d f h
    width_nd_index <- (section_col - 1) * 4 + ifelse(section_row == 1, 1, 2)
    sections[[i]]$stomata_width_letters <- c(
        'ND' = width_panel_letters[width_nd_index],
        'SD' = width_panel_letters[width_nd_index + 2]
    )
}

as_grob <- function(plot_obj) {
    if (inherits(plot_obj, "patchwork")) {
        return(patchwork::patchworkGrob(plot_obj))
    }
    ggplot2::ggplotGrob(plot_obj)
}

build_section_plots <- function() {
    lapply(sections, function(section) {
        blended_plot <- build_blended_section_plot(
            blended_data, section$buffer, section$treatment, blended_image_paths,
            y_limits = blended_y_limits,
            photoperiod_letters = section$photoperiod_letters
        )

        stomata_plot <- build_stomata_section_plot(
            aperture_ratios, section$buffer, section$treatment, arabidopsis_image_paths,
            y_limits = stomata_area_y_limits,
            photoperiod_letters = section$stomata_photoperiod_letters
        )

        list(blended = as_grob(blended_plot), stomata = as_grob(stomata_plot))
    })
}

build_width_section_plots <- function() {
    lapply(sections, function(section) {
        width_plot <- build_stomata_width_section_plot(
            aperture_width, section$buffer, section$treatment, arabidopsis_width_image_paths,
            y_limits = stomata_width_y_limits,
            photoperiod_letters = section$stomata_width_letters
        )
        as_grob(width_plot)
    })
}

figure_width <- figure_page_width_in
figure_height <- figure_page_height_in

draw_figure_1 <- function(section_plots) {
    grid::grid.newpage()

    grid::pushViewport(grid::viewport(width = 0.96, height = 0.96))
    outer_layout <- grid::grid.layout(
        nrow = 2, ncol = 2,
        widths = grid::unit(c(1, 1), 'null'),
        heights = grid::unit(c(1, 1), 'null')
    )
    grid::pushViewport(grid::viewport(layout = outer_layout))

    for (i in seq_along(sections)) {
        draw_section(sections[[i]], section_plots[[i]])
    }

    grid::popViewport(2)
}

section_plots <- build_section_plots()
figure_1_path <- open_pdf_device(file.path(plots_dir, 'figure_1.pdf'),
                                 width = figure_width, height = figure_height)
draw_figure_1(section_plots)
grDevices::dev.off()

# --------------------------------------------------------------------------------------------
# SUPPLEMENTARY: STOMATA WIDTH
# --------------------------------------------------------------------------------------------

draw_stomata_width_section <- function(section, plot_grob) {
    grid::pushViewport(grid::viewport(
        layout.pos.row = section$row,
        layout.pos.col = section$col
    ))

    section_layout <- grid::grid.layout(
        nrow = 2, ncol = 1,
        heights = grid::unit(c(0.08, 0.92), 'npc')
    )
    grid::pushViewport(grid::viewport(layout = section_layout))

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::grid.text(
        label = section$label,
        x = grid::unit(0.5, 'npc'),
        y = grid::unit(0.5, 'npc'),
        just = c('center', 'center'),
        gp = grid::gpar(fontface = 'bold', fontsize = 20)
    )
    grid::popViewport()

    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 0.96, height = 0.95))
    grid::grid.draw(plot_grob)
    grid::popViewport(2)

    grid::grid.rect(gp = grid::gpar(col = 'black', fill = NA, lwd = 1))
    grid::popViewport(2)
}

figure_s1_width <- figure_page_width_in
figure_s1_height <- figure_page_height_in

draw_figure_s1_stomata_width <- function(section_plots) {
    grid::grid.newpage()

    grid::pushViewport(grid::viewport(width = 0.96, height = 0.96))
    outer_layout <- grid::grid.layout(
        nrow = 2, ncol = 2,
        widths = grid::unit(c(1, 1), 'null'),
        heights = grid::unit(c(1, 1), 'null')
    )
    grid::pushViewport(grid::viewport(layout = outer_layout))

    for (i in seq_along(sections)) {
        draw_stomata_width_section(sections[[i]], section_plots[[i]])
    }

    grid::popViewport(2)
}

width_section_plots <- build_width_section_plots()
figure_s1_path <- open_pdf_device(file.path(plots_dir, 'figure_S1_stomata_width.pdf'),
                                  width = figure_s1_width, height = figure_s1_height)
draw_figure_s1_stomata_width(width_section_plots)
grDevices::dev.off()

# --------------------------------------------------------------------------------------------
# FIGURE 2
# --------------------------------------------------------------------------------------------

live_dead_file <- file.path('data', 'Live dead staining', 'DalosoRebuttal_LiveDead_Overall_Data.xls')
live_dead_data <- load_live_dead_data(live_dead_file, sheet = 'Daloso_Rebuttal_Overall_Data')

live_dead_filtered <- filter_live_dead_data(
    live_dead_data,
    buffer_levels = c('MES-BTP Buffer', 'NaOH Buffer'),
    additive_levels = c('Mock', 'Mannitol'),
    growing_levels = c('ND', 'SD'),
    blending_levels = c('30s+30s')
)

live_dead_example_paths <- build_live_dead_example_images()
live_dead_y_limits <- c(0, 125)

build_live_dead_section_plot <- function(data, buffer_value, photoperiod_value,
                                         image_paths = NULL, y_limits = NULL,
                                         panel_label = NULL) {
    section_data <- data %>%
        dplyr::filter(buffer == buffer_value, growing_condition == photoperiod_value) %>%
        dplyr::mutate(
            buffer = factor(buffer, levels = buffer_value),
            growing_condition = factor(growing_condition, levels = photoperiod_value),
            additives = factor(additives, levels = c('Mock', 'Mannitol'))
        )

    section_stats <- compare_additives_by_group(section_data)

    section_images <- NULL
    if (!is.null(image_paths)) {
        section_images <- build_image_paths_subset(image_paths, photoperiod_value, buffer_value)
    }

    plot_live_dead_viability(
        section_data,
        stats = section_stats,
        y_limits = y_limits,
        image_paths = section_images,
        facet_buffer = FALSE,
        facet_growing_condition = FALSE,
        panel_label = panel_label
    )
}

live_dead_sections <- list(
    list(buffer = 'MES-BTP Buffer', photoperiod = 'ND',
         row = 1, col = 1, label = '(a)', title = 'MES-BTP Buffer - ND'),
    list(buffer = 'NaOH Buffer', photoperiod = 'ND',
         row = 1, col = 2, label = '(c)', title = 'NaOH Buffer - ND'),
    list(buffer = 'MES-BTP Buffer', photoperiod = 'SD',
         row = 2, col = 1, label = '(b)', title = 'MES-BTP Buffer - SD'),
    list(buffer = 'NaOH Buffer', photoperiod = 'SD',
         row = 2, col = 2, label = '(d)', title = 'NaOH Buffer - SD')
)

live_dead_plots <- lapply(live_dead_sections, function(section) {
    plot_obj <- build_live_dead_section_plot(
        live_dead_filtered,
        section$buffer,
        section$photoperiod,
        image_paths = live_dead_example_paths,
        y_limits = live_dead_y_limits,
        panel_label = section$label
    )
    as_grob(plot_obj)
})

# Load and create leaf starch plot
source(file.path('code', 'leaf_starch.R'))
leaf_starch_result <- analyze_leaf_starch(y_limits = c(0, 12), panel_label = '(e)')
leaf_starch_grob <- as_grob(leaf_starch_result$plot)

draw_live_dead_section <- function(section, plot_grob) {
    grid::pushViewport(grid::viewport(
        layout.pos.row = section$row,
        layout.pos.col = section$col
    ))

    section_layout <- grid::grid.layout(
        nrow = 2, ncol = 1,
        heights = grid::unit(c(0.08, 0.92), 'npc')
    )
    grid::pushViewport(grid::viewport(layout = section_layout))

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::grid.text(
        label = section$title,
        x = grid::unit(0.5, 'npc'),
        y = grid::unit(0.5, 'npc'),
        just = c('center', 'center'),
        gp = grid::gpar(fontface = 'bold', fontsize = 16)
    )
    grid::popViewport()

    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 0.96, height = 0.95))
    grid::grid.draw(plot_grob)
    grid::popViewport(2)

    grid::grid.rect(gp = grid::gpar(col = 'black', fill = NA, lwd = 1))
    grid::popViewport(2)
}

figure_2_width <- figure_page_width_in
figure_2_height <- figure_page_height_in

draw_figure_2 <- function() {
    grid::grid.newpage()

    grid::pushViewport(grid::viewport(width = 0.92, height = 0.92))
    outer_layout <- grid::grid.layout(
        nrow = 3, ncol = 2,
        widths = grid::unit(c(1, 1), 'null'),
        heights = grid::unit(c(1, 1, 1), 'null')
    )
    grid::pushViewport(grid::viewport(layout = outer_layout))

    # Draw the 4 live/dead sections
    for (i in seq_along(live_dead_sections)) {
        draw_live_dead_section(live_dead_sections[[i]], live_dead_plots[[i]])
    }

    # Draw leaf starch in row 3, spanning both columns
    grid::pushViewport(grid::viewport(
        layout.pos.row = 3,
        layout.pos.col = 1:2
    ))

    leaf_layout <- grid::grid.layout(
        nrow = 2, ncol = 1,
        heights = grid::unit(c(0.08, 0.92), 'npc')
    )
    grid::pushViewport(grid::viewport(layout = leaf_layout))

    grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid::grid.text(
        label = 'Leaf Starch Content',
        x = grid::unit(0.5, 'npc'),
        y = grid::unit(0.5, 'npc'),
        just = c('center', 'center'),
        gp = grid::gpar(fontface = 'bold', fontsize = 16)
    )
    grid::popViewport()

    grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
    grid::pushViewport(grid::viewport(x = 0.5, y = 0.5, width = 0.5, height = 0.95))
    grid::grid.draw(leaf_starch_grob)
    grid::popViewport(2)

    grid::grid.rect(gp = grid::gpar(col = 'black', fill = NA, lwd = 1))
    grid::popViewport(2)

    grid::popViewport(2)
}

figure_2_path <- open_pdf_device(file.path(plots_dir, 'figure_2.pdf'),
                                 width = figure_2_width, height = figure_2_height)
draw_figure_2()
grDevices::dev.off()

# nolint end
