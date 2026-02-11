# Live Dead Staining Analysis Module
#
# This module provides functions to analyze guard cell viability data and create
# reproducible plots for live/dead staining experiments.
#
# ============================================================================================
# USAGE EXAMPLES
# ============================================================================================
#
# Load data:
#   data <- load_live_dead_data(
#       file.path("data", "Live dead staining", "LiveDead_R.xlsx"),
#       sheet = "R"
#   )
#
# Filter and plot:
#   data_filtered <- filter_live_dead_data(
#       data,
#       buffer_levels = c("MES-BTP Buffer", "NaOH Buffer"),
#       additive_levels = c("Mock", "Mannitol"),
#       growing_levels = c("ND", "SD")
#   )
#
#   stats <- compare_additives_by_group(data_filtered)
#   plot <- plot_live_dead_viability(data_filtered, stats = stats, facet_growing_condition = TRUE)
#
# nolint start
library('tidyverse')
library('readxl')
library('ggimage')
library('patchwork')

additive_colors <- c('Mock' = 'gray', 'Mannitol' = 'gray')
additive_labels <- c('Mock' = 'Mock', 'Mannitol' = 'Mannitol')

daloso_plot_theme <- function() {
    theme_bw() +
        theme(
            panel.background = element_blank(),
            panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
            panel.spacing = unit(1.2, 'lines'),
            legend.position = 'none',
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text = element_text(size = 16),
            strip.text = element_text(size = 16, face = 'italic'),
            strip.background = element_rect(fill = 'grey98'),
            plot.title = element_blank(),
            plot.margin = margin(t = 10, r = 20, b = 40, l = 10, unit = 'pt')
        )
}

# ============================================================================================
# HELPERS
# ============================================================================================

normalize_names <- function(names) {
    names <- iconv(names, to = "ASCII//TRANSLIT")
    names[is.na(names)] <- ""
    names <- tolower(names)
    names <- gsub("[^a-z0-9]+", "_", names)
    names <- gsub("^_|_$", "", names)
    names
}

get_significance_symbol <- function(p_value) {
    vapply(p_value, function(value) {
        if (is.na(value)) return('NA')
        if (value < 0.001) return('***')
        if (value < 0.01) return('**')
        if (value < 0.05) return('*')
        'ns'
    }, character(1))
}

normalize_photoperiod_value <- function(value) {
    raw <- tolower(trimws(as.character(value)))
    raw[is.na(raw)] <- ""
    digits <- gsub("\\D", "", raw)
    vapply(seq_along(raw), function(i) {
        item <- raw[i]
        if (item %in% c("sd", "shortday", "short")) return("SD")
        if (item %in% c("nd", "neutralday", "neutral")) return("ND")
        if (digits[i] %in% c("8", "816")) return("SD")
        if (digits[i] %in% c("12", "1212")) return("ND")
        if (item == "") return(NA_character_)
        value[i]
    }, character(1))
}

normalize_buffer_value <- function(value) {
    raw <- tolower(trimws(as.character(value)))
    raw[is.na(raw)] <- ""
    vapply(seq_along(raw), function(i) {
        item <- raw[i]
        if (item %in% c("santelia buffer", "santelia", "mes-btp", "mes btp", "mesbtp", "mes-btp buffer")) {
            return("MES-BTP Buffer")
        }
        if (item %in% c("daloso buffer", "daloso", "naoh", "naoh buffer")) {
            return("NaOH Buffer")
        }
        if (item == "") return(NA_character_)
        value[i]
    }, character(1))
}

normalize_additive_value <- function(value) {
    raw <- tolower(trimws(as.character(value)))
    raw[is.na(raw)] <- ""
    vapply(seq_along(raw), function(i) {
        item <- raw[i]
        if (item %in% c("null", "mock", "control")) return("Mock")
        if (item %in% c("mannitol")) return("Mannitol")
        if (item == "") return(NA_character_)
        value[i]
    }, character(1))
}

safe_wilcox_test <- function(values, groups) {
    if (length(unique(groups)) < 2) return(NA_real_)
    tryCatch(
        stats::wilcox.test(values ~ groups, exact = FALSE)$p.value,
        error = function(e) NA_real_
    )
}

safe_left_join <- function(x, y, by) {
    if (length(by) == 0) {
        return(dplyr::cross_join(x, y))
    }
    dplyr::left_join(x, y, by = by)
}

filter_violin_data <- function(data, response_col, panel_vars, x_var) {
    data <- data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(panel_vars, x_var)))) %>%
        dplyr::filter(dplyr::n() > 1, dplyr::n_distinct(.data[[response_col]]) > 1) %>%
        dplyr::ungroup()

    if (length(panel_vars) == 0) {
        if (dplyr::n_distinct(data[[x_var]]) < 2) {
            return(data[0, ])
        }
        return(data)
    }

    data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(panel_vars))) %>%
        dplyr::filter(dplyr::n_distinct(.data[[x_var]]) > 1) %>%
        dplyr::ungroup()
}

# ============================================================================================
# DATA LOADING
# ============================================================================================

build_live_dead_example_images <- function(example_dir = file.path('data',
                                                                  'Live dead staining',
                                                                  'examples')) {
    if (!dir.exists(example_dir)) return(NULL)

    files <- list.files(
        example_dir,
        pattern = "\\.(png|jpg|jpeg|tif|tiff)$",
        full.names = TRUE,
        ignore.case = TRUE
    )
    if (length(files) == 0) return(NULL)

    image_paths <- list()

    for (path in files) {
        stem <- tolower(tools::file_path_sans_ext(basename(path)))
        stem <- gsub("[^a-z0-9]+", "_", stem)
        stem <- gsub("^_|_$", "", stem)
        parts <- strsplit(stem, "_")[[1]]
        if (length(parts) < 2) next

        photoperiod <- NA_character_
        if (parts[1] %in% c("nd", "sd")) {
            photoperiod <- toupper(parts[1])
        } else if ("nd" %in% parts) {
            photoperiod <- "ND"
        } else if ("sd" %in% parts) {
            photoperiod <- "SD"
        }

        buffer <- NA_character_
        if ("naoh" %in% parts || "daloso" %in% parts) {
            buffer <- "NaOH Buffer"
        }
        if ("mes" %in% parts || "santelia" %in% parts || ("mes" %in% parts && "btp" %in% parts)) {
            buffer <- "MES-BTP Buffer"
        }

        additives <- NA_character_
        if ("mock" %in% parts || "null" %in% parts) {
            additives <- "Mock"
        }
        if ("mannitol" %in% parts) {
            additives <- "Mannitol"
        }

        if (any(is.na(c(photoperiod, buffer, additives)))) next

        if (is.null(image_paths[[photoperiod]])) image_paths[[photoperiod]] <- list()
        if (is.null(image_paths[[photoperiod]][[buffer]])) image_paths[[photoperiod]][[buffer]] <- list()
        image_paths[[photoperiod]][[buffer]][[additives]] <- path
    }

    if (length(image_paths) == 0) return(NULL)
    image_paths
}

load_live_dead_data <- function(file_path, sheet = "R", growing_condition_value = NULL) {
    data <- readxl::read_excel(file_path, sheet = sheet)
    names(data) <- normalize_names(names(data))

    if (!"growing_condition" %in% names(data) && "growth_condition" %in% names(data)) {
        data <- data %>% dplyr::rename(growing_condition = growth_condition)
    }

    required_cols <- c("buffer", "additives", "alive_gc")
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    has_growing <- "growing_condition" %in% names(data)
    if (!has_growing && is.null(growing_condition_value)) {
        stop("Missing required columns: growing_condition")
    }

    data %>%
        dplyr::mutate(
            buffer = normalize_buffer_value(buffer),
            additives = normalize_additive_value(additives),
            growing_condition = if (!is.null(growing_condition_value)) {
                normalize_photoperiod_value(growing_condition_value)
            } else if (has_growing) {
                normalize_photoperiod_value(growing_condition)
            } else {
                NA_character_
            },
            blending_time = if ("blending_time" %in% names(data)) as.character(blending_time) else NA_character_,
            alive_gc = as.numeric(alive_gc),
            dead_gc = if ("dead_gc" %in% names(data)) as.numeric(dead_gc) else NA_real_
        )
}

filter_live_dead_data <- function(data,
                                  buffer_levels = NULL,
                                  additive_levels = NULL,
                                  growing_levels = NULL,
                                  blending_levels = NULL) {
    filtered <- data

    if (!is.null(buffer_levels)) {
        filtered <- filtered %>% dplyr::filter(buffer %in% buffer_levels)
        filtered <- filtered %>% dplyr::mutate(buffer = factor(buffer, levels = buffer_levels))
    } else {
        filtered <- filtered %>% dplyr::mutate(buffer = factor(buffer))
    }

    if (!is.null(additive_levels)) {
        filtered <- filtered %>% dplyr::filter(additives %in% additive_levels)
        filtered <- filtered %>% dplyr::mutate(additives = factor(additives, levels = additive_levels))
    } else {
        filtered <- filtered %>% dplyr::mutate(additives = factor(additives))
    }

    if (!is.null(growing_levels)) {
        filtered <- filtered %>% dplyr::filter(growing_condition %in% growing_levels)
        filtered <- filtered %>% dplyr::mutate(growing_condition = factor(growing_condition,
                                                                          levels = growing_levels))
    } else if ("growing_condition" %in% names(filtered)) {
        filtered <- filtered %>% dplyr::mutate(growing_condition = factor(growing_condition))
    }

    if (!is.null(blending_levels) && "blending_time" %in% names(filtered)) {
        filtered <- filtered %>% dplyr::filter(blending_time %in% blending_levels)
        filtered <- filtered %>% dplyr::mutate(blending_time = factor(blending_time,
                                                                      levels = blending_levels))
    } else if ("blending_time" %in% names(filtered)) {
        filtered <- filtered %>% dplyr::mutate(blending_time = factor(blending_time))
    }

    filtered
}

# ============================================================================================
# STATISTICS
# ============================================================================================

compare_additives_by_group <- function(data,
                                       response_col = "alive_gc",
                                       group_vars = c("buffer", "growing_condition")) {
    group_vars <- group_vars[group_vars %in% names(data)]

    data %>%
        dplyr::filter(is.finite(.data[[response_col]]), !is.na(additives)) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
        dplyr::summarise(
            p_value = safe_wilcox_test(.data[[response_col]], additives),
            counts = list(table(additives)),
            .groups = 'drop'
        ) %>%
        dplyr::mutate(
            n_mock = vapply(counts, function(x) {
                if ("Mock" %in% names(x)) unname(x[["Mock"]]) else 0L
            }, integer(1)),
            n_mannitol = vapply(counts, function(x) {
                if ("Mannitol" %in% names(x)) unname(x[["Mannitol"]]) else 0L
            }, integer(1)),
            significance = get_significance_symbol(p_value)
        ) %>%
        dplyr::select(-counts)
}

compare_buffers_by_group <- function(data,
                                     response_col = "alive_gc",
                                     group_vars = c("additives", "growing_condition")) {
    group_vars <- group_vars[group_vars %in% names(data)]

    data %>%
        dplyr::filter(!is.na(.data[[response_col]]), !is.na(buffer)) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
        dplyr::summarise(
            p_value = safe_wilcox_test(.data[[response_col]], buffer),
            counts = list(table(buffer)),
            .groups = 'drop'
        ) %>%
        dplyr::mutate(
            n_mes_btp_buffer = vapply(counts, function(x) {
                if ("MES-BTP Buffer" %in% names(x)) unname(x[["MES-BTP Buffer"]]) else 0L
            }, integer(1)),
            n_naoh_buffer = vapply(counts, function(x) {
                if ("NaOH Buffer" %in% names(x)) unname(x[["NaOH Buffer"]]) else 0L
            }, integer(1)),
            significance = get_significance_symbol(p_value)
        ) %>%
        dplyr::select(-counts)
}

prepare_buffer_labels <- function(data,
                                  stats,
                                  response_col = "alive_gc",
                                  group_vars = c("additives", "growing_condition")) {
    if (is.null(stats) || nrow(stats) == 0) return(NULL)
    group_vars <- group_vars[group_vars %in% names(stats)]

    y_max <- data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
        dplyr::summarise(label_y = max(.data[[response_col]], na.rm = TRUE) * 1.05,
                         .groups = 'drop')

    safe_left_join(stats, y_max, by = group_vars)
}

prepare_additive_labels <- function(data,
                                    stats,
                                    response_col = "alive_gc",
                                    group_vars = c("buffer", "growing_condition")) {
    if (is.null(stats) || nrow(stats) == 0) return(NULL)
    group_vars <- group_vars[group_vars %in% names(stats)]

    y_max <- data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
        dplyr::summarise(label_y = max(.data[[response_col]], na.rm = TRUE) * 1.1,
                         .groups = 'drop')

    safe_left_join(stats, y_max, by = group_vars)
}

# ============================================================================================
# PLOTTING
# ============================================================================================

build_live_dead_image_strip <- function(plot_data, image_paths, facet_formula, facet_vars) {
    if (is.null(image_paths)) return(NULL)
    if (is.null(plot_data) || nrow(plot_data) == 0) return(NULL)
    if (!"additives" %in% names(plot_data)) return(NULL)

    additive_levels <- levels(plot_data$additives)
    if (is.null(additive_levels) || length(additive_levels) == 0) return(NULL)

    panel_map <- plot_data %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(facet_vars))) %>%
        tidyr::crossing(additives = factor(additive_levels, levels = additive_levels))

    image_rows <- list()
    for (growing in names(image_paths)) {
        for (buffer in names(image_paths[[growing]])) {
            for (additive in names(image_paths[[growing]][[buffer]])) {
                path <- image_paths[[growing]][[buffer]][[additive]]
                if (is.null(path)) next

                row <- data.frame(
                    growing_condition = growing,
                    buffer = buffer,
                    additives = additive,
                    image_path = path,
                    stringsAsFactors = FALSE
                )
                image_rows[[length(image_rows) + 1]] <- row
            }
        }
    }

    image_data <- dplyr::bind_rows(image_rows)
    if (nrow(image_data) == 0) return(NULL)

    for (key in c("buffer", "growing_condition", "additives")) {
        if (key %in% names(image_data) && key %in% names(plot_data)) {
            allowed <- unique(as.character(plot_data[[key]]))
            image_data <- image_data %>% dplyr::filter(.data[[key]] %in% allowed)
        }
    }

    if (nrow(image_data) == 0) return(NULL)

    if ("buffer" %in% names(panel_map) && "buffer" %in% names(plot_data)) {
        image_data$buffer <- factor(image_data$buffer, levels = levels(plot_data$buffer))
    }
    if ("growing_condition" %in% names(panel_map) && "growing_condition" %in% names(plot_data)) {
        image_data$growing_condition <- factor(image_data$growing_condition,
                                               levels = levels(plot_data$growing_condition))
    }
    image_data$additives <- factor(image_data$additives, levels = additive_levels)

    p <- ggplot(panel_map, aes(x = additives, y = 0.5)) +
        geom_blank() +
        ggimage::geom_image(
            data = image_data,
            aes(image = image_path),
            inherit.aes = TRUE,
            size = 0.92
        ) +
        scale_x_discrete(labels = additive_labels) +
        scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
        theme_void() +
        theme(
            axis.text.x = element_text(size = 14, color = 'gray30'),
            axis.ticks.x = element_blank(),
            strip.text = element_blank(),
            strip.background = element_blank(),
            panel.background = element_blank(),
            panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
            panel.spacing = unit(1.2, 'lines'),
            plot.margin = margin(t = 0, r = 20, b = 10, l = 10, unit = 'pt')
        )

    if (!is.null(facet_formula) && length(facet_vars) > 0) {
        p <- p + facet_grid(facet_formula, drop = FALSE, as.table = TRUE,
                            labeller = labeller(.cols = function(x) x))
    }

    p
}

plot_live_dead_viability <- function(data,
                                     stats = NULL,
                                     response_col = "alive_gc",
                                     y_limits = NULL,
                                     image_paths = NULL,
                                     facet_buffer = TRUE,
                                     facet_growing_condition = TRUE,
                                     plot_tag = NULL,
                                     panel_label = NULL) {

    plot_data <- data %>%
        dplyr::filter(is.finite(.data[[response_col]]), !is.na(buffer), !is.na(additives))

    facet_vars <- character(0)
    if (facet_buffer && "buffer" %in% names(plot_data)) {
        facet_vars <- c(facet_vars, 'buffer')
    }
    if (facet_growing_condition && "growing_condition" %in% names(plot_data)) {
        facet_vars <- c(facet_vars, 'growing_condition')
    }

    facet_formula <- NULL
    if (length(facet_vars) > 0) {
        facet_formula <- if (facet_buffer && facet_growing_condition &&
            "buffer" %in% names(plot_data) &&
            "growing_condition" %in% names(plot_data)) {
            buffer ~ growing_condition
        } else if (facet_buffer && "buffer" %in% names(plot_data)) {
            buffer ~ .
        } else if (facet_growing_condition && "growing_condition" %in% names(plot_data)) {
            . ~ growing_condition
        } else {
            NULL
        }
    }

    y_max_data <- max(plot_data[[response_col]], na.rm = TRUE)
    y_min_data <- min(plot_data[[response_col]], na.rm = TRUE)
    y_range <- y_max_data - y_min_data
    extra_top <- ifelse(is.finite(y_range) && y_range > 0,
                        y_range * 0.35,
                        abs(y_max_data) * 0.2 + 1)

    y_lower <- y_min_data * 0.95
    y_upper <- y_max_data + extra_top
    y_upper <- max(y_upper, y_max_data + 12)

    if (!is.null(y_limits) && length(y_limits) == 2 && all(is.finite(y_limits))) {
        y_lower <- y_limits[1]
        y_upper <- y_limits[2]
    }
    extended_y_max <- y_upper

    # Set y-axis breaks for display (0-100 only)
    y_breaks <- seq(0, 100, by = 20)

    violin_data <- filter_violin_data(plot_data, response_col, facet_vars, 'additives')

    p <- ggplot(plot_data, aes(x = additives, y = .data[[response_col]], fill = additives))
    if (nrow(violin_data) > 0) {
        p <- p + geom_violin(data = violin_data, color = 'black', width = 0.9)
    }
    p <- p +
        geom_jitter(color = 'black', width = 0.05, size = 0.5, alpha = 0.7) +
        stat_summary(fun = mean, geom = 'point',
                     shape = 23, size = 3, fill = 'white', color = 'white') +
        scale_fill_manual(values = additive_colors, labels = additive_labels) +
        scale_x_discrete(labels = additive_labels) +
        scale_y_continuous(breaks = y_breaks) +
        labs(x = NULL, y = 'alive guard cells (%) \n') +
        daloso_plot_theme() +
        coord_cartesian(ylim = c(y_lower, y_upper), clip = 'off') +
        NULL

    if (!is.null(facet_formula)) {
        p <- p + facet_grid(facet_formula, drop = FALSE, as.table = TRUE,
                            labeller = labeller(.cols = function(x) x))
    }

    sample_sizes <- plot_data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(facet_vars, 'additives')))) %>%
        dplyr::summarise(
            n = sum(!is.na(.data[[response_col]])),
            .groups = 'drop'
        )

    y_positions <- plot_data %>%
        dplyr::filter(!is.na(.data[[response_col]])) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(facet_vars))) %>%
        dplyr::summarise(
            y_min = min(.data[[response_col]], na.rm = TRUE),
            y_max = max(.data[[response_col]], na.rm = TRUE),
            .groups = 'drop'
        ) %>%
        dplyr::mutate(
            y_range = y_max - y_min,
            gap = ifelse(is.finite(y_range) & y_range > 0,
                         y_range * 0.16,
                         pmax(abs(y_max) * 0.06, 1)),
            label_y = ifelse(is.finite(y_max), y_max + gap, NA_real_)
        )

    sample_sizes <- sample_sizes %>%
        safe_left_join(y_positions, by = facet_vars) %>%
        dplyr::mutate(label = paste0('n=', n)) %>%
        dplyr::filter(!is.na(label_y))

    sample_sizes <- sample_sizes %>%
        dplyr::mutate(label_y = pmin(label_y, extended_y_max - (extended_y_max - y_lower) * 0.06))

    if (nrow(sample_sizes) > 0) {
        p <- p + geom_text(
            data = sample_sizes,
            aes(x = additives, y = label_y, label = label),
            inherit.aes = FALSE,
            size = 3,
            color = 'gray30'
        )
    }

    sample_top <- NULL
    if (nrow(sample_sizes) > 0) {
        sample_top <- sample_sizes %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(facet_vars))) %>%
            dplyr::summarise(sample_top = max(label_y), .groups = 'drop')
    }

    if (!is.null(stats)) {
        group_vars <- intersect(c("buffer", "growing_condition"), names(stats))
        group_vars <- intersect(group_vars, names(y_positions))
        label_data <- prepare_additive_labels(plot_data, stats,
                                              response_col = response_col,
                                              group_vars = group_vars)
        if (!is.null(label_data) && nrow(label_data) > 0) {
            label_data <- label_data %>%
                safe_left_join(dplyr::select(y_positions, -label_y), by = group_vars) %>%
                safe_left_join(sample_top, by = group_vars) %>%
                dplyr::mutate(
                    y_range = y_max - y_min,
                    label_y = pmin(label_y, extended_y_max - (extended_y_max - y_lower) * 0.12),
                    gap_min = ifelse(is.finite(y_range) & y_range > 0,
                                     y_range * 0.08,
                                     pmax(abs(y_max) * 0.05, 1)),
                    gap_sep = ifelse(is.finite(y_range) & y_range > 0,
                                     y_range * 0.12,
                                     pmax(abs(y_max) * 0.06, 1.5))
                )

            if (!"sample_top" %in% names(label_data)) {
                label_data$sample_top <- NA_real_
            }

            label_data <- label_data %>%
                dplyr::mutate(
                    label_y = ifelse(!is.na(sample_top),
                                     pmin(label_y, sample_top - gap_sep),
                                     label_y),
                    min_y = y_max + gap_min,
                    max_y_allowed = ifelse(!is.na(sample_top), sample_top - gap_sep,
                                           extended_y_max - (extended_y_max - y_lower) * 0.12),
                    max_y_allowed = ifelse(is.finite(min_y) & is.finite(max_y_allowed) & max_y_allowed < min_y,
                                           min_y,
                                           max_y_allowed),
                    label_y = pmin(pmax(label_y, min_y), max_y_allowed),
                    y_segment = label_y * 0.93
                )

            label_sig <- label_data %>% dplyr::filter(significance != 'ns')
            label_ns <- label_data %>% dplyr::filter(significance == 'ns')

            if (nrow(label_ns) > 0) {
                p <- p + geom_segment(
                    data = label_ns,
                    aes(x = 1, xend = 2, y = y_segment, yend = y_segment),
                    inherit.aes = FALSE,
                    color = 'black',
                    linewidth = 0.3
                )
            }

            if (nrow(label_sig) > 0) {
                p <- p + geom_segment(
                    data = label_sig,
                    aes(x = 1, xend = 2, y = y_segment, yend = y_segment),
                    inherit.aes = FALSE,
                    color = 'black',
                    linewidth = 0.5
                )
            }

            # Add label expression with p-value
            label_data <- label_data %>%
                dplyr::mutate(
                    p_formatted = ifelse(p_value < 0.001,
                                        sprintf('%.2e', p_value),
                                        sprintf('%.3f', p_value)),
                    label_expr = sprintf('"%s"~(italic(p)==%s)', significance, p_formatted)
                )

            p <- p + geom_text(
                data = label_data,
                aes(x = 1.5, y = label_y, label = label_expr),
                inherit.aes = FALSE,
                size = 4,
                parse = TRUE
            )
        }
    }

    if (!is.null(panel_label)) {
        use_limits <- !is.null(y_limits) && length(y_limits) == 2 && all(is.finite(y_limits))
        y_top <- if (use_limits) y_limits[2] else max(plot_data[[response_col]], na.rm = TRUE)
        range_y <- if (use_limits) diff(y_limits) else diff(range(plot_data[[response_col]], na.rm = TRUE))
        if (!is.finite(range_y) || range_y <= 0) range_y <- 1
        label_y <- y_top - range_y * 0.06

        x_left <- min(as.numeric(plot_data$additives), na.rm = TRUE)
        if (!is.finite(x_left)) x_left <- 1
        label_x <- x_left - 0.46

        p <- p + annotate(
            'text',
            x = label_x,
            y = label_y,
            label = panel_label,
            fontface = 'bold',
            hjust = 0,
            vjust = 0,
            size = 5
        )
    }

    if (!is.null(plot_tag)) {
        p <- p +
            labs(tag = plot_tag) +
            theme(plot.tag = element_text(size = 16, face = 'bold'),
                  plot.tag.position = c(0.02, 0.98))
    }

    if (is.null(image_paths)) return(p)

    image_strip <- build_live_dead_image_strip(plot_data, image_paths, facet_formula, facet_vars)
    if (is.null(image_strip)) return(p)

    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   plot.margin = margin(t = 10, r = 20, b = 5, l = 10, unit = 'pt'))

    p / image_strip + patchwork::plot_layout(heights = c(1, 0.34))
}

# nolint end
