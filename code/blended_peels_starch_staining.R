# Blended Peels Starch Staining Analysis Module
#
# This module provides functions to analyze starch staining data from blended peels.
# It standardizes metadata from folder names, supports simple timepoint comparisons,
# and produces publication-ready violin plots.
#
# ============================================================================================
# USAGE EXAMPLES
# ============================================================================================
#
# Load datasets:
#   data_12 <- load_blended_peels_starch_data(
#       file.path("data", "BL blended peels starch staining", "12-12-summary", "1212-summary.csv"),
#       photoperiod = "ND"
#   )
#   data_8 <- load_blended_peels_starch_data(
#       file.path("data", "BL blended peels starch staining", "8-16-summary", "Daloso-816-summary.csv"),
#       photoperiod = "SD"
#   )
#   data <- dplyr::bind_rows(data_12, data_8)
#
# Compare timepoints and plot:
#   stats <- compare_timepoints_by_group(data)
#   plot <- plot_blended_peels_starch(data, stats = stats, facet_photoperiod = TRUE)
#
# nolint start
library('tidyverse')
library('ggimage')
library('patchwork')

timepoint_colors <- c('EoN' = 'gray', '1hBL' = 'steelblue')
timepoint_labels <- c('EoN' = 'EoN', '1hBL' = '1h BL')

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

build_blended_example_images <- function(example_dir = file.path('data',
                                                                'BL blended peels starch staining',
                                                                'examples')) {
    if (!dir.exists(example_dir)) return(NULL)
    files <- list.files(example_dir, pattern = "\\.(tif|png)$", full.names = TRUE, ignore.case = TRUE)
    if (length(files) == 0) return(NULL)

    image_paths <- list()

    for (path in files) {
        name <- tolower(tools::file_path_sans_ext(basename(path)))
        parts <- strsplit(name, "_")[[1]]
        if (length(parts) < 4) next

        photoperiod <- if (parts[1] == "nd") {
            "ND"
        } else if (parts[1] == "sd") {
            "SD"
        } else {
            NA_character_
        }

        condition <- if (parts[2] == "mock") {
            "Mock"
        } else if (parts[2] == "mannitol") {
            "Mannitol"
        } else {
            NA_character_
        }

        timepoint <- if (parts[3] == "eon") {
            "EoN"
        } else if (parts[3] == "1hbl") {
            "1hBL"
        } else {
            NA_character_
        }

        buffer <- if (parts[4] == "mes") {
            "MES-BTP Buffer"
        } else if (parts[4] == "naoh") {
            "NaOH Buffer"
        } else {
            NA_character_
        }

        if (any(is.na(c(photoperiod, condition, timepoint, buffer)))) next

        if (is.null(image_paths[[photoperiod]])) image_paths[[photoperiod]] <- list()
        if (is.null(image_paths[[photoperiod]][[buffer]])) {
            image_paths[[photoperiod]][[buffer]] <- list()
        }
        if (is.null(image_paths[[photoperiod]][[buffer]][[condition]])) {
            image_paths[[photoperiod]][[buffer]][[condition]] <- list()
        }

        image_paths[[photoperiod]][[buffer]][[condition]][[timepoint]] <- path
    }

    image_paths
}

parse_blended_peels_metadata <- function(folder_name) {
    clean <- toupper(as.character(folder_name))
    clean <- gsub("_", " ", clean)
    clean <- gsub("\\s+", " ", clean)
    clean <- trimws(clean)

    prefix <- sub("^([0-9]+(?:\\.[0-9]+)?).*$", "\\1", clean)

    buffer <- dplyr::case_when(
        grepl("^1", prefix) ~ "MES-BTP Buffer",
        grepl("^2", prefix) ~ "NaOH Buffer",
        TRUE ~ NA_character_
    )

    timepoint <- dplyr::case_when(
        grepl("EON", clean) ~ "EoN",
        grepl("BL", clean) ~ "1hBL",
        TRUE ~ NA_character_
    )

    condition <- dplyr::case_when(
        prefix %in% c("1", "2") ~ "Mock",
        # prefix %in% c("1.1", "2.1") ~ "Direct BL",
        prefix %in% c("1.2", "2.2") ~ "Mannitol",
        TRUE ~ NA_character_
    )

    list(buffer = buffer, timepoint = timepoint, condition = condition, folder_prefix = prefix)
}

safe_wilcox_test <- function(values, groups) {
    if (length(unique(groups)) < 2) return(NA_real_)
    tryCatch(
        stats::wilcox.test(values ~ groups, exact = FALSE)$p.value,
        error = function(e) NA_real_
    )
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

load_blended_peels_starch_data <- function(file_path, photoperiod = NA_character_) {
    data <- readr::read_csv(file_path, show_col_types = FALSE)
    names(data) <- normalize_names(names(data))

    required_cols <- c("folder_name", "starch_guard_cell_area_ratio")
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
        stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
    }

    meta <- parse_blended_peels_metadata(data$folder_name)

    data %>%
        dplyr::mutate(
            buffer = factor(meta$buffer, levels = c("MES-BTP Buffer", "NaOH Buffer")),
            condition = factor(meta$condition,
                               levels = c("Mock",
                                          # "Direct BL",
                                          "Mannitol")),
            timepoint = factor(meta$timepoint, levels = c("EoN", "1hBL")),
            photoperiod = if (is.na(photoperiod)) NA_character_
            else normalize_photoperiod_value(photoperiod),
            starch_guard_cell_area_ratio = as.numeric(starch_guard_cell_area_ratio)
        ) %>%
        dplyr::filter(!is.na(condition))
}

# ============================================================================================
# STATISTICS
# ============================================================================================

compare_timepoints_by_group <- function(data,
                                        response_col = "starch_guard_cell_area_ratio",
                                        group_vars = c("buffer", "condition", "photoperiod")) {
    group_vars <- group_vars[group_vars %in% names(data)]

    data %>%
        dplyr::filter(is.finite(.data[[response_col]]), !is.na(timepoint)) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
        dplyr::summarise(
            p_value = safe_wilcox_test(.data[[response_col]], timepoint),
            counts = list(table(timepoint)),
            .groups = 'drop'
        ) %>%
        dplyr::mutate(
            n_eon = vapply(counts, function(x) {
                if ("EoN" %in% names(x)) unname(x[["EoN"]]) else 0L
            }, integer(1)),
            n_1hbl = vapply(counts, function(x) {
                if ("1hBL" %in% names(x)) unname(x[["1hBL"]]) else 0L
            }, integer(1)),
            significance = get_significance_symbol(p_value)
        ) %>%
        dplyr::select(-counts)
}

prepare_timepoint_labels <- function(data,
                                     stats,
                                     response_col = "starch_guard_cell_area_ratio",
                                     group_vars = c("buffer", "condition", "photoperiod")) {
    if (is.null(stats) || nrow(stats) == 0) return(NULL)
    group_vars <- group_vars[group_vars %in% names(stats)]

    y_max <- data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
        dplyr::summarise(label_y = max(.data[[response_col]], na.rm = TRUE) * 1.1,
                         .groups = 'drop')

    dplyr::left_join(stats, y_max, by = group_vars)
}

build_blended_image_strip <- function(plot_data,
                                      image_paths,
                                      facet_formula,
                                      facet_vars,
                                      use_photoperiod) {
    if (is.null(image_paths)) return(NULL)
    if (is.null(plot_data) || nrow(plot_data) == 0) return(NULL)
    if (!"timepoint" %in% names(plot_data)) return(NULL)

    timepoint_levels <- levels(plot_data$timepoint)
    if (is.null(timepoint_levels) || length(timepoint_levels) == 0) return(NULL)

    panel_map <- plot_data %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(facet_vars))) %>%
        tidyr::crossing(timepoint = factor(timepoint_levels, levels = timepoint_levels))

    image_rows <- list()
    for (photoperiod in names(image_paths)) {
        for (buffer in names(image_paths[[photoperiod]])) {
            for (condition in names(image_paths[[photoperiod]][[buffer]])) {
                for (tp in names(image_paths[[photoperiod]][[buffer]][[condition]])) {
                    path <- image_paths[[photoperiod]][[buffer]][[condition]][[tp]]
                    if (is.null(path)) next

                    row <- data.frame(
                        buffer = buffer,
                        condition = condition,
                        timepoint = tp,
                        image_path = path,
                        stringsAsFactors = FALSE
                    )
                    if (use_photoperiod) row$photoperiod <- photoperiod
                    image_rows[[length(image_rows) + 1]] <- row
                }
            }
        }
    }

    image_data <- dplyr::bind_rows(image_rows)
    if (nrow(image_data) == 0) return(NULL)

    # Keep only images matching the data shown in the plot (even if those vars are not faceted).
    for (key in c("buffer", "photoperiod", "condition")) {
        if (key %in% names(image_data) && key %in% names(plot_data)) {
            allowed <- unique(as.character(plot_data[[key]]))
            image_data <- image_data %>% dplyr::filter(.data[[key]] %in% allowed)
        }
    }

    if (nrow(image_data) == 0) return(NULL)

    # Match the facet factor levels so empty panels still appear.
    if ("buffer" %in% names(panel_map) && "buffer" %in% names(plot_data)) {
        image_data$buffer <- factor(image_data$buffer, levels = levels(plot_data$buffer))
    }
    if ("condition" %in% names(panel_map) && "condition" %in% names(plot_data)) {
        image_data$condition <- factor(image_data$condition, levels = levels(plot_data$condition))
    }
    if (use_photoperiod && "photoperiod" %in% names(panel_map) && "photoperiod" %in% names(plot_data)) {
        image_data$photoperiod <- factor(image_data$photoperiod, levels = levels(plot_data$photoperiod))
    }
    image_data$timepoint <- factor(image_data$timepoint, levels = timepoint_levels)

    ggplot(panel_map, aes(x = timepoint, y = 0.5)) +
        geom_blank() +
        ggimage::geom_image(
            data = image_data,
            aes(image = image_path),
            inherit.aes = TRUE,
            size = 0.85,
            asp = 1
        ) +
        facet_grid(facet_formula, drop = FALSE, as.table = TRUE,
                   labeller = labeller(.cols = function(x) x)) +
        scale_x_discrete(labels = timepoint_labels) +
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
}

build_panel_letter_strip <- function(letter_data,
                                     facet_formula,
                                     facet_labeller,
                                     panel_spacing = unit(1.2, 'lines'),
                                     plot_margin = margin(t = 6, r = 20, b = 0, l = 10, unit = 'pt')) {
    if (is.null(letter_data) || nrow(letter_data) == 0) return(NULL)

    ggplot(letter_data, aes(x = 0.5, y = 0.5, label = panel_label)) +
        geom_text(fontface = 'bold', size = 5, hjust = 0.5, vjust = 0.5) +
        facet_grid(facet_formula, drop = FALSE, as.table = TRUE,
                   labeller = facet_labeller) +
        scale_x_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
        scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)) +
        theme_void() +
        theme(
            strip.text = element_blank(),
            strip.background = element_blank(),
            panel.border = element_blank(),
            panel.spacing = panel_spacing,
            plot.margin = plot_margin
        )
}

# ============================================================================================
# PLOTTING
# ============================================================================================

plot_blended_peels_starch <- function(data,
                                      stats = NULL,
                                      response_col = "starch_guard_cell_area_ratio",
                                      y_limits = NULL,
                                      facet_photoperiod = TRUE,
                                      facet_condition = TRUE,
                                      image_paths = NULL,
                                      facet_buffer = TRUE,
                                      plot_tag = NULL,
                                      facet_labeller = NULL,
                                      panel_letters = NULL,
                                      panel_letter_mode = c('inside', 'outside'),
                                      strip_text_face = 'italic') {
    plot_data <- data %>%
        dplyr::filter(is.finite(.data[[response_col]]), !is.na(timepoint))

    use_photoperiod <- facet_photoperiod &&
        "photoperiod" %in% names(plot_data) &&
        any(!is.na(plot_data$photoperiod))

    if (!use_photoperiod && !is.null(image_paths) && length(image_paths) > 1) {
        plot_photoperiods <- unique(na.omit(as.character(plot_data$photoperiod)))
        if (length(plot_photoperiods) == 1 && plot_photoperiods %in% names(image_paths)) {
            image_paths <- image_paths[plot_photoperiods]
        }
    }

    if (use_photoperiod) {
        plot_data <- plot_data %>%
            dplyr::mutate(photoperiod = factor(photoperiod, levels = c('ND', 'SD')))
        facet_vars <- c()
        if (facet_buffer) facet_vars <- c(facet_vars, 'buffer')
        facet_vars <- c(facet_vars, 'photoperiod')
        if (facet_condition) facet_vars <- c(facet_vars, 'condition')

        row_term <- if (facet_buffer) 'buffer' else '.'
        col_terms <- c('photoperiod')
        if (facet_condition) col_terms <- c(col_terms, 'condition')
        facet_formula <- stats::as.formula(paste(row_term, "~", paste(col_terms, collapse = " + ")))
    } else {
        facet_vars <- c()
        if (facet_buffer) facet_vars <- c(facet_vars, 'buffer')
        if (facet_condition) facet_vars <- c(facet_vars, 'condition')

        row_term <- if (facet_buffer) 'buffer' else '.'
        col_terms <- if (facet_condition) 'condition' else '.'
        facet_formula <- stats::as.formula(paste(row_term, "~", col_terms))
    }

    y_max_data <- max(plot_data[[response_col]], na.rm = TRUE)
    y_min_data <- min(plot_data[[response_col]], na.rm = TRUE)
    y_range <- y_max_data - y_min_data
    extra_top <- ifelse(is.finite(y_range) && y_range > 0,
                        y_range * 0.35,
                        abs(y_max_data) * 0.2 + 0.01)

    y_lower <- y_min_data * 0.95
    y_upper <- y_max_data + extra_top
    if (!is.null(y_limits) && length(y_limits) == 2 && all(is.finite(y_limits))) {
        y_lower <- y_limits[1]
        y_upper <- y_limits[2]
    }
    extended_y_max <- y_upper

    violin_data <- filter_violin_data(plot_data, response_col, facet_vars, 'timepoint')

    if (is.null(facet_labeller)) {
        facet_labeller <- labeller(.cols = function(x) x)
    }

    panel_letter_mode <- match.arg(panel_letter_mode)
    panel_letter_data <- NULL

    p <- ggplot(plot_data,
                aes(x = timepoint, y = .data[[response_col]], fill = timepoint))
    if (nrow(violin_data) > 0) {
        p <- p + geom_violin(data = violin_data, color = 'black', width = 0.9)
    }
    p <- p +
        geom_jitter(color = 'black', width = 0.05, size = 0.5, alpha = 0.7) +
        stat_summary(fun = mean, geom = 'point',
                     shape = 23, size = 3, fill = 'white', color = 'white') +
        scale_fill_manual(values = timepoint_colors, labels = timepoint_labels) +
        scale_x_discrete(labels = timepoint_labels) +
        labs(x = NULL, y = 'starch guard cell area ratio (%) \n') +
        daloso_plot_theme() +
        theme(strip.text = element_text(size = 16, face = strip_text_face)) +
        facet_grid(facet_formula, drop = FALSE, as.table = TRUE,
                   labeller = facet_labeller) +
        coord_cartesian(ylim = c(y_lower, y_upper), clip = 'off')

    if (!is.null(panel_letters)) {
        letter_values <- unname(panel_letters)
        letter_names <- names(panel_letters)
        letter_var <- NULL

        if (!is.null(letter_names)) {
            if ("condition" %in% names(plot_data) &&
                all(letter_names %in% levels(plot_data$condition))) {
                letter_var <- "condition"
            } else if ("photoperiod" %in% names(plot_data) &&
                       all(letter_names %in% levels(plot_data$photoperiod))) {
                letter_var <- "photoperiod"
            }
        }

        if (is.null(letter_var)) {
            if ("condition" %in% names(plot_data)) {
                letter_var <- "condition"
                letter_names <- levels(plot_data$condition)
            } else if ("photoperiod" %in% names(plot_data)) {
                letter_var <- "photoperiod"
                letter_names <- levels(plot_data$photoperiod)
            }
        }

        if (!is.null(letter_var)) {
            if (length(letter_values) < length(letter_names)) {
                letter_names <- letter_names[seq_len(length(letter_values))]
            } else if (length(letter_values) > length(letter_names)) {
                letter_values <- letter_values[seq_len(length(letter_names))]
            }

            letter_data <- data.frame(
                panel_key = letter_names,
                panel_label = paste0("(", letter_values, ")"),
                stringsAsFactors = FALSE
            )
            letter_data$panel_key <- factor(letter_data$panel_key,
                                            levels = levels(plot_data[[letter_var]]))
            names(letter_data)[names(letter_data) == "panel_key"] <- letter_var

            facet_keys <- setdiff(facet_vars, letter_var)
            if (length(facet_keys) > 0) {
                panel_map <- plot_data %>% dplyr::distinct(dplyr::across(dplyr::all_of(facet_keys)))
                letter_data <- tidyr::crossing(panel_map, letter_data)
            }

            range_y <- y_upper - y_lower
            offset_y <- if (is.finite(range_y) && range_y > 0) range_y * 0.06 else 0.6
            x_left <- 0.5
            letter_data$x <- x_left
            if (panel_letter_mode == 'inside') {
                letter_data$label_y <- y_upper - offset_y
                p <- p +
                    geom_text(
                        data = letter_data,
                        aes(x = x, y = label_y, label = panel_label),
                        inherit.aes = FALSE,
                        fontface = 'bold',
                        size = 5,
                        hjust = 0,
                        vjust = 0,
                        nudge_x = 0.04,
                        nudge_y = 0
                    )
            } else {
                panel_letter_data <- letter_data
            }
        }
    }

    sample_sizes <- plot_data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(facet_vars, 'timepoint')))) %>%
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
                         y_range * 0.18,
                         pmax(abs(y_max) * 0.06, 0.3)),
            label_y = ifelse(is.finite(y_max), y_max + gap, NA_real_)
        )

    sample_sizes <- sample_sizes %>%
        dplyr::left_join(y_positions, by = facet_vars) %>%
        dplyr::mutate(label = paste0('n=', n)) %>%
        dplyr::filter(!is.na(label_y)) %>%
        dplyr::mutate(label_y = pmin(label_y, extended_y_max * 0.9))

    if (nrow(sample_sizes) > 0) {
        p <- p + geom_text(
            data = sample_sizes,
            aes(x = timepoint, y = label_y, label = label),
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
        group_vars <- intersect(c("buffer", "condition", "photoperiod"), names(stats))
        group_vars <- intersect(group_vars, names(y_positions))
        label_data <- prepare_timepoint_labels(plot_data, stats,
                                               response_col = response_col,
                                               group_vars = group_vars)
        if (!is.null(label_data) && nrow(label_data) > 0) {
            label_data <- label_data %>%
                dplyr::left_join(dplyr::select(y_positions, -label_y), by = group_vars) %>%
                dplyr::left_join(sample_top, by = group_vars) %>%
                dplyr::mutate(
                    y_range = y_max - y_min,
                    label_y = pmin(label_y, extended_y_max * 0.8),
                    gap_min = ifelse(is.finite(y_range) & y_range > 0,
                                     y_range * 0.08,
                                     pmax(abs(y_max) * 0.05, 0.25)),
                    gap_sep = ifelse(is.finite(y_range) & y_range > 0,
                                     y_range * 0.12,
                                     pmax(abs(y_max) * 0.06, 0.35))
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
                    max_y_allowed = ifelse(!is.na(sample_top), sample_top - gap_sep, extended_y_max * 0.8),
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

    if (!is.null(plot_tag)) {
        p <- p +
            labs(tag = plot_tag) +
            theme(plot.tag = element_text(size = 16, face = 'bold'),
                  plot.tag.position = c(0.02, 0.98))
    }

    has_outside_letters <- !is.null(panel_letter_data) && panel_letter_mode == 'outside'
    if (has_outside_letters) {
        p <- p + theme(plot.margin = margin(t = 0, r = 20, b = 40, l = 10, unit = 'pt'))
    }

    if (is.null(image_paths)) {
        if (!has_outside_letters) return(p)
        letter_strip <- build_panel_letter_strip(panel_letter_data, facet_formula, facet_labeller)
        if (is.null(letter_strip)) return(p)
        return(letter_strip / p + patchwork::plot_layout(heights = c(0.08, 1)))
    }

    image_strip <- build_blended_image_strip(plot_data,
                                             image_paths,
                                             facet_formula,
                                             facet_vars,
                                             use_photoperiod)
    if (is.null(image_strip)) {
        if (!has_outside_letters) return(p)
        letter_strip <- build_panel_letter_strip(panel_letter_data, facet_formula, facet_labeller)
        if (is.null(letter_strip)) return(p)
        return(letter_strip / p + patchwork::plot_layout(heights = c(0.08, 1)))
    }

    p <- p + theme(axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   plot.margin = margin(t = if (has_outside_letters) 0 else 10,
                                        r = 20, b = 5, l = 10, unit = 'pt'))

    if (!has_outside_letters) {
        return(p / image_strip + patchwork::plot_layout(heights = c(1, 0.32)))
    }

    letter_strip <- build_panel_letter_strip(panel_letter_data, facet_formula, facet_labeller)
    if (is.null(letter_strip)) {
        return(p / image_strip + patchwork::plot_layout(heights = c(1, 0.32)))
    }

    letter_strip / p / image_strip + patchwork::plot_layout(heights = c(0.08, 1, 0.32))
}

# nolint end
