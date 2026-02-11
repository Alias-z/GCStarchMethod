# Stomata Aperture Analysis Module
# 
# This module provides functions to analyze stomata aperture dynamics in plants.
# Supports both Arabidopsis and Species (Cowpea/Tobacco) data with:
#   - Pooled analysis (across runs) with example images
#   - Run-faceted analysis (per run) without images
#   - Both area ratio and width metrics
#
# ============================================================================================
# USAGE EXAMPLES
# ============================================================================================
#
# ARABIDOPSIS (Pooled):
#   data <- load_aperture_data(arabidopsis_aperture_file)
#   ratios <- compute_ratios(data, type = 'arabidopsis')
#   stats <- pairwise_comparisons(ratios, type = 'arabidopsis', metric = 'area')
#   plot <- plot_violin(ratios, stats, type = 'arabidopsis', metric = 'area', image_paths = img)
#
# ARABIDOPSIS (By Run):
#   plot <- plot_violin_by_run(ratios, stats, type = 'arabidopsis', metric = 'area')
#
# SPECIES (Pooled):
#   data <- load_aperture_data_species(species_aperture_file)
#   ratios <- compute_ratios(data, type = 'species')
#   stats <- pairwise_comparisons(ratios, type = 'species', metric = 'area')
#   plot <- plot_violin(ratios, stats, type = 'species', metric = 'area', image_paths = img)
#
# SPECIES (By Run):
#   plot <- plot_violin_by_run(ratios, stats, type = 'species_run', metric = 'area')

# nolint start
library('tidyverse')
library('ARTool')
library('ggthemes')
library('ggimage')
library('patchwork')

data_dir <- file.path('data')
aperture_data_dir <- file.path(data_dir, 'Stomata aperture')
arabidopsis_aperture_file <- file.path(aperture_data_dir, 'stomata_aperture_results_pool.xlsx')
species_aperture_file <- file.path(aperture_data_dir, 'stomata_aperture_cowpea_tobacco.xlsx')

# ============================================================================================
# DATA LOADING FUNCTIONS
# ============================================================================================

load_aperture_data <- function(file_path) {
    # Load Arabidopsis aperture data
    readxl::read_xlsx(file_path, na = 'NA', skip = 0) %>% 
        dplyr::mutate(
            `Buffer` = factor(normalize_buffer_value(`Buffer`),
                              levels = c('MES-BTP Buffer', 'NaOH Buffer')),
            dplyr::across(c(`image_name`, `group`, `object_category`), factor),
            `Run` = factor(`Run`),
            `Timepoint` = factor(normalize_timepoint_value(`Timepoint`),
                                 levels = c('EoN', '1hBL')),
            `Treatment` = factor(`Treatment`, levels = c('Mock', 'Mannitol')),
            `Photoperiod` = factor(normalize_photoperiod_value(`Photoperiod`),
                                   levels = c('ND', 'SD'))
        ) %>% 
        dplyr::select(dplyr::where(is.factor), `area  (μm²)`, `length (μm)`, `width (μm)`) %>% 
        dplyr::mutate_if(is.character, as.numeric)
}

load_aperture_data_species <- function(file_path) {
    # Load Species (Cowpea/Tobacco) aperture data
    readxl::read_xlsx(file_path, na = 'NA', skip = 0) %>% 
        dplyr::mutate(
            dplyr::across(c(`image_name`, `group`, `object_category`), factor),
            `Run` = factor(`Run`),
            `Species` = factor(`Species`, levels = c('V. unguiculata', 'N. tabacum')),
            `Timepoint` = factor(normalize_timepoint_value(`Timepoint`),
                                 levels = c('EoN', '1hBL')),
            `Photoperiod` = factor(normalize_photoperiod_value(`Photoperiod`),
                                   levels = c('ND', 'SD'))
        ) %>% 
        dplyr::select(dplyr::where(is.factor), `area  (μm²)`, `length (μm)`, `width (μm)`) %>% 
        dplyr::mutate_if(is.character, as.numeric)
}

# Locale-robust replacements for measurement-column loading.
resolve_measurement_columns <- function(data) {
    names_raw <- names(data)
    names_ascii <- iconv(names_raw, to = "ASCII//TRANSLIT")
    names_ascii[is.na(names_ascii)] <- names_raw[is.na(names_ascii)]
    names_norm <- tolower(gsub("[^a-z0-9]+", "_", names_ascii))

    find_measure_col <- function(prefix) {
        idx <- which(startsWith(names_norm, paste0(prefix, "_")) | names_norm == prefix)
        if (length(idx) == 0) idx <- which(grepl(paste0("^", prefix), names_norm))
        if (length(idx) == 0) return(NA_character_)
        names_raw[idx[1]]
    }

    area_col <- find_measure_col("area")
    length_col <- find_measure_col("length")
    width_col <- find_measure_col("width")

    resolved <- c(area = area_col, length = length_col, width = width_col)
    if (any(is.na(resolved))) {
        stop(
            "Could not resolve measurement columns (area/length/width). Found columns: ",
            paste(names_raw, collapse = ", ")
        )
    }

    list(area = area_col, length = length_col, width = width_col)
}

load_aperture_data <- function(file_path) {
    raw_data <- readxl::read_xlsx(file_path, na = 'NA', skip = 0)
    measurement_cols <- resolve_measurement_columns(raw_data)

    raw_data %>%
        dplyr::mutate(
            `Buffer` = factor(normalize_buffer_value(`Buffer`),
                              levels = c('MES-BTP Buffer', 'NaOH Buffer')),
            dplyr::across(c(`image_name`, `group`, `object_category`), factor),
            `Run` = factor(`Run`),
            `Timepoint` = factor(normalize_timepoint_value(`Timepoint`),
                                 levels = c('EoN', '1hBL')),
            `Treatment` = factor(`Treatment`, levels = c('Mock', 'Mannitol')),
            `Photoperiod` = factor(normalize_photoperiod_value(`Photoperiod`),
                                   levels = c('ND', 'SD'))
        ) %>%
        dplyr::rename(
            area_um2 = dplyr::all_of(measurement_cols$area),
            length_um = dplyr::all_of(measurement_cols$length),
            width_um = dplyr::all_of(measurement_cols$width)
        ) %>%
        dplyr::select(dplyr::where(is.factor), area_um2, length_um, width_um) %>%
        dplyr::mutate(dplyr::across(c(area_um2, length_um, width_um), as.numeric))
}

load_aperture_data_species <- function(file_path) {
    raw_data <- readxl::read_xlsx(file_path, na = 'NA', skip = 0)
    measurement_cols <- resolve_measurement_columns(raw_data)

    raw_data %>%
        dplyr::mutate(
            dplyr::across(c(`image_name`, `group`, `object_category`), factor),
            `Run` = factor(`Run`),
            `Species` = factor(`Species`, levels = c('V. unguiculata', 'N. tabacum')),
            `Timepoint` = factor(normalize_timepoint_value(`Timepoint`),
                                 levels = c('EoN', '1hBL')),
            `Photoperiod` = factor(normalize_photoperiod_value(`Photoperiod`),
                                   levels = c('ND', 'SD'))
        ) %>%
        dplyr::rename(
            area_um2 = dplyr::all_of(measurement_cols$area),
            length_um = dplyr::all_of(measurement_cols$length),
            width_um = dplyr::all_of(measurement_cols$width)
        ) %>%
        dplyr::select(dplyr::where(is.factor), area_um2, length_um, width_um) %>%
        dplyr::mutate(dplyr::across(c(area_um2, length_um, width_um), as.numeric))
}

# ============================================================================================
# RATIO COMPUTATION FUNCTIONS
# ============================================================================================

compute_ratios <- function(aperture_data, type = c('arabidopsis', 'species')) {
    # Compute area and width ratios for paired stoma/outer ledge measurements
    type <- match.arg(type)
    
    base_data <- aperture_data %>%
        dplyr::filter(`object_category` %in% c('outer ledge', 'stoma')) %>%
        dplyr::filter(!is.na(`Timepoint`)) %>%
        dplyr::group_by(`image_name`, `group`) %>%
        dplyr::filter(n() == 2, dplyr::n_distinct(`object_category`) == 2) %>%
        dplyr::summarise(
            `Run` = dplyr::first(`Run`),
            `Timepoint` = dplyr::first(`Timepoint`),
            `Photoperiod` = dplyr::first(`Photoperiod`),
            outer_ledge_area = dplyr::first(`area  (μm²)`[`object_category` == 'outer ledge']),
            outer_ledge_width = dplyr::first(`width (μm)`[`object_category` == 'outer ledge']),
            stoma_area = dplyr::first(`area  (μm²)`[`object_category` == 'stoma']),
            stoma_width = dplyr::first(`width (μm)`[`object_category` == 'stoma']),
            # Type-specific columns
            `Buffer` = if (type == 'arabidopsis') dplyr::first(`Buffer`) else NA,
            `Treatment` = if (type == 'arabidopsis') dplyr::first(`Treatment`) else NA,
            `Species` = if (type == 'species') dplyr::first(`Species`) else NA,
            .groups = 'drop'
        ) %>%
        dplyr::mutate(
            area_ratio_percent = outer_ledge_area / stoma_area * 100,
            width_ratio_percent = outer_ledge_width / stoma_width * 100
        )
    
    if (type == 'arabidopsis') {
        base_data %>% dplyr::select(-`Species`)
    } else {
        base_data %>% dplyr::select(-`Buffer`, -`Treatment`)
    }
}

extract_width_data <- function(aperture_data, type = c('arabidopsis', 'species')) {
    # Extract outer ledge width measurements
    type <- match.arg(type)
    
    base_data <- aperture_data %>%
        dplyr::filter(`object_category` == 'outer ledge') %>%
        dplyr::filter(!is.na(`Timepoint`)) %>%
        dplyr::rename(outer_ledge_width = `width (μm)`)
    
    if (type == 'arabidopsis') {
        base_data %>% dplyr::select(`Run`, `Buffer`, `Treatment`, `Timepoint`, `Photoperiod`, 
                                     `image_name`, `group`, outer_ledge_width)
    } else {
        base_data %>% dplyr::select(`Run`, `Species`, `Timepoint`, `Photoperiod`,
                                     `image_name`, `group`, outer_ledge_width)
    }
}

# Locale-robust replacements for derived measurement calculations.
compute_ratios <- function(aperture_data, type = c('arabidopsis', 'species')) {
    type <- match.arg(type)

    base_data <- aperture_data %>%
        dplyr::filter(`object_category` %in% c('outer ledge', 'stoma')) %>%
        dplyr::filter(!is.na(`Timepoint`)) %>%
        dplyr::group_by(`image_name`, `group`) %>%
        dplyr::filter(n() == 2, dplyr::n_distinct(`object_category`) == 2) %>%
        dplyr::summarise(
            `Run` = dplyr::first(`Run`),
            `Timepoint` = dplyr::first(`Timepoint`),
            `Photoperiod` = dplyr::first(`Photoperiod`),
            outer_ledge_area = dplyr::first(area_um2[`object_category` == 'outer ledge']),
            outer_ledge_width = dplyr::first(width_um[`object_category` == 'outer ledge']),
            stoma_area = dplyr::first(area_um2[`object_category` == 'stoma']),
            stoma_width = dplyr::first(width_um[`object_category` == 'stoma']),
            `Buffer` = if (type == 'arabidopsis') dplyr::first(`Buffer`) else NA,
            `Treatment` = if (type == 'arabidopsis') dplyr::first(`Treatment`) else NA,
            `Species` = if (type == 'species') dplyr::first(`Species`) else NA,
            .groups = 'drop'
        ) %>%
        dplyr::mutate(
            area_ratio_percent = outer_ledge_area / stoma_area * 100,
            width_ratio_percent = outer_ledge_width / stoma_width * 100
        )

    if (type == 'arabidopsis') {
        base_data %>% dplyr::select(-`Species`)
    } else {
        base_data %>% dplyr::select(-`Buffer`, -`Treatment`)
    }
}

extract_width_data <- function(aperture_data, type = c('arabidopsis', 'species')) {
    type <- match.arg(type)

    base_data <- aperture_data %>%
        dplyr::filter(`object_category` == 'outer ledge') %>%
        dplyr::filter(!is.na(`Timepoint`)) %>%
        dplyr::rename(outer_ledge_width = width_um)

    if (type == 'arabidopsis') {
        base_data %>% dplyr::select(`Run`, `Buffer`, `Treatment`, `Timepoint`, `Photoperiod`,
                                    `image_name`, `group`, outer_ledge_width)
    } else {
        base_data %>% dplyr::select(`Run`, `Species`, `Timepoint`, `Photoperiod`,
                                    `image_name`, `group`, outer_ledge_width)
    }
}

# ============================================================================================
# STATISTICAL FUNCTIONS
# ============================================================================================

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

normalize_timepoint_value <- function(value) {
    raw <- tolower(trimws(as.character(value)))
    raw[is.na(raw)] <- ""
    raw <- gsub("\\s+", "", raw)
    vapply(seq_along(raw), function(i) {
        item <- raw[i]
        if (item %in% c("eon", "endofnight")) return("EoN")
        if (grepl("^1h?bl", item)) return("1hBL")
        # if (grepl("^3h?bl", item)) return("3hBL")
        if (item == "") return(NA_character_)
        value[i]
    }, character(1))
}

normalize_buffer_value <- function(value) {
    raw <- tolower(trimws(as.character(value)))
    raw[is.na(raw)] <- ""
    vapply(seq_along(raw), function(i) {
        item <- raw[i]
        if (item %in% c("santelia", "santelia buffer", "mes-btp", "mes btp", "mesbtp", "mes-btp buffer")) {
            return("MES-BTP Buffer")
        }
        if (item %in% c("daloso", "daloso buffer", "naoh", "naoh buffer")) {
            return("NaOH Buffer")
        }
        if (item == "") return(NA_character_)
        value[i]
    }, character(1))
}

pairwise_comparisons <- function(data, type = c('arabidopsis', 'species', 'arabidopsis_run', 'species_run'),
                                  metric = c('area', 'width')) {
    # Perform pairwise Wilcoxon tests for timepoints within groups
    type <- match.arg(type)
    metric <- match.arg(metric)
    
    metric_col <- if (metric == 'area') 'area_ratio_percent' else 'outer_ledge_width'
    
    # Define grouping based on type
    if (type == 'arabidopsis') {
        group_vars <- c('Buffer', 'Treatment', 'Photoperiod')
    } else if (type == 'species') {
        group_vars <- c('Species', 'Photoperiod')
    } else if (type == 'arabidopsis_run') {
        group_vars <- c('Run', 'Buffer', 'Treatment', 'Photoperiod')
    } else {  # species_run
        group_vars <- c('Species', 'Run', 'Photoperiod')
    }
    
    combinations <- data %>% dplyr::select(dplyr::all_of(group_vars)) %>% dplyr::distinct()
    pairwise_results <- list()
    
    cat('--- PAIRWISE TIMEPOINT COMPARISONS (', toupper(metric), ') ---\n', sep = '')
    
    for (i in 1:nrow(combinations)) {
        # Build filter condition
        filter_expr <- TRUE
        group_name_parts <- c()
        for (var in group_vars) {
            val <- as.character(combinations[[var]][i])
            filter_expr <- filter_expr & (data[[var]] == val)
            group_name_parts <- c(group_name_parts, val)
        }
        group_name <- paste(group_name_parts, collapse = '_')
        
        group_data <- data[filter_expr, ]
        timepoints <- unique(group_data$Timepoint)
        
        if (length(timepoints) >= 2) {
            cat('Group:', paste(group_name_parts, collapse = ' × '), '\n')
            pairs <- combn(timepoints, 2, simplify = FALSE)
            group_pairwise <- list()
            
            for (pair in pairs) {
                tp1 <- pair[1]; tp2 <- pair[2]
                d1 <- group_data[group_data$Timepoint == tp1, ][[metric_col]]
                d2 <- group_data[group_data$Timepoint == tp2, ][[metric_col]]
                n1 <- sum(!is.na(d1))
                n2 <- sum(!is.na(d2))
                
                test <- wilcox.test(d1, d2, exact = FALSE)
                comp_name <- paste(tp1, 'vs', tp2)
                group_pairwise[[comp_name]] <- list(
                    statistic = unname(test$statistic),
                    p_value = test$p.value,
                    timepoint_1 = as.character(tp1),
                    timepoint_2 = as.character(tp2),
                    n_1 = n1,
                    n_2 = n2
                )
                
                cat('  ', comp_name, ': n =', n1, '/', n2,
                    ', W =', round(test$statistic, 3),
                    ', p =', format.pval(test$p.value), '\n')
            }
            pairwise_results[[group_name]] <- group_pairwise
            cat('\n')
        }
    }
    return(pairwise_results)
}


# ============================================================================================
# ANNOTATION HELPER FUNCTIONS
# ============================================================================================

prepare_annotations <- function(data, pairwise_results, type, metric = 'area') {
    # Prepare significance annotations for plotting
    if (length(pairwise_results) == 0) return(NULL)
    
    metric_col <- if (metric == 'area') 'area_ratio_percent' else 'outer_ledge_width'
    
    # Parse group variables based on type
    parse_group <- function(group_name, type) {
        parts <- strsplit(group_name, '_')[[1]]
        if (type == 'arabidopsis') {
            list(Buffer = parts[1], Treatment = parts[2], Photoperiod = parts[3])
        } else if (type == 'species') {
            photoperiod_val <- NA_character_
            if (length(parts) >= 2 && parts[length(parts)] %in% c('SD', 'ND')) {
                photoperiod_val <- parts[length(parts)]
                species_val <- paste(parts[1:(length(parts) - 1)], collapse = '_')
            } else {
                species_val <- group_name
            }
            list(Species = species_val, Photoperiod = photoperiod_val)
        } else if (type == 'arabidopsis_run') {
            list(Run = parts[1], Buffer = parts[2], Treatment = parts[3], Photoperiod = parts[4])
        } else {  # species_run - format is "Species_Run X"
            split_pos <- regexpr('_Run ', group_name)
            list(Species = substr(group_name, 1, split_pos - 1),
                 Run = paste0('Run ', substr(group_name, split_pos + 5, nchar(group_name))))
        }
    }
    
    annotations <- data.frame()
    tp_positions <- c('EoN' = 1, '1hBL' = 2)
    
    for (group_name in names(pairwise_results)) {
        parsed <- parse_group(group_name, type)
        
        # Filter data for this group
        filter_cond <- rep(TRUE, nrow(data))
        for (nm in names(parsed)) {
            if (!is.na(parsed[[nm]])) {
                filter_cond <- filter_cond & (data[[nm]] == parsed[[nm]])
            }
        }
        group_data <- data[filter_cond, ]
        max_y <- max(group_data[[metric_col]], na.rm = TRUE)
        
        comparisons <- pairwise_results[[group_name]]
        comp_count <- 0
        
        for (comp_name in names(comparisons)) {
            comparison <- comparisons[[comp_name]]
            tps <- c(comparison$timepoint_1, comparison$timepoint_2)
            sig <- get_significance_symbol(comparison$p_value)
            
            row <- data.frame(
                timepoint1 = tps[1], timepoint2 = tps[2],
                x_start = tp_positions[tps[1]], x_end = tp_positions[tps[2]],
                p_value = comparison$p_value, significance = sig,
                y_position = max_y * 1.15 + comp_count * (max_y * 0.12),
                max_y = max_y,
                stringsAsFactors = FALSE
            )
            # Add group columns
            for (nm in names(parsed)) row[[nm]] <- parsed[[nm]]
            
            annotations <- rbind(annotations, row)
            comp_count <- comp_count + 1
        }
    }

    # Ensure facet variables in annotation layers keep the same ordering/levels as the main data.
    # Otherwise ggplot2 may re-order panels alphabetically when it merges facet values across layers.
    for (key in intersect(names(annotations), names(data))) {
        if (is.factor(data[[key]])) {
            annotations[[key]] <- factor(annotations[[key]], levels = levels(data[[key]]))
        }
    }
    return(annotations)
}

prepare_sample_size_labels <- function(data, metric_col, group_vars, x_var) {
    if (length(group_vars) == 0) return(NULL)

    count_vars <- unique(c(group_vars, x_var))
    sample_sizes <- data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(count_vars))) %>%
        dplyr::summarise(
            n = sum(!is.na(.data[[metric_col]])),
            .groups = 'drop'
        )

    y_positions <- data %>%
        dplyr::filter(!is.na(.data[[metric_col]])) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
        dplyr::summarise(
            y_min = min(.data[[metric_col]], na.rm = TRUE),
            y_max = max(.data[[metric_col]], na.rm = TRUE),
            .groups = 'drop'
        ) %>%
        dplyr::mutate(
            y_range = y_max - y_min,
            gap = ifelse(is.finite(y_range) & y_range > 0,
                         y_range * 0.25,
                         pmax(abs(y_max) * 0.08, 0.4)),
            label_y = ifelse(is.finite(y_max), y_max + gap, NA_real_)
        )

    sample_sizes %>%
        dplyr::left_join(y_positions, by = group_vars) %>%
        dplyr::mutate(label = paste0('n=', n)) %>%
        dplyr::filter(!is.na(label_y))
}

filter_violin_data <- function(data, metric_col, panel_vars, x_var) {
    data <- data %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(c(panel_vars, x_var)))) %>%
        dplyr::filter(dplyr::n() > 1, dplyr::n_distinct(.data[[metric_col]]) > 1) %>%
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
# PLOTTING FUNCTIONS - POOLED (with images)
# ============================================================================================

plot_violin <- function(data, pairwise_results = NULL, type = c('arabidopsis', 'species'),
                        metric = c('area', 'width'), image_paths = NULL, photoperiod_filter = NULL,
                        facet_buffer = TRUE, facet_photoperiod = TRUE, facet_treatment = TRUE,
                        plot_tag = NULL,
                        y_limits = NULL,
                        facet_labeller = NULL,
                        panel_letters = NULL,
                        panel_letter_mode = c('inside', 'outside'),
                        strip_text_face = 'italic') {
    # Create violin plot for pooled data with optional images
    type <- match.arg(type)
    metric <- match.arg(metric)
    
    metric_col <- if (metric == 'area') 'area_ratio_percent' else 'outer_ledge_width'
    y_label <- if (metric == 'area') 'outer ledge / stoma area ratio (%) \n' else 'outer ledge width (um) \n'
    
    # Filter by photoperiod if specified
    if (!is.null(photoperiod_filter)) {
        data <- data %>% dplyr::filter(`Photoperiod` == photoperiod_filter)
    }

    plot_data <- data %>%
        dplyr::filter(is.finite(.data[[metric_col]]), !is.na(.data[['Timepoint']]))
    
    # Calculate y-axis range
    y_max_data <- max(plot_data[[metric_col]], na.rm = TRUE)
    y_min_data <- min(plot_data[[metric_col]], na.rm = TRUE)
    y_range <- y_max_data - y_min_data
    extra_top <- ifelse(is.finite(y_range) && y_range > 0,
                        y_range * ifelse(type == 'species', 0.6, 0.35),
                        abs(y_max_data) * 0.2 + 0.01)

    y_lower <- y_min_data * 0.95
    y_upper <- y_max_data + extra_top
    if (!is.null(y_limits) && length(y_limits) == 2 && all(is.finite(y_limits))) {
        y_lower <- y_limits[1]
        y_upper <- y_limits[2]
    }
    extended_y_max <- y_upper
    
    # Timepoint colors
    tp_colors <- c('EoN' = 'gray', '1hBL' = 'steelblue')
    tp_labels <- c('EoN' = 'EoN', '1hBL' = '1h BL')

    photoperiod_levels <- NULL
    if ('Photoperiod' %in% names(plot_data)) {
        photoperiod_levels <- levels(plot_data$Photoperiod)
    }

    facet_vars <- c()
    facet_formula <- NULL
    if (type == 'arabidopsis') {
        if (facet_buffer) facet_vars <- c(facet_vars, 'Buffer')
        if (facet_photoperiod) facet_vars <- c(facet_vars, 'Photoperiod')
        if (facet_treatment) facet_vars <- c(facet_vars, 'Treatment')

        row_term <- if (facet_buffer) 'Buffer' else '.'
        col_terms <- c()
        if (facet_photoperiod) col_terms <- c(col_terms, 'Photoperiod')
        if (facet_treatment) col_terms <- c(col_terms, 'Treatment')
        if (length(col_terms) == 0) col_terms <- '.'
        facet_formula <- stats::as.formula(paste(row_term, "~", paste(col_terms, collapse = " + ")))
    } else {
        facet_vars <- c(facet_vars, 'Species')
        if (facet_photoperiod) facet_vars <- c(facet_vars, 'Photoperiod')

        col_term <- if (facet_photoperiod) 'Photoperiod' else '.'
        facet_formula <- stats::as.formula(paste('Species', "~", col_term))
    }

    # Keep ordering consistent with the manuscript: Mock before Mannitol.
    if (type == 'arabidopsis' && "Treatment" %in% names(plot_data)) {
        plot_data$Treatment <- factor(plot_data$Treatment, levels = c('Mock', 'Mannitol'))
    }
    sample_sizes <- prepare_sample_size_labels(plot_data, metric_col, facet_vars, 'Timepoint')
    violin_data <- filter_violin_data(plot_data, metric_col, facet_vars, 'Timepoint')
    
    if (is.null(facet_labeller)) {
        facet_labeller <- labeller(.cols = function(x) x)
    }

    panel_letter_mode <- match.arg(panel_letter_mode)
    panel_letter_data <- NULL

    # Base plot
    p <- ggplot(plot_data, aes(x = `Timepoint`, y = .data[[metric_col]], fill = `Timepoint`))
    if (nrow(violin_data) > 0) {
        p <- p + geom_violin(data = violin_data, color = 'black', width = 0.9)
    }
    p <- p +
        geom_jitter(color = 'black', width = 0.05, size = 0.5, alpha = 0.7) +
        stat_summary(fun = mean, geom = 'point', shape = 23, size = 3, fill = 'white', color = 'white') +
        scale_fill_manual(values = tp_colors, labels = tp_labels) +
        scale_x_discrete(labels = tp_labels) +
        labs(y = y_label) +
        theme_bw() +
        theme(
            panel.background = element_blank(),
            panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
            panel.spacing = unit(1.2, 'lines'),
            legend.position = 'none',
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text = element_text(size = 16),
            strip.text = element_text(size = 16, face = strip_text_face),
            strip.background = element_rect(fill = 'grey98'),
            plot.title = element_blank(),
            plot.margin = margin(t = 10, r = 20, b = 40, l = 10, unit = 'pt')
        ) +
        coord_cartesian(ylim = c(y_lower, y_upper), clip = 'off')

    if (!is.null(panel_letters)) {
        letter_values <- unname(panel_letters)
        letter_names <- names(panel_letters)
        letter_var <- NULL

        if (!is.null(letter_names)) {
            if ("Treatment" %in% names(plot_data) &&
                all(letter_names %in% levels(plot_data$Treatment))) {
                letter_var <- "Treatment"
            } else if ("Photoperiod" %in% names(plot_data) &&
                       all(letter_names %in% levels(plot_data$Photoperiod))) {
                letter_var <- "Photoperiod"
            }
        }

        if (is.null(letter_var)) {
            if ("Treatment" %in% names(plot_data)) {
                letter_var <- "Treatment"
                letter_names <- levels(plot_data$Treatment)
            } else if ("Photoperiod" %in% names(plot_data)) {
                letter_var <- "Photoperiod"
                letter_names <- levels(plot_data$Photoperiod)
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
    
    p <- p + facet_grid(facet_formula, drop = FALSE, as.table = TRUE,
                        labeller = facet_labeller)

    if (!is.null(sample_sizes) && nrow(sample_sizes) > 0) {
        p <- p + geom_text(
            data = sample_sizes,
            aes(x = .data[['Timepoint']], y = label_y, label = label),
            inherit.aes = FALSE,
            size = 3,
            color = 'gray30'
        )
    }

    sample_tops <- NULL
    if (!is.null(sample_sizes) && nrow(sample_sizes) > 0) {
        sample_tops <- sample_sizes %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(facet_vars))) %>%
            dplyr::summarise(sample_top = max(label_y), .groups = 'drop')
    }

    # Add significance annotations
    if (!is.null(pairwise_results)) {
        p <- add_significance_annotations(p, plot_data, pairwise_results, type, metric,
                                          extended_y_max, sample_tops, facet_vars)
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

    image_strip <- build_timepoint_image_strip(plot_data,
                                               type = type,
                                               image_paths = image_paths,
                                               facet_formula = facet_formula,
                                               facet_vars = facet_vars,
                                               tp_labels = tp_labels)
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

build_timepoint_image_strip <- function(plot_data,
                                        type,
                                        image_paths,
                                        facet_formula,
                                        facet_vars,
                                        tp_labels) {
    if (is.null(image_paths)) return(NULL)
    if (is.null(plot_data) || nrow(plot_data) == 0) return(NULL)
    if (!"Timepoint" %in% names(plot_data)) return(NULL)

    timepoint_levels <- levels(plot_data$Timepoint)
    if (is.null(timepoint_levels) || length(timepoint_levels) == 0) return(NULL)

    panel_map <- plot_data %>%
        dplyr::distinct(dplyr::across(dplyr::all_of(facet_vars))) %>%
        tidyr::crossing(Timepoint = factor(timepoint_levels, levels = timepoint_levels))

    image_rows <- list()
    if (type == 'arabidopsis') {
        for (photoperiod in names(image_paths)) {
            for (buffer in names(image_paths[[photoperiod]])) {
                for (treatment in names(image_paths[[photoperiod]][[buffer]])) {
                    for (tp in names(image_paths[[photoperiod]][[buffer]][[treatment]])) {
                        path <- image_paths[[photoperiod]][[buffer]][[treatment]][[tp]]
                        if (is.null(path)) next
                        image_rows[[length(image_rows) + 1]] <- data.frame(
                            Photoperiod = photoperiod,
                            Buffer = buffer,
                            Treatment = treatment,
                            Timepoint = tp,
                            image_path = path,
                            stringsAsFactors = FALSE
                        )
                    }
                }
            }
        }
    } else {
        for (species in names(image_paths)) {
            for (tp in names(image_paths[[species]])) {
                path <- image_paths[[species]][[tp]]
                if (is.null(path)) next
                image_rows[[length(image_rows) + 1]] <- data.frame(
                    Species = species,
                    Timepoint = tp,
                    image_path = path,
                    stringsAsFactors = FALSE
                )
            }
        }
    }

    image_data <- dplyr::bind_rows(image_rows)
    if (nrow(image_data) == 0) return(NULL)

    # Keep only images matching the data shown in the plot (even if those vars are not faceted).
    for (key in c("Buffer", "Photoperiod", "Treatment", "Species")) {
        if (key %in% names(image_data) && key %in% names(plot_data)) {
            allowed <- unique(as.character(plot_data[[key]]))
            image_data <- image_data %>% dplyr::filter(.data[[key]] %in% allowed)
        }
    }

    if (nrow(image_data) == 0) return(NULL)

    # Match facet factor levels so empty panels still appear.
    for (key in facet_vars) {
        if (key %in% names(plot_data) && key %in% names(image_data) && is.factor(plot_data[[key]])) {
            image_data[[key]] <- factor(image_data[[key]], levels = levels(plot_data[[key]]))
        }
    }
    image_data$Timepoint <- factor(image_data$Timepoint, levels = timepoint_levels)

    ggplot(panel_map, aes(x = Timepoint, y = 0.5)) +
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
        scale_x_discrete(labels = tp_labels) +
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

add_significance_annotations <- function(p, data, pairwise_results, type, metric, extended_y_max,
                                         sample_tops = NULL, facet_vars = NULL) {
    # Add significance annotations to plot
    metric_col <- if (metric == 'area') 'area_ratio_percent' else 'outer_ledge_width'
    annotations <- prepare_annotations(data, pairwise_results, type, metric)
    
    if (is.null(annotations) || nrow(annotations) == 0) return(p)

    facet_join <- character(0)
    if (!is.null(facet_vars)) {
        facet_join <- intersect(facet_vars, names(annotations))
    }
    if (length(facet_join) == 0) {
        facet_join <- intersect(c("Buffer", "Photoperiod", "Treatment", "Species"), names(annotations))
    }

    y_positions <- data %>%
        dplyr::filter(is.finite(.data[[metric_col]])) %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(facet_join))) %>%
        dplyr::summarise(
            y_min_data = min(.data[[metric_col]], na.rm = TRUE),
            y_max_data = max(.data[[metric_col]], na.rm = TRUE),
            y_range = y_max_data - y_min_data,
            .groups = 'drop'
        )

    annotations <- annotations %>%
        dplyr::left_join(y_positions, by = facet_join) %>%
        dplyr::mutate(
            y_position = pmin(y_position, extended_y_max * 0.8),
            gap_min = ifelse(is.finite(y_range) & y_range > 0,
                             y_range * 0.08,
                             pmax(abs(y_max_data) * 0.05, 0.25)),
            gap_sep = ifelse(is.finite(y_range) & y_range > 0,
                             y_range * 0.12,
                             pmax(abs(y_max_data) * 0.06, 0.35))
        )

    if (!is.null(sample_tops) && nrow(sample_tops) > 0) {
        join_vars <- intersect(names(sample_tops), names(annotations))
        if (length(join_vars) > 0) {
            annotations <- annotations %>%
                dplyr::left_join(sample_tops, by = join_vars) %>%
                dplyr::mutate(
                    y_position = ifelse(!is.na(sample_top),
                                        pmin(y_position, sample_top - gap_sep),
                                        y_position)
                )
        }
    }

    if (!"sample_top" %in% names(annotations)) {
        annotations$sample_top <- NA_real_
    }

    annotations <- annotations %>%
        dplyr::mutate(
            min_y = y_max_data + gap_min,
            y_position = ifelse(is.finite(max_y),
                                pmax(y_position, max_y + gap_min * 0.4),
                                y_position),
            y_segment = y_position * 0.98
        )

    annotations <- annotations %>%
        dplyr::mutate(
            max_y_allowed = ifelse(!is.na(sample_top), sample_top - gap_sep, extended_y_max * 0.8),
            max_y_allowed = ifelse(is.finite(min_y) & is.finite(max_y_allowed) & max_y_allowed < min_y,
                                   min_y,
                                   max_y_allowed),
            y_position = pmin(pmax(y_position, min_y), max_y_allowed),
            y_segment = y_position * 0.93
        )
    
    # Split significant vs non-significant and add label expressions
    ann_sig <- annotations %>%
        dplyr::filter(significance != 'ns') %>%
        dplyr::mutate(
            p_formatted = ifelse(p_value < 0.001,
                                sprintf('%.2e', p_value),
                                sprintf('%.3f', p_value)),
            label_expr = sprintf('"%s"~(italic(p)==%s)', significance, p_formatted)
        )
    ann_ns <- annotations %>%
        dplyr::filter(significance == 'ns') %>%
        dplyr::mutate(
            p_formatted = ifelse(p_value < 0.001,
                                sprintf('%.2e', p_value),
                                sprintf('%.3f', p_value)),
            label_expr = sprintf('"%s"~(italic(p)==%s)', significance, p_formatted)
        )

    # Add non-significant (gray)
    if (nrow(ann_ns) > 0) {
        p <- p +
            geom_segment(data = ann_ns, aes(x = x_start, xend = x_end,
                        y = y_segment, yend = y_segment),
                        inherit.aes = FALSE, color = 'black', linewidth = 0.3) +
            geom_text(data = ann_ns, aes(x = (x_start + x_end) / 2, y = y_position, label = label_expr),
                     inherit.aes = FALSE, size = 4, color = 'black', parse = TRUE)
    }

    # Add significant (black, bold)
    if (nrow(ann_sig) > 0) {
        p <- p +
            geom_segment(data = ann_sig, aes(x = x_start, xend = x_end,
                        y = y_segment, yend = y_segment),
                        inherit.aes = FALSE, color = 'black', linewidth = 0.5) +
            geom_text(data = ann_sig, aes(x = (x_start + x_end) / 2, y = y_position, label = label_expr),
                     inherit.aes = FALSE, size = 4, color = 'black', parse = TRUE)
    }
    
    return(p)
}


# ============================================================================================
# PLOTTING FUNCTIONS - BY RUN (no images)
# ============================================================================================

plot_violin_by_run <- function(data, pairwise_results = NULL, type = c('arabidopsis_run', 'species_run'),
                                metric = c('area', 'width'), photoperiod_filter = NULL,
                                y_limits = NULL) {
    # Create violin plot faceted by Run (no images)
    type <- match.arg(type)
    metric <- match.arg(metric)
    
    metric_col <- if (metric == 'area') 'area_ratio_percent' else 'outer_ledge_width'
    y_label <- if (metric == 'area') 'outer ledge / stoma area ratio (%) \n' else 'outer ledge width (um) \n'
    
    # Filter by photoperiod if specified
    if (!is.null(photoperiod_filter)) {
        data <- data %>% dplyr::filter(`Photoperiod` == photoperiod_filter)
    }

    plot_data <- data %>%
        dplyr::filter(is.finite(.data[[metric_col]]), !is.na(.data[['Timepoint']]))

    # Calculate y-axis range
    y_max_data <- max(plot_data[[metric_col]], na.rm = TRUE)
    y_min_data <- min(plot_data[[metric_col]], na.rm = TRUE)
    y_range <- y_max_data - y_min_data

    y_lower <- y_min_data * 0.95
    y_upper <- y_max_data + (y_range * 1.5)
    if (!is.null(y_limits) && length(y_limits) == 2 && all(is.finite(y_limits))) {
        y_lower <- y_limits[1]
        y_upper <- y_limits[2]
    }
    extended_y_max <- y_upper
    
    # Timepoint colors
    tp_colors <- c('EoN' = 'gray', '1hBL' = 'steelblue')
    tp_labels <- c('EoN' = 'EoN', '1hBL' = '1h BL')
    
    if (type == 'arabidopsis_run') {
        facet_vars <- c('Buffer', 'Run', 'Photoperiod', 'Treatment')
    } else {
        facet_vars <- c('Run', 'Photoperiod', 'Species')
    }
    violin_data <- filter_violin_data(plot_data, metric_col, facet_vars, 'Timepoint')
    
    # Base plot
    p <- ggplot(plot_data, aes(x = `Timepoint`, y = .data[[metric_col]], fill = `Timepoint`))
    if (nrow(violin_data) > 0) {
        p <- p + geom_violin(data = violin_data, color = 'black', width = 0.9)
    }
    p <- p +
        geom_jitter(color = 'black', width = 0.05, size = 0.5, alpha = 0.7) +
        stat_summary(fun = mean, geom = 'point', shape = 23, size = 3, fill = 'white', color = 'white') +
        scale_fill_manual(values = tp_colors, labels = tp_labels) +
        scale_x_discrete(labels = tp_labels) +
        labs(y = y_label) +
        theme_bw() +
        theme(
            panel.background = element_blank(),
            panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
            panel.spacing = unit(1.2, 'lines'),
            legend.position = 'none',
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text = element_text(size = 16),
            strip.text = element_text(size = 12, face = 'italic'),
            strip.background = element_rect(fill = 'grey98'),
            plot.title = element_blank(),
            plot.margin = margin(t = 10, r = 20, b = 40, l = 10, unit = 'pt')
        ) +
        coord_cartesian(ylim = c(y_lower, y_upper))
    
    # Add faceting based on type
    if (type == 'arabidopsis_run') {
        p <- p + facet_grid(`Buffer` + `Run` ~ `Photoperiod` + `Treatment`,
                            drop = TRUE, as.table = TRUE)
    } else {  # species_run
        p <- p + facet_grid(`Run` ~ `Photoperiod` + `Species`, drop = FALSE, as.table = TRUE)
    }

    sample_sizes <- prepare_sample_size_labels(plot_data, metric_col, facet_vars, 'Timepoint')
    if (!is.null(sample_sizes) && nrow(sample_sizes) > 0) {
        p <- p + geom_text(
            data = sample_sizes,
            aes(x = .data[['Timepoint']], y = label_y, label = label),
            inherit.aes = FALSE,
            size = 3,
            color = 'gray30'
        )
    }
    
    sample_tops <- NULL
    if (!is.null(sample_sizes) && nrow(sample_sizes) > 0) {
        sample_tops <- sample_sizes %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(facet_vars))) %>%
            dplyr::summarise(sample_top = max(label_y), .groups = 'drop')
    }

    # Add significance annotations
    if (!is.null(pairwise_results)) {
        p <- add_significance_annotations_by_run(p, plot_data, pairwise_results, type, metric,
                                                 extended_y_max, sample_tops)
    }
    
    return(p)
}

add_significance_annotations_by_run <- function(p, data, pairwise_results, type, metric, extended_y_max,
                                                sample_tops = NULL) {
    # Add significance annotations for by-run plots with proper facet matching
    metric_col <- if (metric == 'area') 'area_ratio_percent' else 'outer_ledge_width'
    tp_positions <- c('EoN' = 1, '1hBL' = 2)
    
    annotations <- data.frame()
    
    for (group_name in names(pairwise_results)) {
        # Parse group name based on type
        if (type == 'arabidopsis_run') {
            # Format: "Run X_Buffer_Treatment_Photoperiod"
            parts <- strsplit(group_name, '_')[[1]]
            run_val <- parts[1]
            buffer_val <- parts[2]
            treatment_val <- parts[3]
            photoperiod_val <- parts[4]
            
            # Filter data for this group
            group_data <- data %>% 
                dplyr::filter(Run == run_val, Buffer == buffer_val, 
                             Treatment == treatment_val, Photoperiod == photoperiod_val)
        } else {  # species_run
            # Format: "Species_Run X"
            split_pos <- regexpr('_Run ', group_name)
            species_val <- substr(group_name, 1, split_pos - 1)
            after_run <- substr(group_name, split_pos + 5, nchar(group_name))
            run_parts <- strsplit(after_run, '_')[[1]]
            run_val <- paste0('Run ', run_parts[1])
            photoperiod_val <- if (length(run_parts) >= 2 && run_parts[2] %in% c('SD', 'ND')) {
                run_parts[2]
            } else {
                NA_character_
            }
            
            group_data <- data %>% dplyr::filter(Species == species_val, Run == run_val)
            if (!is.na(photoperiod_val)) {
                group_data <- group_data %>% dplyr::filter(Photoperiod == photoperiod_val)
            }
        }
        
        if (nrow(group_data) == 0) next
        
        max_y <- max(group_data[[metric_col]], na.rm = TRUE)
        comparisons <- pairwise_results[[group_name]]
        comp_count <- 0
        
        for (comp_name in names(comparisons)) {
            comparison <- comparisons[[comp_name]]
            tps <- c(comparison$timepoint_1, comparison$timepoint_2)
            sig <- get_significance_symbol(comparison$p_value)
            
            row <- data.frame(
                x_start = tp_positions[tps[1]], x_end = tp_positions[tps[2]],
                significance = sig,
                p_value = comparison$p_value,
                y_position = max_y * 1.15 + comp_count * (max_y * 0.12),
                max_y = max_y,
                stringsAsFactors = FALSE
            )
            
            if (type == 'arabidopsis_run') {
                row$Buffer <- buffer_val
                row$Run <- run_val
                row$Photoperiod <- photoperiod_val
                row$Treatment <- treatment_val
            } else {
                row$Species <- species_val
                row$Run <- run_val
                row$Photoperiod <- photoperiod_val
            }
            
            annotations <- rbind(annotations, row)
            comp_count <- comp_count + 1
        }
    }
    
    if (nrow(annotations) == 0) return(p)
    
    # Adjust y positions
    annotations <- annotations %>%
        dplyr::mutate(y_position = pmin(y_position, extended_y_max * 0.65))

    if (!is.null(sample_tops) && nrow(sample_tops) > 0) {
        join_vars <- intersect(names(sample_tops), names(annotations))
        if (length(join_vars) > 0) {
            gap <- extended_y_max * 0.05
            annotations <- annotations %>%
                dplyr::left_join(sample_tops, by = join_vars) %>%
                dplyr::mutate(
                    y_position = ifelse(!is.na(sample_top),
                                        pmin(y_position, sample_top - gap),
                                        y_position),
                    y_position = ifelse(is.finite(max_y),
                                        pmax(y_position, max_y + gap * 0.5),
                                        y_position)
                )
        }
    }
    
    # Split significant vs non-significant and add label expressions
    ann_sig <- annotations %>%
        dplyr::filter(significance != 'ns') %>%
        dplyr::mutate(
            p_formatted = ifelse(p_value < 0.001,
                                sprintf('%.2e', p_value),
                                sprintf('%.3f', p_value)),
            label_expr = sprintf('"%s"~(italic(p)==%s)', significance, p_formatted)
        )
    ann_ns <- annotations %>%
        dplyr::filter(significance == 'ns') %>%
        dplyr::mutate(
            p_formatted = ifelse(p_value < 0.001,
                                sprintf('%.2e', p_value),
                                sprintf('%.3f', p_value)),
            label_expr = sprintf('"%s"~(italic(p)==%s)', significance, p_formatted)
        )

    # Add non-significant (gray)
    if (nrow(ann_ns) > 0) {
        p <- p +
            geom_segment(data = ann_ns, aes(x = x_start, xend = x_end,
                        y = y_position * 0.98, yend = y_position * 0.98),
                        inherit.aes = FALSE, color = 'black', linewidth = 0.3) +
            geom_text(data = ann_ns, aes(x = (x_start + x_end) / 2, y = y_position, label = label_expr),
                     inherit.aes = FALSE, size = 4, color = 'black', parse = TRUE)
    }

    # Add significant (black, bold)
    if (nrow(ann_sig) > 0) {
        p <- p +
            geom_segment(data = ann_sig, aes(x = x_start, xend = x_end,
                        y = y_position * 0.98, yend = y_position * 0.98),
                        inherit.aes = FALSE, color = 'black', linewidth = 0.5) +
            geom_text(data = ann_sig, aes(x = (x_start + x_end) / 2, y = y_position, label = label_expr),
                     inherit.aes = FALSE, size = 4, color = 'black', parse = TRUE)
    }
    
    return(p)
}

# ============================================================================================
# CONVENIENCE WRAPPER FUNCTIONS
# ============================================================================================

# Arabidopsis - Pooled Analysis
analyze_arabidopsis_pooled <- function(file_path = arabidopsis_aperture_file,
                                       image_paths = NULL, photoperiod_filter = NULL,
                                       area_y_limits = NULL, width_y_limits = NULL) {
    data <- load_aperture_data(file_path)
    ratios <- compute_ratios(data, type = 'arabidopsis')
    width_data <- extract_width_data(data, type = 'arabidopsis')
    
    area_stats <- pairwise_comparisons(ratios, type = 'arabidopsis', metric = 'area')
    width_stats <- pairwise_comparisons(width_data, type = 'arabidopsis', metric = 'width')
    
    area_plot <- plot_violin(ratios, area_stats, type = 'arabidopsis', metric = 'area',
                             image_paths = image_paths, photoperiod_filter = photoperiod_filter,
                             y_limits = area_y_limits)
    width_plot <- plot_violin(width_data, width_stats, type = 'arabidopsis', metric = 'width',
                              image_paths = image_paths, photoperiod_filter = photoperiod_filter,
                              y_limits = width_y_limits)
    
    list(data = data, ratios = ratios, width_data = width_data,
         area_stats = area_stats, width_stats = width_stats,
         area_plot = area_plot, width_plot = width_plot)
}

# Arabidopsis - By Run Analysis
analyze_arabidopsis_by_run <- function(file_path = arabidopsis_aperture_file,
                                       photoperiod_filter = NULL,
                                       area_y_limits = NULL, width_y_limits = NULL) {
    data <- load_aperture_data(file_path)
    ratios <- compute_ratios(data, type = 'arabidopsis')
    width_data <- extract_width_data(data, type = 'arabidopsis')
    
    area_stats <- pairwise_comparisons(ratios, type = 'arabidopsis_run', metric = 'area')
    width_stats <- pairwise_comparisons(width_data, type = 'arabidopsis_run', metric = 'width')
    
    area_plot <- plot_violin_by_run(ratios, area_stats, type = 'arabidopsis_run', metric = 'area',
                                    photoperiod_filter = photoperiod_filter,
                                    y_limits = area_y_limits)
    width_plot <- plot_violin_by_run(width_data, width_stats, type = 'arabidopsis_run', metric = 'width',
                                     photoperiod_filter = photoperiod_filter,
                                     y_limits = width_y_limits)
    
    list(data = data, ratios = ratios, width_data = width_data,
         area_stats = area_stats, width_stats = width_stats,
         area_plot = area_plot, width_plot = width_plot)
}

# Species - Pooled Analysis
analyze_species_pooled <- function(file_path = species_aperture_file,
                                   image_paths = NULL, photoperiod_filter = NULL,
                                   area_y_limits = NULL, width_y_limits = NULL) {
    data <- load_aperture_data_species(file_path)
    ratios <- compute_ratios(data, type = 'species')
    width_data <- extract_width_data(data, type = 'species')
    
    area_stats <- pairwise_comparisons(ratios, type = 'species', metric = 'area')
    width_stats <- pairwise_comparisons(width_data, type = 'species', metric = 'width')
    
    area_plot <- plot_violin(ratios, area_stats, type = 'species', metric = 'area',
                             image_paths = image_paths, photoperiod_filter = photoperiod_filter,
                             y_limits = area_y_limits)
    width_plot <- plot_violin(width_data, width_stats, type = 'species', metric = 'width',
                              image_paths = image_paths, photoperiod_filter = photoperiod_filter,
                              y_limits = width_y_limits)
    
    list(data = data, ratios = ratios, width_data = width_data,
         area_stats = area_stats, width_stats = width_stats,
         area_plot = area_plot, width_plot = width_plot)
}

# Species - By Run Analysis
analyze_species_by_run <- function(file_path = species_aperture_file,
                                   photoperiod_filter = NULL,
                                   area_y_limits = NULL, width_y_limits = NULL) {
    data <- load_aperture_data_species(file_path)
    ratios <- compute_ratios(data, type = 'species')
    width_data <- extract_width_data(data, type = 'species')
    
    area_stats <- pairwise_comparisons(ratios, type = 'species_run', metric = 'area')
    width_stats <- pairwise_comparisons(width_data, type = 'species_run', metric = 'width')
    
    area_plot <- plot_violin_by_run(ratios, area_stats, type = 'species_run', metric = 'area',
                                    photoperiod_filter = photoperiod_filter,
                                    y_limits = area_y_limits)
    width_plot <- plot_violin_by_run(width_data, width_stats, type = 'species_run', metric = 'width',
                                     photoperiod_filter = photoperiod_filter,
                                     y_limits = width_y_limits)
    
    list(data = data, ratios = ratios, width_data = width_data,
         area_stats = area_stats, width_stats = width_stats,
         area_plot = area_plot, width_plot = width_plot)
}

# nolint end
