# Leaf Starch Analysis Module
#
# This module analyzes leaf starch content under different photoperiod conditions.
# Data: Leaf starch measurements (μmol/g FW) for ND (12/12) vs SD (8/16) photoperiods.

# nolint start
library('tidyverse')
library('readxl')

leaf_starch_file <- file.path('data', 'Leaf starch', 'Leaf_starch_finals.xlsx')

#' Load leaf starch data from Excel file
#'
#' @param file_path Path to the Excel file
#' @return Tibble with columns: photoperiod, starch_content
#' @export
load_leaf_starch_data <- function(file_path = leaf_starch_file) {
    raw_data <- readxl::read_excel(file_path, sheet = 'Sheet1')

    # Reshape from wide to long format
    data <- tibble::tibble(
        photoperiod = c(
            rep('ND', nrow(raw_data)),
            rep('SD', nrow(raw_data))
        ),
        starch_content = c(
            raw_data[[1]],  # 12/12 column
            raw_data[[2]]   # 8/16 column
        )
    ) %>%
        dplyr::mutate(
            photoperiod = factor(photoperiod, levels = c('ND', 'SD'))
        )

    data
}

#' Compare leaf starch between photoperiods using Wilcoxon test
#'
#' @param data Data frame with photoperiod and starch_content columns
#' @return Tibble with comparison statistics
#' @export
compare_photoperiods <- function(data) {
    nd_data <- data %>% dplyr::filter(photoperiod == 'ND') %>% dplyr::pull(starch_content)
    sd_data <- data %>% dplyr::filter(photoperiod == 'SD') %>% dplyr::pull(starch_content)

    test_result <- wilcox.test(nd_data, sd_data, alternative = 'two.sided')

    tibble::tibble(
        comparison = 'ND vs SD',
        p_value = test_result$p.value,
        w_statistic = test_result$statistic,
        n_nd = length(nd_data),
        n_sd = length(sd_data),
        mean_nd = mean(nd_data, na.rm = TRUE),
        mean_sd = mean(sd_data, na.rm = TRUE),
        median_nd = median(nd_data, na.rm = TRUE),
        median_sd = median(sd_data, na.rm = TRUE)
    )
}

#' Plot leaf starch content by photoperiod
#'
#' @param data Data frame with photoperiod and starch_content columns
#' @param stats Statistics tibble from compare_photoperiods()
#' @param y_limits Y-axis limits (default: NULL for automatic)
#' @param panel_label Optional panel label to draw in the top-left corner
#' @return ggplot object
#' @export
plot_leaf_starch <- function(data, stats = NULL, y_limits = NULL, panel_label = NULL) {
    photoperiod_colors <- c('ND' = 'gray', 'SD' = 'gray')
    photoperiod_labels <- c('ND' = 'ND', 'SD' = 'SD')

    # Create base plot
    p <- ggplot(data, aes(x = photoperiod, y = starch_content, fill = photoperiod)) +
        geom_violin(color = 'black', width = 0.9) +
        geom_jitter(width = 0.05, height = 0, size = 0.5, alpha = 0.7, color = 'black') +
        stat_summary(fun = mean, geom = 'point',
                     shape = 23, size = 3, fill = 'white', color = 'white') +
        scale_fill_manual(values = photoperiod_colors, labels = photoperiod_labels) +
        scale_x_discrete(labels = photoperiod_labels) +
        labs(
            x = NULL,
            y = 'leaf starch content (umol/g FW)'
        ) +
        theme_bw() +
        theme(
            panel.background = element_blank(),
            panel.border = element_rect(color = 'black', fill = NA, linewidth = 1),
            panel.spacing = unit(1.2, 'lines'),
            legend.position = 'none',
            axis.title.x = element_blank(),
            axis.title.y = element_text(size = 18),
            axis.text = element_text(size = 16),
            plot.title = element_blank(),
            plot.margin = margin(t = 10, r = 20, b = 40, l = 10, unit = 'pt')
        )

    # Apply y-axis limits if provided
    if (!is.null(y_limits)) {
        p <- p + coord_cartesian(ylim = y_limits)
    }

    # Add sample size labels (n) for each photoperiod group
    sample_sizes <- data %>%
        dplyr::group_by(photoperiod) %>%
        dplyr::summarise(
            n = sum(!is.na(starch_content)),
            y_min = min(starch_content, na.rm = TRUE),
            y_max = max(starch_content, na.rm = TRUE),
            .groups = 'drop'
        ) %>%
        dplyr::mutate(
            y_range = y_max - y_min,
            gap = ifelse(is.finite(y_range) & y_range > 0,
                         y_range * 0.18,
                         pmax(abs(y_max) * 0.06, 0.3)),
            label_y = y_max + gap,
            label = paste0('n=', n)
        )

    if (!is.null(y_limits) && length(y_limits) == 2 && all(is.finite(y_limits))) {
        sample_sizes <- sample_sizes %>%
            dplyr::mutate(label_y = pmin(label_y, y_limits[2] * 0.9))
    }

    if (nrow(sample_sizes) > 0) {
        p <- p + geom_text(
            data = sample_sizes,
            aes(x = photoperiod, y = label_y, label = label),
            inherit.aes = FALSE,
            size = 3,
            color = 'gray30'
        )
    }

    # Add significance annotation if stats provided
    if (!is.null(stats) && nrow(stats) > 0) {
        p_val <- stats$p_value[1]
        significance <- if (p_val < 0.001) {
            '***'
        } else if (p_val < 0.01) {
            '**'
        } else if (p_val < 0.05) {
            '*'
        } else {
            'ns'
        }

        # Calculate y position for annotation
        # Position above the tallest data point across ALL groups
        y_max_per_group <- data %>%
            dplyr::group_by(photoperiod) %>%
            dplyr::summarise(max_val = max(starch_content, na.rm = TRUE), .groups = 'drop')

        y_max_data <- max(y_max_per_group$max_val, na.rm = TRUE)

        # Add 10% above the max
        label_y <- y_max_data * 1.1

        # Segment position slightly below label
        y_segment <- label_y * 0.93

        # Constrain to y-axis limits if provided
        if (!is.null(y_limits)) {
            # Only constrain if the natural position would exceed limits
            if (label_y > y_limits[2]) {
                label_y <- y_limits[2] * 0.95
                y_segment <- label_y * 0.93
            }
        }

        # Add segment (different style for ns vs significant)
        if (significance == 'ns') {
            p <- p + annotate(
                'segment',
                x = 1, xend = 2,
                y = y_segment, yend = y_segment,
                color = 'black',
                linewidth = 0.3
            )
        } else {
            p <- p + annotate(
                'segment',
                x = 1, xend = 2,
                y = y_segment, yend = y_segment,
                color = 'black',
                linewidth = 0.5
            )
        }

        # Add significance label with p-value (italic p)
        p_formatted <- ifelse(p_val < 0.001,
                             sprintf('%.2e', p_val),
                             sprintf('%.3f', p_val))
        label_expr <- sprintf('"%s"~(italic(p)==%s)', significance, p_formatted)
        p <- p + annotate(
            'text',
            x = 1.5,
            y = label_y,
            label = label_expr,
            size = 4,
            parse = TRUE
        )
    }

    if (!is.null(panel_label)) {
        use_limits <- !is.null(y_limits) && length(y_limits) == 2 && all(is.finite(y_limits))
        y_top <- if (use_limits) y_limits[2] else max(data$starch_content, na.rm = TRUE)
        range_y <- if (use_limits) diff(y_limits) else diff(range(data$starch_content, na.rm = TRUE))
        if (!is.finite(range_y) || range_y <= 0) range_y <- 1
        label_y <- y_top - range_y * 0.06

        x_left <- min(as.numeric(data$photoperiod), na.rm = TRUE)
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

    p
}

#' Analyze leaf starch data (main workflow)
#'
#' @param y_limits Y-axis limits for plot
#' @param panel_label Optional panel label to draw in the plot
#' @return List with data, stats, and plot
#' @export
analyze_leaf_starch <- function(y_limits = NULL, panel_label = NULL) {
    data <- load_leaf_starch_data()
    stats <- compare_photoperiods(data)
    plot <- plot_leaf_starch(data, stats = stats, y_limits = y_limits,
                             panel_label = panel_label)

    list(
        data = data,
        stats = stats,
        plot = plot
    )
}

# nolint end
