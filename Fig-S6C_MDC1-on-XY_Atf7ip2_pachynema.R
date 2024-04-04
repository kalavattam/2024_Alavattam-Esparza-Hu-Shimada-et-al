#!/usr/bin/env Rscript

#  Fig-S6C_MDC1-on-XY_Atf7ip2_pachynema.R
#  KA

#  Determine the directory the script is being run from
dir_script <- if (
    interactive() && requireNamespace("rstudioapi", quietly = TRUE)
) {
    dirname(rstudioapi::getSourceEditorContext()$path)
} else {
    args <- commandArgs(trailingOnly = FALSE)
    arg_script <- args[grep("^--file=", args)]
    path_script <- sub("^--file=", "", arg_script[length(arg_script)]) 
    normalizePath(path_script, mustWork = FALSE)
    dirname(path_script)
}


#  Define function(s)
#  Define the function to plot with error bars
plot_with_error_bars <- function(summary_data, error_type = "sd") {
    # Initialize the plot with genotype and mean percentage
    plot <- ggplot2::ggplot(
        summary_data, 
        aes(x = genotype, y = mean_percentage, fill = genotype)
    ) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    labs(x = "Genotype", y = "Mean Percentage") +
    theme_minimal()
    
    #  Determine which error metric to use and add error bars if applicable
    if (error_type == "sd") {
        #  Filter out groups where SD is 0 to avoid adding error bars for them
        data_with_error <- filter(summary_data, sd > 0)
        if (nrow(data_with_error) > 0) {
            plot <- plot + geom_errorbar(
                data = data_with_error,
                aes(ymin = mean_percentage - sd, ymax = mean_percentage + sd),
                position = position_dodge(width = 0.7),
                width = 0.25
            )
        }
    } else if (error_type == "sem") {
        #  Filter out groups where SEM is 0 to avoid adding error bars for them
        data_with_error <- filter(summary_data, sem > 0)
        if (nrow(data_with_error) > 0) {
            plot <- plot + geom_errorbar(
                data = data_with_error,
                aes(ymin = mean_percentage - sem, ymax = mean_percentage + sem),
                position = position_dodge(width = 0.7),
                width = 0.25
            )
        }
    }
    
    return(plot)
}


#  Load necessary libraries
library(tidyverse)


#  Read the TSV-file observational data into a data frame, and wrangle it a bit
file_path <- file.path(
    dir_script, "tsvs", "Fig-S6C_MDC1-on-XY_Atf7ip2_pachynema.tsv"
)
data <- readr::read_tsv(file_path, show_col_types = FALSE) %>%
    #  Filter out rows where the total number of observations is zero
    dplyr::filter(S6C_observed_pachynema > 0) %>%
    #  Calculate the percentage
    dplyr::mutate(
        percentage =
            `S6C_MDC1-on-XY_Atf7ip2_pachynema` /
            S6C_observed_pachynema * 100,
        genotype = factor(genotype, levels = c("WT", "KO"))
    )

#  Group data by genotype and replicate, then summarize
grouped_data <- data %>%
    dplyr::group_by(genotype, replicate) %>%
    dplyr::summarize(mean_percentage = mean(percentage, na.rm = TRUE))

#  Perform an unpaired t-test between WT and KO
#  Extract the mean percentages for WT and KO
wt_means <- grouped_data %>%
    dplyr::filter(genotype == "WT") %>%
    dplyr::pull(mean_percentage)

ko_means <- grouped_data %>%
    dplyr::filter(genotype == "KO") %>%
    dplyr::pull(mean_percentage)

#  Perform unpaired two-tailed t-test
t_test_result <- t.test(wt_means, ko_means, var.equal = TRUE)

#  Output the p-value from the t-test
t_test_result$p.value  # [1] 0.2722284

#  Summarize data across replicates for each genotype to calculate overall
#+ mean, SD, and SEM
summary_data <- data %>%
    #  Ensure observations are valid
    dplyr::filter(S6C_observed_pachynema > 0) %>%
    dplyr::mutate(
        percentage =
            `S6C_MDC1-on-XY_Atf7ip2_pachynema` /
            S6C_observed_pachynema * 100
    ) %>%
    dplyr::group_by(genotype) %>%
    dplyr::summarize(
        mean_percentage = mean(percentage, na.rm = TRUE),
        sd = sd(percentage, na.rm = TRUE),
        n = sum(S6C_observed_pachynema), .groups = 'drop'
    ) %>%
    dplyr::mutate(
        sem = sd / sqrt(n)
    )

#  Check if 'sd' and 'sem' are calculated correctly
print(summary_data)

#  Plot mean ± SD
plot_sd <- plot_with_error_bars(summary_data, "sd")

#  Plot mean ± SEM
plot_sem <- plot_with_error_bars(summary_data, "sem")

#  Optionally, print the plots to the console
print(plot_sd)
print(plot_sem)

#  Save the plots if needed
ggplot2::ggsave(
    file.path(
        dir_script, "plots", "Fig-S6C_MDC1-on-XY_Atf7ip2_pachynema_SD.pdf"
    ),
    plot = plot_sd, width = 8, height = 6, dpi = 300
)

ggplot2::ggsave(
    file.path(
        dir_script, "plots", "Fig-S6C_MDC1-on-XY_Atf7ip2_pachynema_SEM.pdf"
    ),
    plot = plot_sem, width = 8, height = 6, dpi = 300
)
