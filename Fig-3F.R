#!/usr/bin/env Rscript

#  Fig-3F.R
#  KA

#  Initialize variable 'dir_script', which holds the directory path of the
#+ current script; the path will be used to locate and save files below
dir_script <- if (
    #  'interactive()' checks if R is running in an interactive environment 
    #+ (like RStudio); if true, then the code proceeds to check if the
    #+ 'rstudioapi' package  is available and can be loaded quietly (without
    #+ printing messages)
    interactive() && requireNamespace("rstudioapi", quietly = TRUE)
) {
    #  If running interactively in RStudio, use 'rstudioapi' to get the path of 
    #+ the current script and then use function 'dirname()' to extract the
    #+ directory part of the path; 'rstudioapi::getSourceEditorContext()$path'
    #+ gets the full path of the script being run in RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
} else {
    #  If not running interactively (like when running a script from the
    #+ command line), get the command line arguments
    args <- commandArgs(trailingOnly = FALSE)
    
    #  From the command line arguments, find the argument that starts with
    #+ '--file='; this is typically used to specify the script file being
    #+ executed
    arg_script <- args[grep("^--file=", args)]
    
    #  Extract the script path from the '--file=' argument
    path_script <- sub("^--file=", "", arg_script[length(arg_script)]) 
    
    #  Normalize the path (convert to absolute path) and handle any errors;
    #+ 'normalizePath' turns a relative path into an absolute path, and
    #+ 'mustWork = FALSE' means the function won't throw an error if the path
    #+ doesn't exist
    normalizePath(path_script, mustWork = FALSE)
    
    #  Extract the directory part of the script path
    dirname(path_script)
}

#  Load necessary libraries
library(tidyverse)
library(scales)

#  Load in data from TSV file; files of interest are stored in the base
#+ directory of the repo, which is the same location
data <- readr::read_tsv(
    file.path(
        dir_script, "tsvs",
        "wrangle-assess_blinded-scoring_Atf7ip2.SYCP1_XY_observations_run-3.tsv"
    ),
    show_col_types = FALSE
) %>%
    
    #  Convert the proportions to numeric format (if it's not already)
    dplyr::mutate(proportion = as.numeric(gsub("%", "", proportion)) / 100)

#  Further wrangle the data, calculating total sample sizes along with the
#+ summary statistics SD and SEM
aggregated_data <- data %>%
    dplyr::filter(class == 1) %>%
    dplyr::select(c(genotype, replicate, proportion)) %>%
    dplyr::group_by(genotype) %>%
    dplyr::summarize(
        mean_prop = mean(proportion, na.rm = TRUE),
        sd = sd(proportion, na.rm = TRUE),
        sem = sd / sqrt(n()),
        .groups = 'drop'
    ) %>%
    dplyr::mutate(genotype = factor(genotype, levels = c("WT", "KO")))

#  Plot with SD bars
plot_sd <- ggplot2::ggplot(
    aggregated_data, aes(x = genotype, y = mean_prop, fill = genotype)
) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(
        aes(ymin = mean_prop - sd, ymax = mean_prop + sd),
        width = 0.25,
        position = position_dodge(0.9)
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Mean proportion (%)", x = "Genotype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

#  Plot with SEM bars
plot_sem <- ggplot2::ggplot(
    aggregated_data, aes(x = genotype, y = mean_prop, fill = genotype)
) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(
        aes(ymin = mean_prop - sem, ymax = mean_prop + sem),
        width = 0.25,
        position = position_dodge(0.9)
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Mean proportion (%)", x = "Genotype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

print_plots <- TRUE
if (isTRUE(print_plots)) {
    print(plot_sd)
    print(plot_sem)
}

save_plots <- TRUE
if (isTRUE(save_plots)) {
    #  Save the plot with SD bars (will be saved to $HOME directory)
    ggplot2::ggsave(
        "Fig-3F_prop-SYCP1-XY_SD_run-3.pdf",
        plot = plot_sd,
        width = 8.5,
        height = 6,
        dpi = 300
    )
    
    #  Save the plot with SEM bars (will be saved to $HOME directory)
    ggplot2::ggsave(
        "Fig-3F_prop-SYCP1-XY_SEM_run-3.pdf",
        plot = plot_sem,
        width = 8.5,
        height = 6,
        dpi = 300
    )
}

#  Run t-tests for each substage
substage_levels <- unique(data$substage)  # Get unique substages
t_test_results <- list()  # Initialize a list to store the results

for (sub in substage_levels) {
    #  Filter data for WT and KO within the current substage
    data_wt <- data %>%
        dplyr::filter(class == 1 & substage == sub & genotype == "WT") %>%
        dplyr::pull(proportion)  # Extract the sample_size column
  
    data_ko <- data %>%
        dplyr::filter(class == 1 & substage == sub & genotype == "KO") %>%
        dplyr::pull(proportion)  # Extract the sample_size column
  
  #  Perform the t-test
  t_test <- t.test(data_wt, data_ko, var.equal = TRUE)

  #  Store the result
  t_test_results[[sub]] <- t_test
}

#  Assess p-values
check_p_values <- TRUE
if (isTRUE(check_p_values)) {
    print(t_test_results$P$p.value)
}

#  The p-values from t-tests are as follows:
#+ - P 0.4852286 (run 1)
#+ - P 0.1959663 (run 2)
#+ - P 0.001089535 (run 3)
