#!/usr/bin/env Rscript

#  Fig-S4D.R
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


#  Load necessary libraries
library(tidyverse)
library(scales)


#  Define function(s)
#  Function to read data and perform Fisher's exact test
run_fishers_exact_test <- function(file_path) {
    #  Read the table from the TSV file
    data <- readr::read_tsv(file_path, show_col_types = FALSE)
    
    #  Create a contingency table from the data
    contingency_table <- xtabs(Count ~ Genotype + Status, data = data)
    
    #  Run Fisher's exact test
    fishers_test_result <- fisher.test(contingency_table)
    
    #  Return the result
    return(fishers_test_result)
}


#  Function to read data, calculate proportions, and plot them
plot_proportional_stacked_bars <- function(file_path) {
    # file_path <- pachytene_file_path
    # file_path <- diplotene_file_path
    
    #  Read the table from the TSV file
    data <- readr::read_tsv(file_path, show_col_types = FALSE)
    
    #  Adjust 'Genotype' column to a factor with levels set to ensure "WT" is
    #+ first, followed by "KO"
    data$Genotype <- factor(data$Genotype, levels = c("WT", "KO"))
    
    #  Calculate proportions
    data <- data %>%
        dplyr::group_by(Genotype) %>%
        dplyr::mutate(Total = sum(Count)) %>%
        ungroup() %>%
        dplyr::mutate(Proportion = Count / Total)
    
    #  Plot the proportional stacked bar chart
    ggplot2::ggplot(data, aes(x = Genotype, y = Proportion, fill = Status)) +
        geom_bar(stat = "identity", position = "fill") +
        scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
        labs(
            x = "Genotype",
            y = "Percentage",
            fill = "Status",
            title = paste(
                "Proportional status in", gsub(".*/|\\.tsv", "", file_path)
            )
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_blank()
        )
}


#  Define variables for paths to the TSV files
pachytene_file_path <- file.path(
    dir_script, "tsvs", "S4D_contingency-table_pachynema.tsv"
)
diplotene_file_path <- file.path(
    dir_script, "tsvs", "S4D_contingency-table_diplonema.tsv"
)


#  Run Fisher's exact tests for each table and print the results
pachytene_results <- run_fishers_exact_test(pachytene_file_path)
diplotene_results <- run_fishers_exact_test(diplotene_file_path)


#  Print the results
cat("Fisher's exact test p-value for pachytene spermatocytes:\n")
print(pachytene_results$p.value)

cat("Fisher's exact test p-value for diplotene spermatocytes:\n")
print(diplotene_results$p.value)


#  Draw plots for diplotene and pachytene spermatocytes
plot_pachytene <- plot_proportional_stacked_bars(pachytene_file_path)
plot_diplotene <- plot_proportional_stacked_bars(diplotene_file_path)

#  Print the plots
print(plot_pachytene)
print(plot_diplotene)

#  Optionally, save the plots to files
ggplot2::ggsave(
    file.path(
        dir_script, "plots", "S4D_proportional-contingency-table_diplonema.pdf"
    ),
    plot_diplotene,
    width = 8,
    height = 6
)
ggplot2::ggsave(
    file.path(
        dir_script, "plots", "S4D_proportional-contingency-table_pachynema.pdf"
    ),
    plot_pachytene,
    width = 8,
    height = 6
)
