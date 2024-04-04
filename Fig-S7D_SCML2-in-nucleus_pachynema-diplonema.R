#!/usr/bin/env Rscript

#  Fig-S7D_SCML2-in-nucleus_pachynema-diplonema.R
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
#  Function to plot data with error bars for specified substages
plot_with_error_bars <- function(summary_data, substages, error_type = "sd") {
    plots <- list()  # Initialize list to store plots for each substage

    for (substage in substages) {
        #  Filter summary_data for the current substage
        substage_data <- filter(summary_data, substage == !!substage)
      
        #  Initialize plot
        plot <- ggplot2::ggplot(
            substage_data, 
            aes(x = genotype, y = mean_percentage, fill = genotype)
        ) +
        geom_bar(
            stat = "identity",
            position = position_dodge(width = 0.7),
            width = 0.6
        ) +
        labs(
            x = "Genotype",
            y = "Mean Percentage",
            title = paste("Substage", substage)
        ) +
        theme_minimal()
    
        #  Determine error metric
        error_metric <- if (error_type == "sd") "sd" else "sem"
        error_column <- sym(error_metric)
      
        #  Add error bars
        plot <- plot + geom_errorbar(
            aes(
                ymin = mean_percentage - !!error_column,
                ymax = mean_percentage + !!error_column
            ),
            position = position_dodge(width = 0.7),
            width = 0.25
        )
        
        # Store plot for the current substage
        plots[[substage]] <- plot
    }

  return(plots)
}


save_plots_to_pdf <- function(plots, substages, error_type, base_dir) {
  for (substage in substages) {
    plot_name <- paste("Fig-S7D_SCML2-in-nucleus_", substage, "_", error_type, ".pdf", sep = "")
    file_path <- file.path(base_dir, "plots", plot_name)
    
    # Retrieve the plot from the list using substage as key
    plot <- plots[[substage]]
    
    # Save the plot
    ggplot2::ggsave(file_path, plot = plot, width = 8, height = 6, dpi = 300)
  }
}


#  Load necessary libraries
library(tidyverse)


#  Read the TSV-file observational data into a data frame, and wrangle it a bit
file_path <- file.path(
    dir_script, "tsvs", "Fig-S7D_SCML2-in-nucleus_pachynema-diplonema.tsv"
)

data <- readr::read_tsv(file_path, show_col_types = FALSE) %>%
    #  No need to filter out all rows in bulk; we'll handle zeros during
    #+ percentage calculation
    dplyr::mutate(
        percentage_PE = ifelse(
            `S7D_observed_pachynema-early` > 0,
            `S7D_SCML2-in-nucleus_pachynema-early` /
            `S7D_observed_pachynema-early` * 100,
            NA
        ),
        percentage_PM = ifelse(
            `S7D_observed_pachynema-mid` > 0,
            `S7D_SCML2-in-nucleus_pachynema-mid` /
            `S7D_observed_pachynema-mid` * 100,
            NA
        ),
        percentage_PL = ifelse(
            `S7D_observed_pachynema-late` > 0,
            `S7D_SCML2-in-nucleus_pachynema-late` /
            `S7D_observed_pachynema-late` * 100,
            NA
        ),
        percentage_D = ifelse(
            `S7D_observed_diplonema` > 0,
            `S7D_SCML2-in-nucleus_diplonema` /
            `S7D_observed_diplonema` * 100,
            NA
        ),
        genotype = factor(genotype, levels = c("WT", "KO"))
    )

#  Group data by genotype and replicate, then summarize
grouped_data <- data %>%
    dplyr::group_by(genotype, replicate) %>%
    dplyr::summarize(
        mean_percentage_PE = mean(percentage_PE, na.rm = TRUE),
        mean_percentage_PM = mean(percentage_PM, na.rm = TRUE),
        mean_percentage_PL = mean(percentage_PL, na.rm = TRUE),
        mean_percentage_D = mean(percentage_D, na.rm = TRUE),
        .groups = 'drop'
    )

#  Perform substage-wise unpaired two-tailed t-tests between WT and KO. Assume
#+ 'grouped_data' already contains mean percentages for each substage as
#+ calculated previously.

# Define substages to loop through for t-tests
substages <- c("PE", "PM", "PL", "D")

# Initialize a list to store t-test results
t_test_results <- list()

# Loop through each substage to perform t-tests
for (substage in substages) {
    # substage <- "PE"
    
    #  Dynamically generate column names for mean percentages
    mean_col <- paste0("mean_percentage_", substage)
    
    # Extract mean percentages for WT and KO for the current substage
    wt_means <- grouped_data %>%
        dplyr::filter(genotype == "WT") %>%
        dplyr::pull(!!sym(mean_col))
    
    ko_means <- grouped_data %>%
        dplyr::filter(genotype == "KO") %>%
        dplyr::pull(!!sym(mean_col))
    
    # Check if there are enough observations to perform a t-test
    if (length(wt_means) > 1 && length(ko_means) > 1) {
        # Perform unpaired two-tailed t-test
        t_test_result <- t.test(wt_means, ko_means, var.equal = TRUE)
        
        # Store the p-value with the substage name as the key
        t_test_results[[substage]] <- t_test_result$p.value
    } else {
        # Store NA if not enough data for t-test
        t_test_results[[substage]] <- NA
    }
}

# Output the p-values from the t-tests
t_test_results
# $PE
# [1] 0.6015905
# 
# $PM
# [1] 0.08867762
# 
# $PL
# [1] 0.6666667
# 
# $D
# [1] 0.3164921

#  Summarize per-substage data across replicates for each genotype to calculate
#+ overall means, SDs, and SEMs
#+ 
#  Assuming 'data' has columns for percentage_PE, percentage_PM, percentage_PL,
#+ and percentage_D, and we want to summarize this data

#  Define substages based on the column names in 'data'
substages <- c("PE", "PM", "PL", "D")

# Initialize an empty list to store the summary data for each substage
summary_data <- list()

# Loop through each substage to summarize data
for (substage in substages) {
    # substage <- "PM"
    percentage_col <- paste0("percentage_", substage)  # Construct column name
  
  summary <- data %>%
    filter(!is.na(!!sym(percentage_col))) %>%  # Ensure valid observations
    group_by(genotype, substage = !!substage) %>%  # Group by genotype and dynamically add substage
    summarize(
      mean_percentage = mean(!!sym(percentage_col), na.rm = TRUE),
      sd = sd(!!sym(percentage_col), na.rm = TRUE),
      n = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      sem = sd / sqrt(n)
    )
  
  # Store the summary for each substage in the list
  summary_data[[substage]] <- summary
}

# Combine all substages' summary data into one dataframe
final_summary_data <- bind_rows(summary_data, .id = "substage")

# Correct the substage column to reflect the actual substages instead of the indices
final_summary_data$substage <- substages[final_summary_data$substage]

#+ Now, final_summary_data contains summarized data for each substage and
#+ genotype
#+ 
#+ Let's check if 'sd' and 'sem' are calculated correctly
print(final_summary_data)
# # A tibble: 8 × 6
#   genotype substage mean_percentage    sd     n   sem
#   <fct>    <chr>              <dbl> <dbl> <int> <dbl>
# 1 WT       PE                  80    44.7     5 20   
# 2 KO       PE                  73.5  43.7    17 10.6 
# 3 WT       PM                  77.8  44.1     9 14.7 
# 4 KO       PM                  92.5  24.5    20  5.47
# 5 WT       PL                  88.9  33.3     9 11.1 
# 6 KO       PL                 100     0       3  0   
# 7 WT       D                   75    45.2    12 13.1 
# 8 KO       D                   60    54.8     5 24.5 

#  Assuming 'final_summary_data' is the dataframe with summary stats for each
#+ substage and genotype and 'substages' is a vector of substages like
#+ c("PE", "PM", "PL", "D")
plots_sd <- plot_with_error_bars(final_summary_data, substages, "sd")
plots_sem <- plot_with_error_bars(final_summary_data, substages, "sem")

#  Optionally, print the plots to the console
print(plots_sd)
print(plots_sem)

#  Save the plots if needed, assuming 'dir_script' is the base directory
#+ and 'substages' is a vector of substages like c("PE", "PM", "PL", "D")

#  Save SD plots
save_plots_to_pdf(plots_sd, substages, "SD", dir_script)

#  Save SEM plots
save_plots_to_pdf(plots_sem, substages, "SEM", dir_script)

data %>%
    dplyr::group_by(genotype) %>%
    dplyr::summarize(
        n_PM = sum(`S7D_observed_pachynema-mid`),
        n_PL = sum(`S7D_observed_pachynema-late`),
        n_D =  sum(`S7D_observed_diplonema`)
    )
# # A tibble: 2 × 4
#   genotype  n_PM  n_PL   n_D
#   <fct>    <dbl> <dbl> <dbl>
# 1 WT          24    21    25
# 2 KO         107     5    10
