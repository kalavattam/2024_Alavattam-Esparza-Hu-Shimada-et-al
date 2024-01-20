#!/usr/bin/env Rscript

#  Fig-4B.R
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

#  Load in data from TSV file; files of interest are stored in the base
#+ directory of the repo, which is the same location
data <- readr::read_tsv(
    file.path(
        dir_script, "tsvs",
        "wrangle-assess_blinded-scoring_Atf7ip2.K9me3_EP_observations.tsv"
    ),
    show_col_types = FALSE
) %>%
    
    #  Convert the proportions to numeric format (if it's not already)
    dplyr::mutate(proportion = as.numeric(gsub("%", "", proportion)) / 100) %>%

    #  Combine substage and genotype into a new factor with a specific order
    dplyr::mutate(
        substage_genotype = factor(
            paste(substage, genotype), levels = c("EP WT", "EP KO")
        ),
        class = factor(class)
    )

#  Aggregate data to calculate mean proportion and other summary statistics for
#+ each combination of substage_genotype and class
aggregated_data <- data %>%
    #  Group the data by both substage_genotype and class, thereby creating
    #+ subsets of data for each combination of substage_genotype and class
    dplyr::group_by(substage_genotype, class) %>%
    
    #  Summarize the data in each group
    #+ - mean_proportion: Calculate the average of proportions, ignoring NA
    #+                    values
    #+ - sd: Calculate the standard deviation of proportions, ignoring NA
    #+       values
    #+ - n: Count the number of observations in each group
    #+ - .groups = 'drop': Prevents the creation of an extra layer of grouping,
    #+                     which would otherwise be carried over to subsequent
    #+                     operations
    dplyr::summarize(
        mean_proportion = mean(proportion, na.rm = TRUE),
        sd = sd(proportion, na.rm = TRUE),
        n = n(),
        .groups = 'drop'
    ) %>%
    
    #  Create new columns based on existing ones.
    #+ - class: Convert the class column into a factor (categorical variable)
    #+ - sem: Calculate the standard error of the mean (SEM); SEM is the
    #+        standard deviation (sd) divided by the square root of the sample
    #+        size (n)
    dplyr::mutate(
        class = factor(class),
        sem = sd / sqrt(n)
    ) %>%
    
    #  Regroup the data by substage_genotype only; this is necessary for the
    #+ cumulative summation operations performed below, as it limits the
    #+ summation operations to individual groups rather than across groups
    dplyr::group_by(substage_genotype) %>%
    
    #  Arrange the data within each substage_genotype group in descending order
    #+ of class; this affects how error bars will be stacked in the plot
    dplyr::arrange(dplyr::desc(class)) %>%
    
    #  Calculate the cumulative sum of mean_proportion within each
    #+ substage_genotype group; this is used for positioning the stacked bars
    #+ and their corresponding error bars
    dplyr::mutate(cumulative_proportion = cumsum(mean_proportion)) %>%
    
    #  Ungroup the data to remove any residual grouping, ensuring that
    #+ subsequent operations are performed on the entire dataset
    dplyr::ungroup()

#  Plot with SD bars:
#  Initialize the plot with the ggplot function, specifying the data source and
#+ aesthetic mappings
#+ 
#+ The aes() function defines the aesthetics (how data are represented in the
#+ plot):
#+ - x = substage_genotype: Defines the x-axis using the 'substage_genotype'
#+                          column
#+ - y = mean_proportion: Sets the y-axis to represent the mean proportion
#+ - fill = class: Specifies that the color fill of the bars should represent
#+                 different 'class' values
plot_sd <- ggplot2::ggplot(
    aggregated_data,
    aes(x = substage_genotype, y = mean_proportion, fill = class)
) +
    
    #  Add bar geometry to the plot:
    #+ stat = "identity" tells ggplot2 that the y-values provided
    #+ (mean_proportion) are the heights of the bars
    geom_bar(stat = "identity") +
    
    #  Add error bars to the plot:
    #+ ymin and ymax define the lower and upper limits of the error bars,
    #+ representing the range of SD
    #+ - cumulative_proportion - sd gives the lower end of the error bar
    #+ - cumulative_proportion + sd gives the upper end of the error bar
    #+ - width = 0.25 sets the width of the horizontal lines at the top and
    #+   bottom of the error bars
    #+ - position = "identity" ensures error bars are drawn exactly at the
    #+   specified y-values without adjustment
    geom_errorbar(
        aes(
            ymin = cumulative_proportion - sd,
            ymax = cumulative_proportion + sd
        ),
        width = 0.25,
        position = "identity"
    ) +
    
    #  Scale the y-axis to display values as percentages:
    #+ labels = scales::percent formats the y-axis labels as percentages
    scale_y_continuous(labels = scales::percent) +
    
    #  Add labels to the axes:
    #+ - y = "Proportion (%)": Label for the y-axis
    #+ - x = "Substage and genotype": Label for the x-axis
    labs(y = "Proportion (%)", x = "Substage and genotype") +
    
    #  Apply a minimal theme to the plot
    theme_minimal() +
    
    # Customize the theme, specifically the x-axis text:
    #+ axis.text.x = element_text(angle = 90, hjust = 1):
    #+ - angle = 90 rotates the x-axis labels by 90 degrees for better
    #+   readability
    #+ - hjust = 1 "right" (top)-justifies the text labels at the x-axis ticks
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

#  Plot with SEM bars
plot_sem <- ggplot2::ggplot(
    aggregated_data,
    aes(x = substage_genotype, y = mean_proportion, fill = class)
) +
    geom_bar(stat = "identity") +
    geom_errorbar(
        aes(
            ymin = cumulative_proportion - sem,
            ymax = cumulative_proportion + sem
        ),
        width = 0.25,
        position = "identity"
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1.2)) +
    labs(y = "Proportion (%)", x = "Substage and genotype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

print_plots <- TRUE
if (isTRUE(print_plots)) {
    print(plot_sem)
    print(plot_sd)
}

save_plots <- TRUE
if (isTRUE(save_plots)) {
    #  Save the plot with SD bars
    ggplot2::ggsave(
        file.path(dir_script, "plots", "Fig-4B_prop-K9me3_EP-classes_SD.pdf"),
        plot = plot_sd,
        width = 8.5,
        height = 6,
        dpi = 300
    )
    
    #  Save the plot with SEM bars
    ggplot2::ggsave(
        file.path(dir_script, "plots", "Fig-4B_prop-K9me3_EP-classes_SEM.pdf"),
        plot = plot_sem,
        width = 8.5,
        height = 6,
        dpi = 300
    )
}

#  Run Fisher's exact and Chi-squared tests
#  Aggregate data to get the total tally for each substage, genotype, and class
aggregated_data <- data %>%
    dplyr::group_by(substage, genotype, class) %>%
    dplyr::summarize(total_tally = sum(total, na.rm = TRUE), .groups = 'drop')

#  Initialize a list to store contingency tables
contingency_tables <- list()

#  Iterate over each substage, creating contingency tables for each
for (sub in unique(aggregated_data$substage)) {
    #  Subset the aggregated data for the current substage
    substage_data <- subset(aggregated_data, substage == sub)
    
    #  Filter out rows where all tally values are zero: This is done by
    #+ reshaping the data to a wide format and then filtering to exclude rows
    #+ with row sums of zero 
    substage_data_wide <- reshape2::dcast(
        substage_data, class ~ genotype, value.var = "total_tally"
    )
    substage_data_wide <- substage_data_wide[
        rowSums(substage_data_wide[, -1]) != 0, 
    ]

    #  Reshape the filtered data back to a long format for creating the
    #+ contingency table
    substage_data_filtered <- reshape2::melt(
        substage_data_wide,
        id.vars = "class",
        variable.name = "genotype",
        value.name = "total_tally"
    )
    
    #  Drop unused factor (categorical variable) levels
    substage_data_filtered <- droplevels(substage_data_filtered)
    
    #  Reshape the data to create a 2x4 contingency table
    contingency_table <- stats::xtabs(
        total_tally ~ class + genotype, data = substage_data_filtered
    )

    #  Store the table using the substage as the name
    contingency_tables[[sub]] <- contingency_table
}

#  Initialize lists to store results
fisher_results <- list()  # Initialize a list to store Fisher's test results
chisq_results <- list()   # Initialize a list to store Chi-squared test results

for (sub in names(contingency_tables)) {
    #  Perform Fisher's exact test; include warning/error capture for the sake
    #+ of completeness
    fisher_test_result <- tryCatch({
        fisher.test(contingency_tables[[sub]])
    }, warning = function(w) {
        return(NA)  # Return NA if there's a warning
    }, error = function(e) {
        return(NA)  # Return NA if there's an error
    })

    #  Perform Chi-squared test with a try-catch to handle potential
    #+ errors/warnings
    chisq_test_result <- tryCatch({
        chisq.test(contingency_tables[[sub]])
    }, warning = function(w) {
        return(NA)  # Return NA if there's a warning (e.g., expected count low)
    }, error = function(e) {
        return(NA)  # Return NA if there's an error
    })

    #  Store the test results
    fisher_results[[sub]] <- fisher_test_result
    chisq_results[[sub]] <- chisq_test_result
}

#  Check p-values from Fisher's exact tests
check_p_values <- TRUE
if (isTRUE(check_p_values)) print(fisher_results$EP$p.value)

#  The p-values from Fisher's exact tests are as follows:
#+ - EP 1.795446e-07

#  Check p-values from Chi-squared tests
check_p_values <- TRUE
if (isTRUE(check_p_values)) print(chisq_results$EP)

#  The p-values from Chi-squared tests are as follows:
#+ - EP 1.243e-06
