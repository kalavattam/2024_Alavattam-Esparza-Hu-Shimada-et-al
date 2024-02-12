#!/usr/bin/env Rscript

#  KA_Fig-3C.R
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
#  Function to plot stacked bar graphs for WT and KO genotypes with percentages
#+ of counts assigned to each category
plot_bar_graphs_categorical <- function(aggregated_data) {
    # aggregated_data <- agg_data
    
    #  Create a new column that combines substage and genotype with the desired
    #+ order; cast it as a factor
    aggregated_data$substage_genotype <- factor(
        paste(aggregated_data$substage, aggregated_data$genotype),
        levels = c(
            "ZEM WT", "ZEM KO", "ZL WT", "ZL KO", "PEM WT", "PEM KO",
            "PL WT", "PL KO", "D WT", "D KO"
        )
    )
    
    #  Calculate the total counts for each substage_genotype combination
    total_counts_per_substage_genotype <- aggregated_data %>%
        group_by(substage_genotype) %>%
        summarize(total_counts = sum(total_tally), .groups = 'drop')
    
    #  Then, join this back to the original data to calculate percentages
    aggregated_data_w_pcts <- dplyr::left_join(
        aggregated_data,
        total_counts_per_substage_genotype,
        by = c("substage_genotype")
    ) %>%
        mutate(percentage = (total_tally / total_counts) * 100)
    
    #  Plotting the data with corrected percentages
    plot <- ggplot2::ggplot(
        aggregated_data_w_pcts,
        aes(x = substage_genotype, y = percentage, fill = as.factor(class))
    ) +
        geom_bar(stat = "identity", position = "fill") +
        scale_y_continuous(
            labels = scales::percent_format()
        ) +
        labs(
            x = "Genotype",
            y = "Percentage", 
            title = "Percentage of counts by class for each genotype in EP",
            fill = "Class"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.title = element_text(size = 12),
            legend.text = element_text(size = 10)
        )
    
    return(plot)
}

#  Load in and slightly wrangle the data from the TSV file
data <- readr::read_tsv(
    file.path(
        dir_script, "tsvs",
        "wrangle-assess_blinded-scoring_Atf7ip2.yH2AX_KA-eval_work-unblinded_sort-practical_etc.tsv"
    ),
    show_col_types = FALSE,
    na = c("", "NA", "#N/A")
) %>%
    dplyr::rename(substage = eval_stage_4, class = eval_class) %>%
    dplyr::select(c(genotype, replicate, substage, class)) %>%
    dplyr::mutate(
        substage_genotype = factor(
            paste(substage, genotype),
            levels = c(
                "ZEM WT", "ZEM KO", "ZL WT", "ZL KO", "PEM WT", "PEM KO",
                "PL WT", "PL KO", "D WT", "D KO"
            )
        ),
        class = factor(class)
    ) %>%
    dplyr::relocate(c(substage_genotype, substage), .before = genotype) %>%
    na.omit()

#  Tally counts for each combination in the original data
data_count <- data %>%
    dplyr::group_by(substage_genotype, replicate, class) %>%
    dplyr::summarize(count = n(), .groups = 'drop')


#  Get all unique combinations of substage_genotype, replicate, and class
comp_combos <- expand.grid(
    substage_genotype = unique(data$substage_genotype),
    replicate = unique(data$replicate),
    class = unique(data$class)
)

#  Join comp_combos with data_count, filling in missing counts with 0
props <- comp_combos %>%
    dplyr::left_join(data_count, by = c("substage_genotype", "replicate", "class")) %>%
    tidyr::replace_na(list(count = 0)) %>%
    dplyr::arrange(substage_genotype, class)

#  Calculate the total counts and proportions for each substage_genotype and
#+ replicate combination
props <- props %>%
    dplyr::group_by(substage_genotype, replicate) %>%
    dplyr::mutate(total = sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(proportion = ifelse(total == 0, 0, count / total))

#  Calculate mean proportions across replicates
agg_props <- props %>%
    dplyr::group_by(substage_genotype, class) %>%
    dplyr::summarize(
        mean_proportion = mean(proportion),
        sd = sd(proportion),
        sem = sd / sqrt(n()),
        .groups = 'drop'
    ) %>%
    dplyr::group_by(substage_genotype) %>%
    dplyr::arrange(substage_genotype, dplyr::desc(class)) %>%
    dplyr::mutate(
        cumulative_proportion = cumsum(mean_proportion)
    ) %>%
    dplyr::ungroup()

#  Plot with SD bars
plot_sd <- ggplot2::ggplot(
    agg_props,
    aes(x = substage_genotype, y = mean_proportion, fill = class)
) +
    geom_bar(stat = "identity") +
    geom_errorbar(
        aes(
            ymin = cumulative_proportion - sd,
            ymax = cumulative_proportion + sd
        ),
        width = 0.25,
        position = "identity"
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Proportion (%)", x = "Substage and genotype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

#  Plot with SEM bars
plot_sem <- ggplot2::ggplot(
    agg_props,
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
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Proportion (%)", x = "Substage and genotype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

print_plots <- TRUE
if (isTRUE(print_plots)) {
    print(plot_sem)
    print(plot_sd)
}

save_plots <- FALSE
if (isTRUE(save_plots)) {
    #  Save the plot with SD bars (will be saved to $HOME directory)
    ggplot2::ggsave(
        file.path(dir_script, "plots", "KA-Fig-3C_eval-stage-4_SD.pdf"),
        plot = plot_sd,
        width = 8.5,
        height = 6,
        dpi = 300
    )
    
    #  Save the plot with SEM bars (will be saved to $HOME directory)
    ggplot2::ggsave(
        file.path(dir_script, "plots", "KA-Fig-3C_eval-stage-4_SEM.pdf"),
        plot = plot_sem,
        width = 8.5,
        height = 6,
        dpi = 300
    )
}

agg_data <- data %>%
    dplyr::group_by(substage, genotype, class) %>%
    dplyr::summarize(total_tally = n(), .groups = 'drop')

#  Plot the categorical data (no consideration of per-replicate means)
plot_categorical <- plot_bar_graphs_categorical(agg_data)

print_plot <- TRUE
if (isTRUE(print_plot)) {
    print(plot_categorical)
}

save_plot <- TRUE
if (isTRUE(save_plot)) {
    ggplot2::ggsave(
        file.path(
            dir_script, "plots", "KA-Fig-3C_eval-stage-4_categorical.pdf"
        ),
        plot = plot_categorical,
        width = 8.5,
        height = 6,
        dpi = 300
    )
}

#  Initialize a list to store contingency tables
contingency_tables <- list()

#  Iterate over each substage, creating contingency tables for each
for (sub in unique(agg_data$substage)) {
    # sub <- "PEM"
    
    #  For each sub, create a dataframe (complete_combinations_sub) that
    #+ contains all combinations of class and genotype
    classes <- unique(agg_data$class)
    genotypes <- unique(agg_data$genotype)
    complete_combinations_sub <- expand.grid(class = classes, genotype = genotypes)
    
    #  Subset the aggregated data for the current substage
    substage_data <- subset(agg_data, substage == sub)
    
    #  Merge substage_data with complete_combinations_sub, replacing NA with 0
    substage_data_comp <- merge(
        complete_combinations_sub, substage_data, by = c("class", "genotype"),
        all.x = TRUE
    )
    substage_data_comp[is.na(substage_data_comp)] <- 0
    
    #  Filter out rows where all tally values are zero: This is done by
    #+ reshaping the data to a wide format and then filtering to exclude rows
    #+ with row sums of zero 
    substage_data_wide <- reshape2::dcast(
        substage_data_comp, class ~ genotype, value.var = "total_tally"
    )
    substage_data_wide <- substage_data_wide[
        rowSums(substage_data_wide[, -1]) != 0, 
    ]

    #  Reshape the filtered data back to a long format for creating the
    #+ contingency table
    substage_data_filt <- reshape2::melt(
        substage_data_wide,
        id.vars = "class",
        variable.name = "genotype",
        value.name = "total_tally"
    )
    
    #  Drop unused factor (categorical variable) levels
    substage_data_filt <- droplevels(substage_data_filt)
    
    #  Reshape the data to create a 2x3 contingency table
    contingency_table <- stats::xtabs(
        total_tally ~ class + genotype, data = substage_data_filt
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
if (isTRUE(check_p_values)) {
    print(fisher_results$ZEM$p.value)
    print(fisher_results$ZL$p.value)
    print(fisher_results$PEM$p.value)
    print(fisher_results$PL$p.value)
    print(fisher_results$D$p.value)
}

#  The p-values from Fisher's exact tests are as follows:
#+ - ZEM 1
#+ - ZL  1
#+ - PEM 6.896717e-33
#+ - PL  0.156128
#+ - D   2.675916e-09

#  Check p-values from Chi-squared tests
check_p_values <- FALSE
if (isTRUE(check_p_values)) {
    print(chisq_results$ZEM)
    print(chisq_results$ZL)
    print(chisq_results$PEM)
    print(chisq_results$PL)
    print(chisq_results$D)
}

#  The p-values from Chi-squared tests are as follows:
#+ - ZEM NA
#+ - ZL  NA
#+ - PEM NA
#+ - PL  NA
#+ - D   NA
