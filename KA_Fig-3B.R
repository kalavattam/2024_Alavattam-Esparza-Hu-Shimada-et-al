#!/usr/bin/env Rscript

#  KA_Fig-3B.R
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
#+ directory of the repo, which is the same location as this script
data <- readr::read_tsv(
    file.path(
        dir_script, "tsvs",
        "wrangle-assess_blinded-scoring_Atf7ip2.yH2AX_KA-eval_work-unblinded_sort-practical_etc.tsv"
    ),
    show_col_types = FALSE,
    na = c("", "NA", "#N/A")
) %>%
    dplyr::select(c(genotype, replicate, eval_stage_4, eval_class)) %>%
    na.omit()

#  Aggregate data, calculating sample sizes and percent proportions
agg_data <- data %>%
    dplyr::group_by(
        substage = eval_stage_4, genotype, replicate, class = eval_class
    ) %>%
    dplyr::summarize(tally = n(), .groups = 'drop') %>%
    dplyr::group_by(substage, genotype, replicate) %>%
    dplyr::mutate(sample_size = sum(tally)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(proportion = scales::percent(tally / sample_size))

#  Further wrangle the data, calculating total sample sizes along with the
#+ summary statistics SD and SEM
agg_data_3B <- agg_data %>%
    dplyr::distinct(
        substage, genotype, replicate, sample_size, .keep_all = TRUE
    ) %>%
    dplyr::mutate(substage_genotype = paste(substage, genotype)) %>%
    dplyr::group_by(substage_genotype, substage, genotype) %>%
    dplyr::summarize(
        mean_size = mean(sample_size) / 100,
        sd = sd(sample_size) / 100,
        sem = sd / sqrt(n()),
        .groups = 'drop'
    ) %>%
    dplyr::mutate(
        substage_genotype = factor(
            substage_genotype,
            levels = c(
                "ZEM WT", "ZEM KO", "ZL WT", "ZL KO", "PEM WT", "PEM KO",
                "PL WT", "PL KO", "D WT", "D KO"
            )
        )
    )

#  Plot with SD bars
plot_sd <- ggplot2::ggplot(
    agg_data_3B, aes(x = substage_genotype, y = mean_size, fill = genotype)
) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(
        aes(ymin = mean_size - sd, ymax = mean_size + sd),
        width = 0.25,
        position = position_dodge(0.9)
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(y = "Mean proportion (%)", x = "Substage and genotype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

#  Plot with SEM bars
plot_sem <- ggplot2::ggplot(
    agg_data_3B, aes(x = substage_genotype, y = mean_size, fill = genotype)
) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(
        aes(ymin = mean_size - sem, ymax = mean_size + sem),
        width = 0.25,
        position = position_dodge(0.9)
    ) +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    labs(y = "Mean proportion (%)", x = "Substage and genotype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

print_plots <- TRUE
if (isTRUE(print_plots)) {
    print(plot_sd)
    print(plot_sem)
}

save_plots <- TRUE
if (isTRUE(save_plots)) {
    #  Save the plot with SD bars
    ggplot2::ggsave(
        file.path(dir_script, "plots", "KA-Fig-3B_eval-stage-4_SD.pdf"),
        plot = plot_sd,
        width = 8.5,
        height = 6,
        dpi = 300
    )
    
    #  Save the plot with SEM bars
    ggplot2::ggsave(
        file.path(dir_script, "plots", "KA-Fig-3B_eval-stage-4_SEM.pdf"),
        plot = plot_sem,
        width = 8.5,
        height = 6,
        dpi = 300
    )
}

#  Run t-tests for each substage
substage_levels <- unique(agg_data$substage)  # Get unique substages
t_test_results <- list()  # Initialize a list to store the results

#  Perform the t-tests with stage-wise proportion summary statistics
for (sub in substage_levels) {
    #  Initialize variables for WT summary statistics
    WT_mean <- agg_data_3B %>%
        dplyr::filter(substage == sub & genotype == "WT") %>%
        dplyr::pull(mean_size)
    
    WT_SEM <- agg_data_3B %>%
        dplyr::filter(substage == sub & genotype == "WT") %>%
        dplyr::pull(sem)
    
    WT_n <- agg_data %>% 
        dplyr::filter(substage == sub & genotype == "WT") %>%
        dplyr::distinct(substage, genotype, replicate, sample_size) %>%
        dplyr::pull(sample_size) %>%
        sum()
    
    #  Initialize variables for KO summary statistics
    KO_mean <- agg_data_3B %>%
        dplyr::filter(substage == sub & genotype == "KO") %>%
        dplyr::pull(mean_size)
    
    KO_SEM <- agg_data_3B %>%
        dplyr::filter(substage == sub & genotype == "KO") %>%
        dplyr::pull(sem)
    
    KO_n <- agg_data %>% 
        dplyr::filter(substage == sub & genotype == "KO") %>%
        dplyr::distinct(substage, genotype, replicate, sample_size) %>%
        dplyr::pull(sample_size) %>%
        sum()
    
    #  Calculate the t-statistic
    t_value <- (WT_mean - KO_mean) / sqrt(WT_SEM^2 + KO_SEM^2)
    
    #  Calculate the degrees of freedom
    df <- ((WT_SEM^2 + KO_SEM^2)^2) / ((WT_SEM^4 / (WT_n - 1)) + (KO_SEM^4 / (KO_n - 1)))
    
    #  Calculate the t-test p-value
    p_value <- 2 * pt(-abs(t_value), df)

    #  Store the result
    t_test_results[[sub]][["p.value"]] <- p_value
}

#  Assess p-values
check_p_values <- TRUE
if (isTRUE(check_p_values)) {
    print(t_test_results$ZEM$p.value)
    print(t_test_results$ZL$p.value)
    print(t_test_results$PEM$p.value)
    print(t_test_results$PL$p.value)
    print(t_test_results$D$p.value)
}

#  The p-values from t-tests are as follows:
#+ - ZEM  0.8104641
#+ - ZL   0.03361754
#+ - PEM  0.4533276
#+ - PL   0.03098515
#+ - D    0.04498352
