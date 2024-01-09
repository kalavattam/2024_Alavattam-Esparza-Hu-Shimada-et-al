#!/usr/bin/env Rscript

#  Fig-3B.R
#  KA

#  Load necessary libraries
library(argparse)
library(googlesheets4)
library(tidyverse)
library(scales)

#  Run authentication workflow for Google Sheets API using the googlesheets4
#+ package
run_gs4_auth <- FALSE  # Only needs to be run the first time
if (isTRUE(run_gs4_auth)) googlesheets4::gs4_auth()

#  Read in Google Sheets sheet "yH2AX_all"
URL_pre <- "https://docs.google.com/spreadsheets/d"
URL_suf <- "1mjOWHt3MS4rwYOlr0aiIeIeHrIBvzUSPBd6GnjuXES4/edit"
URL <- paste(URL_pre, URL_suf, sep = "/")
data <- googlesheets4::read_sheet(URL)

#  Further wrangle the data, calculating total sample sizes along with the
#+ summary statistics SD and SEM
aggregated_data <- data %>%
    dplyr::filter(class == 1) %>%
    dplyr::select(c(factor, substage, genotype, replicate, sample_size)) %>%
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
    aggregated_data, aes(x = substage_genotype, y = mean_size, fill = genotype)
) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(
        aes(ymin = mean_size - sd, ymax = mean_size + sd),
        width = 0.25,
        position = position_dodge(0.9)
    ) +
    scale_y_continuous(labels = scales::percent) +
    labs(y = "Mean proportion (%)", x = "Substage and genotype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

#  Plot with SEM bars
plot_sem <- ggplot2::ggplot(
    aggregated_data, aes(x = substage_genotype, y = mean_size, fill = genotype)
) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(
        aes(ymin = mean_size - sem, ymax = mean_size + sem),
        width = 0.25,
        position = position_dodge(0.9)
    ) +
    scale_y_continuous(labels = scales::percent) +
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
    #  Save the plot with SD bars (will be saved to $HOME directory)
    ggplot2::ggsave(
        "Fig-3B_prop-stages_SD.pdf",
        plot = plot_sd,
        width = 8.5,
        height = 6,
        dpi = 300
    )
    
    #  Save the plot with SEM bars (will be saved to $HOME directory)
    ggplot2::ggsave(
        "Fig-3B_prop-stages_SEM.pdf",
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
        dplyr::pull(sample_size)  # Extract the sample_size column
  
    data_ko <- data %>%
        dplyr::filter(class == 1 & substage == sub & genotype == "KO") %>%
        dplyr::pull(sample_size)  # Extract the sample_size column
  
  #  Perform the t-test
  t_test <- t.test(data_wt, data_ko, var.equal = TRUE)

  #  Store the result
  t_test_results[[sub]] <- t_test
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
#+ - ZEM  [1] 1
#+ - ZL   [1] 0.8874654
#+ - PEM  [1] 0.8507149
#+ - PL   [1] 0.1589297
#+ - D    [1] 0.03323569
