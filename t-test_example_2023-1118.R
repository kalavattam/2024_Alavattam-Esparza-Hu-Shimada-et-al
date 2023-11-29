#!/usr/bin/env Rscript

#  t-test_example_2023-1118.R
#  KA

#  The following code draws on data in tab "K9me3_MP (Kris)" in file
#+ "Mcaf2_Revisions_BlindedQuant_2023-1116.KA-2023-1118.xlsx"

library(dplyr)


#  Perform t-test ------------------------------------------------------------
#  Create data frames for data points for ctrl and KO groups
ctrl_data <- data.frame(
    value = c(100.00, 91.67, 75.00),
    group = "ctrl"
)

ko_data <- data.frame(
    value = c(0.00, 39.13, 8.33, 100.00),
    group = "ko"
)

#  Combine the data into one data frame
combined_data <- rbind(ctrl_data, ko_data)

#  Perform the t-test
t_test_result <- t.test(
    value ~ group,
    data = combined_data,
    var.equal = TRUE  # Student's t-test
)

#  Prints the t-test results
print(t_test_result)


#  Calculate summary statistics per group -------------------------------------
#  Create data frame for per-group summary statistics
combined_stats <- combined_data %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(
        mean = mean(value),
        SD = sd(value),
        SEM = sd(value) / sqrt(n())
    )

#  Print the summary statistics
print(combined_stats)


#  Run Fisher's exact test ----------------------------------------------------
#  Set up contingency table
contingency <- data.frame(
    ctrl = c(86, 16),
    ko = c(28, 41),
    row.names = c("yes", "no")
)

#  Perform the Fisher's exact test
fisher_test_results <- fisher.test(contingency)

#  Print the Fisher's exact test results
print(fisher_test_results)
