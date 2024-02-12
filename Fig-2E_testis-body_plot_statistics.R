
#  Fig-2D_testis-body_plot_statistics.R
#  KA
#  2020-0811-0812, 2022-0803, 2024-0210

#  Measure and visualize testis and body masses at groups of weeks as time-
#+ point categories


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


#  Load in libraries
library("tidyverse")
library("car")


#  Define functions
sem <- function(x) {
    # Calculate the standard error of the mean (SEM) for a vector of numbers
    # :param x: vector of floats <dbl>
    sd(x) / sqrt(length(x))
}


week_month_summarize <- function(tbl, u, v, t) {
    # Function will group weeks and genotype, then generate summary means and
    # SDs, and SEMs for each group
    # :param tbl: a tibble
    # :param u: a non-quoted variable/column within the tibble for genotype
    #           information
    # :param v: a non-quoted variable/column within the tibble, e.g., testis_mg
    # :param t: a non-quoted category of time, weeks or months, for grouping
    #           and summarization
    # 
    # stackoverflow.com/questions/48062213/dplyr-using-column-names-as-function-arguments
    enq_u <- enquo(u)
    enq_v <- enquo(v)
    enq_t <- enquo(t)
    tbl %>%
        dplyr::group_by(!!enq_u, !!enq_t) %>% 
        dplyr::summarize(
            mean = mean(!!enq_v),
            sd = sd(!!enq_v),
            sem = sem(!!enq_v)
        )
}


plot_line3 <- function(x, t, error = "sd") {
    # Function will plot the output of week_summarize() and is set up for
    # analyses of three genotype categories: WT, HZ, and KO
    # :param x: an object containing the output of week_summarize()
    # :param t: a non-quoted category of time, weeks or months, for plotting
    # :param error: type of error to display, "sem" for standard error or "sd"
    #               for standard deviation
    # 
    # stackoverflow.com/questions/57380726/how-to-use-quosure-in-ggplot-to-generate-plotting-function
    enq_t <- enquo(t)
    
    error_width <- if (error == "sd") {
        x$sd
    } else if (error == "sem") {
        x$sem
    } else {
        stop("Error type not supported. Use 'sd' for standard deviation or 'sem' for standard error of the mean.")
    }
    
    x %>%
        ggplot2::ggplot(ggplot2::aes(!!enq_t, mean, group = genotype)) +
        ggplot2::geom_line(aes(linetype = genotype)) +
        ggplot2::geom_point(size = 2, aes(shape = genotype)) +
        ggplot2::geom_errorbar(
            ggplot2::aes(ymin = mean - error_width, ymax = mean + error_width),
            width = 0.05
        ) +
        ggplot2::scale_linetype_manual(
            values = c("longdash", "dotdash", "dotted")
        ) +
        ggplot2::scale_shape_manual(values = c(16, 1, 17)) +
        ggplot2::theme(legend.title = element_blank()) +
        ggplot2::labs(x = enq_t, fill = "Genotype")
}


annotate_stars <- function(vector_p_values) {
    ifelse(
        vector_p_values > 0.05,
        "n.s.",
        ifelse(
            vector_p_values <= 0.001,
            "***",
            ifelse(
                vector_p_values <= 0.01,
                "**",
                ifelse(
                    vector_p_values <= 0.05,
                    "*",
                    NA
                )
            )
        )
    )
}


#  Load data
tsv <- "Atf7ip2_testis-body.txt"
tbl <- read.delim(
    file.path(dir_script, "tsvs", tsv),
    header = TRUE
) %>%
    tibble::as_tibble()
tbl$sac_date <- lubridate::mdy(tbl$sac_date)
tbl$birth_date <- lubridate::mdy(tbl$birth_date)

#  Add columns for testis weights (mg) and testis-to-body-weight ratios
tbl <- tbl %>%
    dplyr::mutate(
        testis_1_mg = testis_1_g * 1000,
        testis_2_mg = testis_2_g * 1000,
        testis_combined_mg = testis_combined_g * 1000,
        t_to_b_1 = testis_1_mg / body_g,
        t_to_b_2 = testis_2_mg / body_g,
        t_to_b_c = testis_combined_mg / body_g
    )

#  Weeks of age (rounded)
tbl$weeks <- ((tbl$sac_date - tbl$birth_date) / 7) %>%
    round() %>%
    as.factor()

levels(tbl$weeks) <- c(
    "4–6", "4–6", "4–6",
    "7–8", "7–8",
    "9–12", "9–12",
    "13–26", "13–26", "13–26",
    "27–52", "27–52", "27–52", "27–52", "27–52",
    "53+", "53+", "53+", "53+"
)


#  Organize the data for distinctions, adjustments, etc.

#  Pivot/gather the testis_#_mg data to remove distinctions between first and
#+ second testes
tbl_mg <- tbl %>%
    tidyr::pivot_longer(
        cols = testis_combined_mg,
        names_to = "testis_meta",
        values_to = "testis_mg"
    )


#  Pivot/gather the t_to_b_# ratio data to remove distinctions between first
#+ and second testes
tbl_rat <- tbl %>%
    tidyr::pivot_longer(
        cols = t_to_b_c,
        names_to = "t_to_b_meta",
        values_to = "t_to_b"
    )


#  Run analyses: Three genotype categories: WT, HZ, and KO --------------------
cats_init <- c("HZ", "KO", "WT")
cats_ord <- c("WT", "HZ", "KO")

tbl$genotype <- factor(tbl$genotype)
levels(tbl$genotype) <- cats_init
tbl$genotype <- tbl$genotype %>% ordered(levels = cats_ord)

tbl_mg$genotype <- factor(tbl_mg$genotype)
levels(tbl_mg$genotype) <- cats_init
tbl_mg$genotype <- tbl_mg$genotype %>% ordered(levels = cats_ord)

tbl_rat$genotype <- factor(tbl_rat$genotype)
levels(tbl_rat$genotype) <- cats_init
tbl_rat$genotype <- tbl_rat$genotype %>% ordered(levels = cats_ord)

tbl$genotype.simple <- ifelse(
    tbl$genotype == "WT" | tbl$genotype == "HZ", "Ctrl", "KO"
) %>%
    factor()

tbl_mg$genotype.simple <- ifelse(
    tbl_mg$genotype == "WT" | tbl_mg$genotype == "HZ", "Ctrl", "KO"
) %>%
    factor()

tbl_rat$genotype.simple <- ifelse(
    tbl_rat$genotype == "WT" | tbl_rat$genotype == "HZ", "Ctrl", "KO"
) %>%
    factor()

# #  Summarize/examine testis weights (mg) by weeks
# s <- week_month_summarize(tbl_mg, testis_mg, weeks)
# p <- plot_line3(s, weeks, "sd")
# p + labs(y = "Testis weights (mg)")

# #  Summarize/examine body weights (g) by weeks
# s <- week_month_summarize(tbl, body_g, weeks)
# p <- plot_line3(s, weeks, "sd")
# p + labs(y = "Body weights (g)")

#  Summarize/examine ratios of testis weights (mg) to body weights (g) by weeks
s <- week_month_summarize(tbl_rat, genotype, t_to_b, weeks)
p <- plot_line3(s, weeks, "sd")
(
    q <- p + ggplot2::labs(
        y = "Ratios of combined testis weights (mg) to body weights (mg)"
    )
)

s <- week_month_summarize(tbl_rat, genotype.simple, t_to_b, weeks) %>%
    dplyr::rename("genotype" = "genotype.simple")
p <- plot_line3(s, weeks, "sd")
(
    q <- p + ggplot2::labs(
        y = "Ratios of combined testis weights (mg) to body weights (mg)"
    )
)

table_weeks_no <- tbl %>%  #MAYBE Write out
    dplyr::group_by(weeks) %>%
    dplyr::summarize(
        no_samples = length(weeks)
    )
table_genotype_no <- tbl %>%  #MAYBE Write out
    dplyr::group_by(genotype) %>%
    dplyr::summarize(
        no_samples = length(genotype)
    )
table_weeks_genotype_no <- tbl %>%  #MAYBE Write out
    dplyr::count(weeks, genotype, name = "no_samples")

table_weeks_no_s <- tbl %>%  #MAYBE Write out
    dplyr::group_by(weeks) %>%
    dplyr::summarize(
        no_samples = length(weeks)
    )
table_genotype_no_s <- tbl %>%  #MAYBE Write out
    dplyr::group_by(genotype.simple) %>%
    dplyr::summarize(
        no_samples = length(genotype)
    )
table_weeks_genotype_no_s <- tbl %>%  #MAYBE Write out
    dplyr::count(weeks, genotype.simple, name = "no_samples")


#  Run statistical analyses

#  Perform pairwise t-tests using 'genotype.simple' and 'weeks' for
#+ stratification
tbl_rat$genotype_by_weeks <- paste0(tbl_rat$weeks, "_", tbl_rat$genotype.simple)
proper <- c(
    "4–6_Ctrl", "4–6_KO", "7–8_Ctrl", "7–8_KO", "9–12_Ctrl", "9–12_KO",
    "13–26_Ctrl", "13–26_KO", "27–52_Ctrl", "27–52_KO", "53+_Ctrl", "53+_KO"
)
tbl_rat$genotype_by_weeks <- factor(tbl_rat$genotype_by_weeks, levels = proper)

pairwise.t.test(
    tbl_rat$t_to_b,
    tbl_rat$genotype_by_weeks,
    p.adjust.method = "BH"
)


#  Set image dimensions, other information ------------------------------------
image_w <- 6
image_h <- 6
ggplot2::ggsave(
    filename = file.path(dir_script, "plots", "Fig-2D_testis-body_SD.pdf"),
    plot = q,
    units = "in",
    width = image_w,
    height = image_h,
    device = "pdf"
)
