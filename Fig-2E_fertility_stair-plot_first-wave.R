
#  Fig-2E_fertility_stair-plot_first-wave.R
#  KA
#  2020-0812, 2024-0210


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


#  Load necessary library
library("tidyverse")


#  Define functions
div <- function(x, y) {
    #  Function for auto-populating values for KO in stair plots
    #+ - x: a vector
    #+ - y: a vector
    seq(min(x), max(x), (max(x) - min(x)) / (length(y) - 1))
}


plot_stairs <- function(x, t, u) {
    #  Function will plot x and y vectors from a tibble or data frame; tuned for
    #+ analyses using two genotype categories: HZ and KO
    #+ - x: a tibble or data frame
    #+ - t: a non-quoted column of x; to be plotted on the x axis
    #+ - u: a non-quoted column of x; to be plotted on the y axis
    #+ 
    #+ stackoverflow.com/questions/57380726/how-to-use-quosure-in-ggplot-to-generate-plotting-function
    enq_t <- dplyr::enquo(t)
    enq_u <- dplyr::enquo(u)
    x %>% ggplot2::ggplot(aes(!!enq_t, !!enq_u)) +
        geom_step(aes(linetype = genotype)) +
        scale_linetype_manual(values = c("dotdash", "dotted")) +
        theme(legend.title = element_blank()) +
        labs(x = enq_t, y = enq_u)
}


#  Load and save an initial tibble
tsv <- "Atf7ip2_fertility_first-wave.tsv"
tbl <- readr::read_tsv(
    file = file.path(dir_script, "tsvs", tsv),
    col_names = TRUE,
    trim_ws = TRUE,
    na = c("", "#N/A"),
    col_types = cols(
        tag = col_factor(),
        genotype = col_factor(),
        date_start = col_date(),
        date_litter = col_date()
    )
)

#  Round the weeks of litter birth
tbl$Weeks <- round(tbl$Weeks)


#  ANALYSIS: Cumulative pups per male

#  Wrangling and cleaning up the data for control males, HZ
#+ 
#+ Create a vector of cumulative sums of pups; process NA as 0
HZ <- tbl %>%
    dplyr::group_by(date_litter) %>%
    subset(genotype == "HZ") %>%
    subset(!is.na(date_litter)) %>%
    dplyr::mutate(litter_pups = tidyr::replace_na(litter_pups, 0))

HZ$litter_pups_acc <- HZ$litter_pups %>% cumsum()

#  Wrangling and cleaning up the data for experiment males, KO
KO <- tbl %>%
    dplyr::group_by(date_litter) %>% 
    subset(genotype == "KO")

#  To visualize a KO line in the step plot, auto-populate cells for week
#+ information
KO$Weeks <- div(HZ$Weeks, KO$Weeks)
KO$litter_pups_acc <- rep(0, length(KO$litter_pups))


#  Plot the stair plot for cumulative pups per male
(
    q <- dplyr::full_join(HZ, KO) %>%
        plot_stairs(t = Weeks, u = litter_pups_acc) +
        labs(x = "Weeks", y = "Cumulative pups per male")
)


#  Set image dimensions, other information, then save the plot
image_w <- 6
image_h <- 4
ggplot2::ggsave(
    filename = file.path(
        dir_script, "plots", "Fig-2E_fertility_stair-plot_first-wave.pdf"
    ),
    plot = q,
    units = "in",
    width = image_w,
    height = image_h,
    device = "pdf"
)
