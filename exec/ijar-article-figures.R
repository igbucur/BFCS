
# 0. Setup ----------------------------------------------------------------

library(ggplot2) # plotting the ROC / PRC curves (part of tidyverse)
library(dplyr) # data set processing (part of tidyverse)
library(ggthemes) # provides nice themes for ggplots
library(precrec) # needed for computing the ROC / PRC values
library(DescTools) # for computing the Brier score
library(ggpubr) # needed for plot layout (ggarrange)
library(reshape2) # needed for reshaping data frames

# off_diagonal <- function(m) m[lower.tri(m) | upper.tri(m)]
# TPR <- function(prob, g) sum(prob * g) / sum(g)
# FPR <- function(prob, g) sum(prob * (1 - g)) / sum(1 - g)
# sG_structure <- as.matrix(read.table('inst/extdata/sG_structure.txt'))
# dG_structure <- as.matrix(read.table('inst/extdata/dG_structure.txt'))

source('R/compute_prior_structures.R')

figures_dir <- "man/figures/"
if (!dir.exists(figures_dir)) dir.create(figures_dir)

prior_DAG <- uniform_prior_GRN_DAG()
prior_DMAG <- uniform_prior_GRN_DMAG()
prior_ADMG <- uniform_prior_GRN_ADMG()

# For use when making ggplots
IJAR_theme <- theme_classic() + theme(
  text = element_text(size = 15), 
  legend.text = element_text(size = 15), 
  legend.position = "bottom",
  legend.key.size = unit(1.5, "lines"), 
  legend.key.width = unit(1.5, "lines"),
  aspect.ratio = 1, plot.margin = unit(c(0, 0, 0, 0), "lines"))



# 1. Figure 3: Consistency of BFCS -------------------------------------------

#' Produce plot that shows the consistency of BFCS using simulation experiments.
#'
#' @param bayes_factors List of computed Bayes factors for various number of samples
#' @param priors List of priors to combine with Bayes factors for deriving posterior
#' over causal structures
#' @param str_name 
#'
#' @return ggplot showing consistency of BFCS in finding the correct causal model
#' 
#' @keywords internal
plot_consistency_data <- function(bayes_factors, priors = list(prior_DAG, prior_DMAG, prior_ADMG), 
                                  str_name = c('DAG', 'DMAG', 'ADMG')) {
  
  df <- do.call('rbind', lapply(1:length(priors), function(i) {
    prob <- sapply(bayes_factors, function(Bfs) {
      Bfs <- Bfs * priors[[i]]
      Bfs[6, ] / colSums(Bfs)
    })
    prob %>%
      melt() %>%
      rename(rep = Var1, nobs = Var2) %>%
      mutate(nobs = nobs_seq[nobs]) %>%
      mutate(str = str_name[i])
  }))
  
  ggplot(df, aes(x = log10(nobs), y = value, group = interaction(nobs, str), color = str)) +
    geom_boxplot(outlier.alpha = 0.1, color = 'darkblue', fill = 'lightblue') + 
    ylim(c(0, 1)) +
    xlab("Number of samples on the base-10 logarithmic scale") +
    ylab("Probability of causal model") +
    labs(color = "Causal structure") + 
    theme_classic() +
    theme(legend.position = "top", text = element_text(size = 20),
          axis.title = element_text(size = 25),
          axis.text = element_text(size = 25))
}

nobs_seq <- c(100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000)

load('data/consistency_simulation_experiment_norm.RData')
load('data/consistency_simulation_experiment_binom.RData')

cases <- paste(rep(c("T", "I", "F"), 2), rep(c("norm", "binom"), each = 3), sep = "_")

for(case in cases) {
  ggsave(paste0(figures_dir, "IJAR_consistency_", case, ".pdf"), 
         plot_consistency_data(
           get(paste0("Bfs_seq_", case)), 
           priors = list(prior_DMAG), str_name = "DMAG"
         ) + guides(color = FALSE),
         width = 10, height = 5
  )
}

# 2. Figures 4-7: Comparing the performance of Trigger and BFCS ------------

off_diagonal <- function(m) m[lower.tri(m) | upper.tri(m)]

suffix <- '_all'
sG_structure <- as.matrix(read.table(paste0('inst/extdata/sG_structure_runif', suffix, '.txt')))
lG_structure <- as.matrix(read.table(paste0('inst/extdata/lG_structure_runif', suffix, '.txt')))
dG_structure <- as.matrix(read.table(paste0('inst/extdata/dG_structure_runif', suffix, '.txt')))
sG_ancestral <- as.matrix(Matrix::expm(sG_structure) > 0) * 1
lG_ancestral <- as.matrix(Matrix::expm(lG_structure) > 0) * 1
dG_ancestral <- as.matrix(Matrix::expm(dG_structure) > 0) * 1

methods <- data.frame(
  name = c('trigger', rep('BFCS', 2), 'BFCS_loc', 'BGe'),
  type = c('avg', 'DAG', 'DMAG', 'DMAG', 'DMAG')
)

run_comparison <- function(data_filenames, method_names, graph_structure, 
                           legend_names = method_names) {
  
  my_colors <- ggthemes::colorblind_pal()(length(method_names))
  # my_colors <- c('black', 'light blue', 'dark blue', 'green', 'orange')
  
  plot_PRROC_curve <- function(df, type = "ROC") {
    
    my_labels <- if (type == "ROC") {
      list(x = "1 - Specificity", y = "Sensitivity")
    } else if (type == "PRC") {
      list(x = "Recall", y = "Precision")
    } else {
      stop("unknown curve type")
    }
    
    require(dplyr) 
    require(ggplot2)
    
    # TODO: make scale_size work somehow
    ggplot(df %>% filter(type == !!type), aes(x = x, y = y)) +
      geom_line(aes(color = modname, linetype = modname, size = modname)) +
      ylim(c(0, 1)) + xlab(my_labels$x) + ylab(my_labels$y) +
      scale_color_manual(name = "", labels = legend_names, values = my_colors) +
      scale_linetype_manual(name = "", labels = legend_names, values = c(2, 1, 1, 1, 2)) +
      scale_size_manual(name = "", labels = legend_names, values = c(0.5, 1, 1, 1, 0.5)) +
      guides(color = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1)) +
      ggtitle(type) + 
      IJAR_theme
  }
  
  calib_line_plot <- function(calib_data) {
    require(dplyr)
    require(ggplot2)
    calib_data %>% 
      filter(Count != 0) %>%
      ggplot(aes(x = midpoint / 100, y = Percent / 100, color = calibModelVar, linetype = calibModelVar, size = calibModelVar)) + 
      geom_line(position = position_dodge(width = 0.1)) +
      geom_errorbar(aes(ymin = Lower / 100, ymax = Upper / 100), width= .05, position = position_dodge(width = 0.1)) +
      geom_abline(slope = 1, intercept = 0, linetype = 3, col = 'gray') +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      scale_color_manual(name = "", labels = legend_names, values = my_colors) +
      scale_linetype_manual(name = "", labels = legend_names, values = c(2, 1, 1, 1, 2)) +
      scale_size_manual(name = "", labels = legend_names, values = c(0.5, 1, 1, 1, 0.5)) +
      xlab("Bin Midpoint Average Estimated Percentage") + 
      ylab("Observed Event Percentage") +
      guides(color = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1)) +
      ggtitle("Calibration") + 
      IJAR_theme
  }
  
  stopifnot(length(data_filenames) == length(method_names))
  
  probs_list <- lapply(data_filenames, function(file) {
    off_diagonal(as.matrix(read.table(file)))
  })
  
  mdat <- precrec::mmdata(probs_list, off_diagonal(graph_structure), modnames = method_names)
  curves <- precrec::evalmod(mdat)
  df <- as.data.frame(curves)
  
  plot_ROC <- plot_PRROC_curve(df, type = "ROC")
  plot_PRC <- plot_PRROC_curve(df, type = "PRC")
  Brier_scores <- sapply(probs_list, DescTools::BrierScore, resp = off_diagonal(graph_structure))
  names(Brier_scores) <- method_names
  
  calib_df <- data.frame(do.call('cbind', probs_list))
  names(calib_df) <- method_names
  calib_df$edge <- as.factor(off_diagonal(graph_structure))
  
  # print(names(my_df))
  
  calib_result <- caret::calibration(as.formula(paste("edge ~", paste(method_names, collapse = " + "))), cuts = 5, class = 1, data = calib_df)
  
  calib_plot <- calib_line_plot(calib_result$data)
  
  
  list(plot_ROC = plot_ROC, plot_PRC = plot_PRC, calibration = Brier_scores,
       AUC = auc(curves), plot_clb = calib_plot)
}


# Figure 4 : Sparser graph, direct relationships

results_folder <- paste0('results/runif', suffix, '/')

sG_n1e2_dir <- run_comparison(
  apply(methods, 1, function(row) paste0(results_folder, "/", row[1], "_sparse_n1e2_", row[2], ".txt")),
  apply(methods, 1, function(row) paste0(row[1], "_", row[2])),
  sG_structure, legend_names = c("trigger", "BFCS DAG", "BFCS DMAG", "BFCS loclink", "BGe")
)

sG_n1e3_dir <- run_comparison(
  apply(methods, 1, function(row) paste0(results_folder, "/", row[1], "_sparse_n1e3_", row[2], ".txt")),
  apply(methods, 1, function(row) paste0(row[1], "_", row[2])),
  sG_structure, legend_names = c("trigger", "BFCS DAG", "BFCS DMAG", "BFCS loclink", "BGe")
)

plot_sG_dir <- ggpubr::ggarrange(
  sG_n1e2_dir$plot_ROC, sG_n1e3_dir$plot_ROC,
  sG_n1e2_dir$plot_PRC, sG_n1e3_dir$plot_PRC,
  sG_n1e2_dir$plot_clb, sG_n1e3_dir$plot_clb,
  nrow = 3, ncol = 2,
  common.legend = TRUE, legend = "bottom"
)
ggsave(paste0(figures_dir, "IJAR_sG_dir_runif", suffix, ".pdf"), plot_sG_dir, onefile = FALSE, width = 10, height = 15)


# Figure 5 : Sparser graph, ancestral relationships

lG_n1e2_dir <- run_comparison(
  apply(methods, 1, function(row) paste0(results_folder, "/", row[1], "_less_dense_n1e2_", row[2], ".txt")),
  apply(methods, 1, function(row) paste0(row[1], "_", row[2])),
  lG_structure, legend_names = c("trigger", "BFCS DAG", "BFCS DMAG", "BFCS loclink", "BGe")
)

lG_n1e3_dir <- run_comparison(
  apply(methods, 1, function(row) paste0(results_folder, "/", row[1], "_less_dense_n1e3_", row[2], ".txt")),
  apply(methods, 1, function(row) paste0(row[1], "_", row[2])),
  lG_structure, legend_names = c("trigger", "BFCS DAG", "BFCS DMAG", "BFCS loclink", "BGe")
)

plot_lG_dir <- ggpubr::ggarrange(
  lG_n1e2_dir$plot_ROC, lG_n1e3_dir$plot_ROC, 
  lG_n1e2_dir$plot_PRC, lG_n1e3_dir$plot_PRC, 
  lG_n1e2_dir$plot_clb, lG_n1e3_dir$plot_clb,
  nrow = 3, ncol = 2,
  common.legend = TRUE, legend = "bottom"
)
ggsave(paste0(figures_dir, "IJAR_lG_dir_runif", suffix, ".pdf"), plot_lG_dir, onefile = FALSE, width = 10, height = 15)


# Figure 6 : Denser graph, direct relationships

sG_n1e2_anc <- run_comparison(
  apply(methods, 1, function(row) paste0(results_folder, "/", row[1], "_sparse_n1e2_", row[2], ".txt")),
  apply(methods, 1, function(row) paste0(row[1], "_", row[2])),
  sG_ancestral, legend_names = c("trigger", "BFCS DAG", "BFCS DMAG", "BFCS loclink", "BGe")
)

sG_n1e3_anc <- run_comparison(
  apply(methods, 1, function(row) paste0(results_folder, "/", row[1], "_sparse_n1e3_", row[2], ".txt")),
  apply(methods, 1, function(row) paste0(row[1], "_", row[2])),
  sG_ancestral, legend_names = c("trigger", "BFCS DAG", "BFCS DMAG", "BFCS loclink", "BGe")
)

plot_sG_anc <- ggpubr::ggarrange(
    sG_n1e2_anc$plot_ROC, sG_n1e3_anc$plot_ROC, 
    sG_n1e2_anc$plot_PRC, sG_n1e3_anc$plot_PRC, 
    sG_n1e2_anc$plot_clb, sG_n1e3_anc$plot_clb,
    nrow = 3, ncol = 2,
    common.legend = TRUE, legend = "bottom"
)
ggsave(paste0(figures_dir, "IJAR_sG_anc_runif", suffix, ".pdf"), plot_sG_anc, onefile = FALSE, width = 10, height = 15)


# Figure 7 : Denser graph, ancestral relationships

lG_n1e2_anc <- run_comparison(
  apply(methods, 1, function(row) paste0(results_folder, "/", row[1], "_less_dense_n1e2_", row[2], ".txt")),
  apply(methods, 1, function(row) paste0(row[1], "_", row[2])),
  lG_ancestral, legend_names = c("trigger", "BFCS DAG", "BFCS DMAG", "BFCS loclink", "BGe")
)

lG_n1e3_anc <- run_comparison(
  apply(methods, 1, function(row) paste0(results_folder, "/", row[1], "_less_dense_n1e3_", row[2], ".txt")),
  apply(methods, 1, function(row) paste0(row[1], "_", row[2])),
  lG_ancestral, legend_names = c("trigger", "BFCS DAG", "BFCS DMAG", "BFCS loclink", "BGe")
)


plot_lG_anc <- ggpubr::ggarrange(
  lG_n1e2_anc$plot_ROC, lG_n1e3_anc$plot_ROC, 
  lG_n1e2_anc$plot_PRC, lG_n1e3_anc$plot_PRC, 
  lG_n1e2_anc$plot_clb, lG_n1e3_anc$plot_clb,
  nrow = 3, ncol = 2,
  common.legend = TRUE, legend = "bottom"
)
ggsave(paste0(figures_dir, "IJAR_lG_anc_runif", suffix, ".pdf"), plot_lG_anc, onefile = FALSE, width = 10, height = 15)



# 3. Figure 8 - Time complexity as a function of network size -------------

load('data/times_all_methods.RData')

times_less_dense <- tibble::as_tibble(do.call('rbind', times_all_methods[2, ]))

time_network_size_plot <- times_less_dense %>%
  dplyr::filter(row_number() <= 30) %>%
  dplyr::mutate(size = row_number() * 5) %>%
  tidyr::pivot_longer(cols = BFCS.elapsed:trigger.elapsed) %>%
  ggplot2::ggplot() +
  geom_line(aes(x = size, y = value, col = name), size = 1) + IJAR_theme + 
  guides(col = guide_legend(title = "Method", labels = c("BFCS", "BGe", "trigger"))) +
  ylab("Time in seconds") +
  xlab("Network size (nodes)") +
  scale_color_manual(labels = c("BFCS", "BGe", "trigger"), values = c(1, 2, 3)) +
  scale_y_log10() + IJAR_theme +
  theme(text = element_text(size = 20), legend.text = element_text(size = 20))

ggsave(paste0(figures_dir, 'IJAR_time_versus_network_size.pdf'), time_network_size_plot)


# 4. Figure 9 - Time complexity as a function of sample size --------------

time_sample_size_plot <- tibble::as_tibble(t(readRDS('data/performance_vs_nobs.rds'))) %>%
  tidyr::pivot_longer(cols = BFCS.elapsed:trigger.elapsed) %>%
  ggplot2::ggplot() +
  geom_line(aes(x = nobs, y = value, col = name), size = 1) + IJAR_theme + 
  guides(col = guide_legend(title = "Method", labels = c("BFCS", "BGe", "trigger"))) +
  ylab("Time in seconds") +
  xlab("Number of observations") +
  scale_color_manual(labels = c("BFCS", "BGe", "trigger"), values = c(1, 2, 3)) +
  scale_x_log10() + scale_y_log10() + IJAR_theme +
  theme(text = element_text(size = 20), legend.text = element_text(size = 20))

ggsave(paste0(figures_dir, 'IJAR_time_versus_sample_size.pdf'), time_sample_size_plot)

