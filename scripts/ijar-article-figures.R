
# 0. Setup ----------------------------------------------------------------

library(ggplot2) # plotting the ROC / PRC curves (part of tidyverse)
library(dplyr) # data set processing (part of tidyverse)
library(ggthemes) # provides nice themes for ggplots
library(precrec) # needed for computing the ROC / PRC values
library(DescTools) # for computing the Brier score
library(ggpubr) # needed for plot layout (ggarrange)
library(reshape2) # needed for reshaping data frames

source('R/compute_prior_structures.R')
source('R/utils.R')

figures_dir <- "figures/"
if (!dir.exists(figures_dir)) dir.create(figures_dir)

prior_DAG <- uniform_prior_GRN_DAG()
prior_DMAG <- uniform_prior_GRN_DMAG()


# 4.1. Consistency of Detecting Local Causal Structures - Figure 3 --------

#' Produce plot that shows the consistency of BFCS using simulation experiments.
#'
#' @param Bayes_factors List of computed Bayes factors for various number of samples
#' @param prior Prior for causal structures over triplets.
#'
#' @return ggplot showing consistency of BFCS in finding the correct causal model
#' 
#' @keywords internal
plot_consistency_data <- function(Bayes_factors, prior = prior_DMAG) {
  
  posterior_probabilities <- sapply(Bayes_factors, function(Bfs) {
    Bfs <- Bfs * prior
    Bfs[6, ] / colSums(Bfs)
  }) %>%
    reshape2::melt() %>%
    dplyr::rename(num_rep = Var1, num_obs = Var2)
  
  ggplot(posterior_probabilities, aes(x = log10(num_obs), y = value, group = log10(num_obs))) +
    geom_boxplot(outlier.alpha = 0.1, color = 'darkblue', fill = 'lightblue') + ylim(c(0, 1)) +
    xlab("Number of samples on the base-10 logarithmic scale") +
    ylab("Probability of causal model") + IJAR_theme
}

load('data/Bayes_factors_consistency.RData')

consistency_plots <- list(
  plot_consistency_data(Bayes_factors_consistency$gaussian$causal),
  plot_consistency_data(Bayes_factors_consistency$gaussian$independent),
  plot_consistency_data(Bayes_factors_consistency$gaussian$full),
  plot_consistency_data(Bayes_factors_consistency$binomial$causal),
  plot_consistency_data(Bayes_factors_consistency$binomial$independent),
  plot_consistency_data(Bayes_factors_consistency$binomial$full)
)

ggsave(paste0(figures_dir, "IJAR_Figure_3.pdf"), 
       ggpubr::ggarrange(plotlist = consistency_plots, nrow = 2, ncol = 3),
       width = 16, height = 9)


# 4.2. Causal Discovery in Gene Regulatory Networks - Figures 4-7 ---------

run_comparison <- function(simulated_GRN_probabilities, graph_structure_type, simulations = c("samples_1e2", "samples_1e3")) {
  
  plot_PRROC_curve <- function(df, type = "ROC") {
    
    my_labels <- if (type == "ROC") {
      list(x = "1 - Specificity", y = "Sensitivity")
    } else if (type == "PRC") {
      list(x = "Recall", y = "Precision")
    } else {
      stop("unknown curve type")
    }
    
    # TODO: make scale_size work somehow
    ggplot(df %>% filter(type == !!type), aes(x = x, y = y)) +
      geom_line(aes(color = modname, linetype = modname, size = modname)) +
      ylim(c(0, 1)) + xlab(my_labels$x) + ylab(my_labels$y) +
      scale_color_manual(name = "", labels = method_names, values = my_colors) +
      scale_linetype_manual(name = "", labels = method_names, values = c(2, 1, 1, 1, 2)) +
      scale_size_manual(name = "", labels = method_names, values = c(0.5, 1, 1, 1, 0.5)) +
      guides(color = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1)) +
      ggtitle(type) + 
      IJAR_theme
  }
  
  calib_line_plot <- function(calib_data) {

    calib_data %>% 
      filter(Count != 0) %>%
      ggplot(aes(x = midpoint / 100, y = Percent / 100, color = calibModelVar, linetype = calibModelVar, size = calibModelVar)) + 
      geom_line(position = position_dodge(width = 0.1)) +
      geom_errorbar(aes(ymin = Lower / 100, ymax = Upper / 100), width= .05, position = position_dodge(width = 0.1)) +
      geom_abline(slope = 1, intercept = 0, linetype = 3, col = 'gray') +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
      scale_color_manual(name = "", labels = method_names, values = my_colors) +
      scale_linetype_manual(name = "", labels = method_names, values = c(2, 1, 1, 1, 2)) +
      scale_size_manual(name = "", labels = method_names, values = c(0.5, 1, 1, 1, 0.5)) +
      xlab("Bin Midpoint Average Estimated Percentage") + 
      ylab("Observed Event Percentage") +
      guides(color = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1)) +
      ggtitle("Calibration") + 
      IJAR_theme
  }
  

  
  if (graph_structure_type == 'direct') {
    graph_structure <- simulated_GRN_probabilities$structure$direct
  } else if (graph_structure_type == 'ancestral') {
    graph_structure <- simulated_GRN_probabilities$structure$ancestral
  } else {
    stop("Wrong graph structure type specified. Possible options are 'direct' and 'ancestral'.")
  }
  

  
  simulated_GRN_plots <- vector("list", length = 3 * length(simulations))
  
  for (i in 1:length(simulations)) {
    
    simulation_probabilities <- simulated_GRN_probabilities[[simulations[i]]]
    
    method_names <- names(simulation_probabilities)
    sanitized_method_names <- gsub(" ", "_", method_names)
    my_colors <- ggthemes::colorblind_pal()(length(method_names))
    
    probabilities_off_diagonal <- lapply(simulation_probabilities, off_diagonal)
    structure_off_diagonal <- off_diagonal(graph_structure)
    
    # Create PRROC curves
    mdat <- precrec::mmdata(probabilities_off_diagonal, structure_off_diagonal, 
                            modnames = sanitized_method_names)
    curves <- precrec::evalmod(mdat)
    PRROC_df <- as.data.frame(curves)
    
    simulated_GRN_plots[[i]] <- plot_PRROC_curve(PRROC_df, type = "ROC")
    simulated_GRN_plots[[i + length(simulations)]] <- plot_PRROC_curve(PRROC_df, type = "PRC")

    
    # Create input data frame for caret::calibration method
    calib_df <- data.frame(do.call('cbind', probabilities_off_diagonal))
    names(calib_df) <- sanitized_method_names
    calib_df$edge <- as.factor(structure_off_diagonal)
    # Compute calibration using caret::calibration
    calib_result <- caret::calibration(
      as.formula(paste("edge ~", paste(sanitized_method_names, collapse = " + "))), 
      cuts = 5, class = 1, data = calib_df)
    
    simulated_GRN_plots[[i + 2 * length(simulations)]] <- calib_line_plot(calib_result$data)
  }

  
  ggpubr::ggarrange(plotlist = simulated_GRN_plots,
    nrow = 3, ncol = 2,
    common.legend = TRUE, legend = "bottom"
  )
}

load('data/simulated_GRN_probabilities_IJAR.RData')


# Figure 4 : Sparser graph, direct relationships
ggsave(paste0(figures_dir, "IJAR_Figure_4.pdf"), 
       run_comparison(simulated_GRN_probabilities_IJAR$sparse_graph, "direct"),
       onefile = FALSE, width = 10, height = 15)
# Figure 4 : Denser graph, direct relationships
ggsave(paste0(figures_dir, "IJAR_Figure_5.pdf"), 
       run_comparison(simulated_GRN_probabilities_IJAR$less_dense_graph, "direct"),
       onefile = FALSE, width = 10, height = 15)
# Figure 4 : Sparser graph, ancestral relationships
ggsave(paste0(figures_dir, "IJAR_Figure_6.pdf"), 
       run_comparison(simulated_GRN_probabilities_IJAR$sparse_graph, "ancestral"),
       onefile = FALSE, width = 10, height = 15)
# Figure 4 : Denser graph, ancestral relationships
ggsave(paste0(figures_dir, "IJAR_Figure_7.pdf"), 
       run_comparison(simulated_GRN_probabilities_IJAR$less_dense_graph, "ancestral"),
       onefile = FALSE, width = 10, height = 15)



# 4.3. Comparing Results from an Experiment on Yeast - Table 2 ------------

# Read posterior probabilities for causal regulatory relationships
load('data/yeast_regulatory_probabilities.RData')

# Top ranked regulated genes are taken from Chen et al. (2008)
NAM9_regulated_Chen <- c("MDM35", "CBP6", "QRI5", "RSM18", "RSM7", "MRPL11", "MRPL25",
                         "DLD2", "YPR126C", "MSS116", "AFG3", "SSC1", "MRPL33", "YPR125W")

# Create data frame containing probabilities
compare_NAM9_regulated_sort_Chen <- data.frame(
  Rank = 1:length(NAM9_regulated_Chen),
  Gene = NAM9_regulated_Chen,
  `Chen et al.` = yeast_trigger_Chen["NAM9", NAM9_regulated_Chen],
  trigger = yeast_trigger_w50k["NAM9", NAM9_regulated_Chen],
  BFCS = yeast_BFCS_DMAG["NAM9", NAM9_regulated_Chen],
  check.names = FALSE
)

# Create, format and save LaTeX table 
knitr::kable(compare_NAM9_regulated_sort_Chen[1:10, ], format = "latex", digits = 3, 
             booktabs = TRUE, row.names = FALSE, align = "l") %>%
  kableExtra::kable_styling() %>%
  kableExtra::save_kable(paste0(figures_dir, "IJAR_Table_2a.pdf"))


# Top ranked regulated genes according to the BFCS DMAG probabilities
NAM9_regulated_BFCS <- order(yeast_BFCS_DMAG["NAM9", ], decreasing = TRUE)[1:length(NAM9_regulated_Chen)]

# Create data frame containing probabilities
compare_NAM9_regulated_sort_BFCS <- data.frame(
  Rank = 1:length(NAM9_regulated_BFCS),
  Gene = NAM9_regulated_BFCS,
  `Chen et al.` = yeast_trigger_Chen["NAM9", NAM9_regulated_BFCS],
  trigger = yeast_trigger_w50k["NAM9", NAM9_regulated_BFCS],
  BFCS = yeast_BFCS_DMAG["NAM9", NAM9_regulated_BFCS],
  check.names = FALSE
)

# Create, format and save LaTeX table 
knitr::kable(compare_NAM9_regulated_sort_BFCS[1:10, ], format = "latex", digits = 3, 
             booktabs = TRUE, row.names = FALSE, align = "l") %>%
  kableExtra::kable_styling() %>%
  kableExtra::save_kable(paste0(figures_dir, "IJAR_Table_2b.pdf"))

# 5. Computational and Time Complexity - Figure 8 and 9 -------------------

load('data/computational_and_time_complexity.RData')

time_network_size_plot <- time_versus_network_size %>%
  dplyr::filter(row_number() <= 30) %>%
  tidyr::pivot_longer(cols = BFCS.elapsed:trigger.elapsed) %>%
  ggplot2::ggplot() +
  geom_line(aes(x = num_var, y = value, col = name), size = 1) + IJAR_theme + 
  guides(col = guide_legend(title = "Method", labels = c("BFCS", "BGe", "trigger"))) +
  ylab("Time in seconds") +
  xlab("Network size (nodes)") +
  scale_color_manual(labels = c("BFCS", "BGe", "trigger"), values = c(1, 2, 3)) +
  scale_y_log10() + IJAR_theme +
  theme(text = element_text(size = 20), legend.text = element_text(size = 20))

ggsave(paste0(figures_dir, 'IJAR_Figure_8.pdf'), time_network_size_plot)


time_sample_size_plot <- time_versus_sample_size %>%
  tidyr::pivot_longer(cols = BFCS.elapsed:trigger.elapsed) %>%
  ggplot2::ggplot() +
  geom_line(aes(x = num_obs, y = value, col = name), size = 1) + IJAR_theme + 
  guides(col = guide_legend(title = "Method", labels = c("BFCS", "BGe", "trigger"))) +
  ylab("Time in seconds") +
  xlab("Number of observations") +
  scale_color_manual(labels = c("BFCS", "BGe", "trigger"), values = c(1, 2, 3)) +
  scale_x_log10() + scale_y_log10() + IJAR_theme +
  theme(text = element_text(size = 20), legend.text = element_text(size = 20))

ggsave(paste0(figures_dir, 'IJAR_Figure_9.pdf'), time_sample_size_plot)



