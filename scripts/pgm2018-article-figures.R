
# 0. SETUP -------------------------------------------------------------------

library(tidyverse) # ggplot2 (plotting), dplyr + tidyr (data manipulation)
library(reshape2) # general data reshaping
library(ggthemes) # provides nice themes for ggplots
library(precrec) # needed for computing the ROC / PRC values
library(DescTools) # for computing the Brier score
library(ggpubr) # needed for plot layout (ggarrange)

source('R/utils.R')

# Create folder where to put article figures
figures_dir <- "figures/"
if (!dir.exists(figures_dir)) dir.create(figures_dir)

prior_DAG <- uniform_prior_GRN_DAG()
prior_DMAG <- uniform_prior_GRN_DMAG()

PGM_theme <- theme_classic() + theme(
  text = element_text(size = 15), 
  legend.text = element_text(size = 15), 
  legend.position = "bottom",
  legend.key.size = unit(1.5, "lines"), 
  legend.key.width = unit(1.5, "lines"),
  aspect.ratio = 1, plot.margin = unit(c(0, 0, 0, 0), "lines"))


# 4.1 Consistency of Detecting Local Causal Structures - Figure 3 ---------

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
    ylab("Probability of causal model") + PGM_theme
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

ggsave(paste0(figures_dir, "PGM_Figure_3.pdf"), 
       ggpubr::ggarrange(plotlist = consistency_plots, nrow = 2, ncol = 3),
       width = 16, height = 9)


# 4.2 Causal Discovery in Gene Regulatory Networks - Figure 4 and Table 2 ----

#' Compare different algorithms by plotting precision-recall and ROC curve
#' and by computing Brier scores.
#'
#' @param simulated_GRN_probabilities List containing the structure of a simulated
#' GRN and the posterior probabilities obtained for all simulations.
#' @param simulations Vector specifying different simulations for each algorithm
#' given the data structure (here the sample size is different).
#'
#' @return List containing ROC, PRC, and Brier scores for each simulation.
#'
#' @keywords internal
compare_algorithms_PGM <- function(
  simulated_GRN_probabilities, simulations = c("samples_1e2", "samples_1e3")
  ) {
  
  graph_structure <- simulated_GRN_probabilities$structure
  
  # Plotting parameters
  my_colors <- ggthemes::colorblind_pal()(8)[c(1:3, 6:7, 4)]
  my_labels <- list(
    "ROC" = list(x = "1 - Specificity", y = "Sensitivity"),
    "PRC" = list(x = "Recall", y = "Precision")
  )
  
  # NOTE: sapply as opposed to lapply automatically sets the names of the list using the first argument
  sapply(simulations, function(simulation) {
    
    # extract posterior probabilities list for the current simulation and remove diagonal
    posterior_probabilities_list <- lapply(simulated_GRN_probabilities[[simulation]], off_diagonal)
    
    # Evaluate PRC and ROC measures given scores and labels (graph structure)
    mdat <- precrec::mmdata(
      posterior_probabilities_list, off_diagonal(graph_structure),
      modnames = names(simulated_GRN_probabilities[[simulation]])
    )
    dat <- as.data.frame(evalmod(mdat))
    
    # Plot ROC and PRC curves using precrec::mmdata structure
    result <- sapply(c("ROC", "PRC"), function(type) {
      
      ggplot(dat %>% filter(type ==  !!type), aes(x = x, y = y)) +
        geom_line(aes(color = modname, linetype = modname), size = 0.85) +
        ggtitle(type) + xlab(my_labels[[type]]$x) + ylab(my_labels[[type]]$y) +
        scale_color_manual(name = "", values = my_colors) + ylim(c(0, 1)) +
        scale_linetype_manual(name = "", values = c(3, 3, 3, 1, 1, 2)) + 
        theme_bw() + theme(legend.text = element_text(size = 14)) +
        guides(color = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1))
    }, simplify = FALSE)
    
    # Compute Brier scores
    result$Brier_scores <- sapply(posterior_probabilities_list, function(probs) {
      DescTools::BrierScore(off_diagonal(graph_structure), probs)
    })
    
    result
    
  }, simplify = FALSE) 
  
}

load('data/simulated_GRN_probabilities_PGM.RData')

# Compare algorithms for both the sparse and dense graph
sparse_graph_PGM_results <- compare_algorithms_PGM(simulated_GRN_probabilities_PGM$sparse_graph)
dense_graph_PGM_results <- compare_algorithms_PGM(simulated_GRN_probabilities_PGM$dense_graph)

# Combine PRROC plots with ggarrange, and save them
PGM_PRROC_plot <- ggpubr::ggarrange(
  sparse_graph_PGM_results$samples_1e2$ROC, sparse_graph_PGM_results$samples_1e2$PRC,
  dense_graph_PGM_results$samples_1e2$ROC, dense_graph_PGM_results$samples_1e2$PRC,
  sparse_graph_PGM_results$samples_1e3$ROC, sparse_graph_PGM_results$samples_1e3$PRC,
  dense_graph_PGM_results$samples_1e3$ROC, dense_graph_PGM_results$samples_1e3$PRC,
  nrow = 2, ncol = 4,
  common.legend = TRUE, legend = "bottom"
)
ggsave(paste0(figures_dir, "PGM_Figure_4.pdf"), PGM_PRROC_plot, onefile = FALSE)

# Collect Brier scores and organize them in a tibble
PGM_Brier_scores_sparse_graph <- 
  tibble::as_tibble(t(sapply(sparse_graph_PGM_results, '[[', 'Brier_scores')), rownames = "samples") %>% mutate(graph = "Sparse")
PGM_Brier_scores_dense_graph <- 
  tibble::as_tibble(t(sapply(dense_graph_PGM_results, '[[', 'Brier_scores')), rownames = "samples") %>% mutate(graph = "Dense")
PGM_Brier_scores <- PGM_Brier_scores_sparse_graph %>%
  dplyr::bind_rows(PGM_Brier_scores_dense_graph) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(trigger = mean(c(`trigger 1`, `trigger 2`, `trigger 3`))) %>%
  dplyr::select(samples, graph, trigger, `BFCS DAG`, `BFCS DMAG`, BGe) %>%
  dplyr::mutate(samples = as.integer(substr(samples, 9, 11)))  %>%
  tidyr::pivot_wider(names_from = graph, values_from = trigger:BGe) %>%
  dplyr::relocate(ends_with("Sparse"), .before = trigger_Dense) 

# Output LaTeX table from the Brier scores stored in tibble
# Table format is different from paper but content is identical
knitr::kable(PGM_Brier_scores, format = "latex", digits = 3, booktabs = T, 
             col.names = gsub("_Sparse|_Dense", "", names(PGM_Brier_scores)), align = "c") %>%
  kableExtra::kable_styling() %>%
  kableExtra::add_header_above(c(" " = 1, "Sparse GRN" = 4, "Dense GRN" = 4)) %>% 
  kableExtra::column_spec(c(1, 5), latex_column_spec = "c|") %>%
  kableExtra::save_kable(paste0(figures_dir, "PGM_Table_2.pdf"))


# 4.3. Comparing Results from an Experiment on Yeast - Table 3 ------------

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
  kableExtra::column_spec(1, latex_column_spec = c("l|")) %>%
  kableExtra::save_kable(paste0(figures_dir, "PMG_Table_3a.pdf"))


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
  kableExtra::column_spec(1, latex_column_spec = c("l|")) %>%
  kableExtra::save_kable(paste0(figures_dir, "PGM_Table_3b.pdf"))

