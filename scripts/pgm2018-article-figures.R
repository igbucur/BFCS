
# 0. SETUP -------------------------------------------------------------------

library(tidyverse) # ggplot2 (plotting), dplyr + tidyr (data manipulation)
library(reshape2) # general data reshaping
library(ggthemes) # provides nice themes for ggplots
library(precrec) # needed for computing the ROC / PRC values
library(DescTools) # for computing the Brier score
library(ggpubr) # needed for plot layout (ggarrange)

off_diagonal <- function(m) m[lower.tri(m) | upper.tri(m)]
# TPR <- function(prob, g) sum(prob * g) / sum(g)
# FPR <- function(prob, g) sum(prob * (1 - g)) / sum(1 - g)
# sG_structure <- as.matrix(read.table('inst/extdata/sG_structure.txt'))
# dG_structure <- as.matrix(read.table('inst/extdata/dG_structure.txt'))

source('R/compute_prior_structures.R')

# Create folder where to put article figures
figures_dir <- "figures/"
if (!dir.exists(figures_dir)) dir.create(figures_dir)

prior_DAG <- uniform_prior_GRN_DAG()
prior_DMAG <- uniform_prior_GRN_DMAG()
# prior_ADMG <- uniform_prior_GRN_ADMG()



# 4.1 Consistency of Detecting Local Causal Structures - Figure 3 ---------

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
plot_consistency_data <- function(Bayes_factors, prior = prior_DMAG) {
  
  posterior_probabilities <- sapply(Bayes_factors, function(Bfs) {
    Bfs <- Bfs * prior
    Bfs[6, ] / colSums(Bfs)
  }) %>%
    reshape2::melt() %>%
    dplyr::rename(rep = Var1, nobs = Var2) %>%
    dplyr::mutate(nobs = nobs_seq[nobs])
  
  ggplot(posterior_probabilities, aes(x = log10(nobs), y = value, group = log10(nobs))) +
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

ggsave(paste0(figures_dir, "PGM_Figure_3.pdf"), 
       ggpubr::ggarrange(plotlist = consistency_plots, nrow = 2, ncol = 3),
       width = 16, height = 9)


# 4.2 Causal Discovery in Gene Regulatory Networks - Figure 4 -------------

#' Plot precision-recall or ROC curve to compare various algorithms.
#'
#' @param df Data frame in long table format containing posterior probabilities 
#' of algorithms under comparison.
#' @param type String specifying type of plot we want to produce (PRC or ROC)
#'
#' @return
#'
#' @keywords internal
plot_PRROC_curves <- function(
  datafile_root = 'inst/extdata/', simulations = c("sG_1e2", "sG_1e3", "dG_1e2", "dG_1e3"),
  algorithms = c(paste0("trigger_", 1:3), "BFCS_DAG", "BFCS_DMAG", "BGe_scaled")
  ) {
  
  library(ggplot2)
  library(precrec)
  library(dplyr)
  
  # if (!(simulation %in% c('sG_1e2', 'dG_1e2', 'sG_1e3', 'dG_1e3'))) {
  #   stop("Invalid simulation")
  # }
  
  off_diagonal <- function(m) m[lower.tri(m) | upper.tri(m)]
  my_colors <- ggthemes::colorblind_pal()(8)[c(1:3, 6:7, 4)]
  
  my_labels <- list(
    "ROC" = list(x = "1 - Specificity", y = "Sensitivity"),
    "PRC" = list(x = "Recall", y = "Precision")
  )
  
  sapply(simulations, function(simulation) {
    posterior_probabilities_list <- lapply(
      algorithms,
      function(algorithm) {
        filename <- paste0(datafile_root, 'PGM_prob_', simulation, '_', algorithm, '.txt')
        off_diagonal(as.matrix(read.table(filename)))
      }
    )
    
    graph_structure <- as.matrix(read.table(paste0(datafile_root, substr(simulation, 1, 2), '_structure.txt')))
    
    mdat <- precrec::mmdata(
      posterior_probabilities_list, off_diagonal(graph_structure),
      modnames = algorithms
    )
    dat <- as.data.frame(evalmod(mdat))
    
    
    sapply(c("ROC", "PRC"), function(type) {
      
      ggplot(dat %>% filter(type ==  !!type), aes(x = x, y = y)) +
        geom_line(aes(color = modname, linetype = modname), size = 0.85) +
        ggtitle(type) + xlab(my_labels[[type]]$x) + ylab(my_labels[[type]]$y) +
        scale_color_manual(name = "", values = my_colors) + ylim(c(0, 1)) +
        scale_linetype_manual(name = "", values = c(3, 3, 3, 1, 1, 2)) + 
        theme_bw() + theme(legend.text = element_text(size = 14)) +
        guides(color = guide_legend(nrow = 1, keywidth = 2.5, keyheight = 1))
    }, simplify = FALSE)
  }, simplify = FALSE) 
}

curves <- plot_PRROC_curves()
combined_plot <- ggpubr::ggarrange(
  curves$sG_1e2$ROC, curves$sG_1e2$PRC,
  curves$dG_1e2$ROC, curves$dG_1e2$PRC,
  curves$sG_1e3$ROC, curves$sG_1e3$PRC,
  curves$dG_1e3$ROC, curves$dG_1e3$PRC,
  nrow = 2, ncol = 4,
  common.legend = TRUE, legend = "bottom"
)

ggsave(paste0(figures_dir, "PGM_Figure_4.pdf"), onefile = FALSE)



# 4.2 Causal Discovery in Gene Regulatory Networks - Table 2 --------------

#' Compute the Brier calibration score to compare the algorithms.
#'
#' @param datafile_root string; root directory for data files
#' @param simulations character vector; simulation acronyms
#' @param algorithms character vector; algorithm acronyms
#' 
#' @return
compute_Brier_scores <- function(
  datafile_root = 'inst/extdata/', simulations = c("sG_1e2", "sG_1e3", "dG_1e2", "dG_1e3"),
  algorithms = c(paste0("trigger_", 1:3), "BFCS_DAG", "BFCS_DMAG", "BGe_scaled")
) {
  
  t(sapply(simulations, function(simulation) {
    
    graph_structure <- as.matrix(read.table(paste0(datafile_root, substr(simulation, 1, 2), '_structure.txt')))
    

    Brier_scores <- sapply(
      algorithms,
      function(algorithm) {
        filename <- paste0(datafile_root, 'PGM_prob_', simulation, '_', algorithm, '.txt')
        DescTools::BrierScore(off_diagonal(graph_structure), off_diagonal(as.matrix(read.table(filename))))
      }
    )
  
  }))
}

Brier_scores <- tibble::as_tibble(compute_Brier_scores(), rownames = "simulation") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(trigger = mean(c(trigger_1, trigger_2, trigger_3))) %>%
  dplyr::select(simulation, trigger, BFCS_DAG, BFCS_DMAG, BGe_scaled) %>%
  dplyr::rename(BGe = BGe_scaled) %>%
  dplyr::mutate(Samples = as.integer(substr(simulation, 4, 6))) %>%
  dplyr::mutate(Graph = ifelse(substr(simulation, 1, 1) == "s", "Sparse", "Dense")) %>%
  dplyr::select(-simulation) %>%
  tidyr::pivot_wider(names_from = Graph, values_from = trigger:BGe) %>%
  dplyr::relocate(ends_with("Sparse"), .before = trigger_Dense) 

table_cols <- gsub("_", " ", gsub("_Sparse|_Dense", "", names(Brier_scores)))

# Table format is different from paper but content is identical
knitr::kable(Brier_scores, format = "latex", digits = 3, booktabs = T, col.names = table_cols, align = "c") %>%
  kableExtra::kable_styling() %>%
  kableExtra::add_header_above(c(" " = 1, "Sparse GRN" = 4, "Dense GRN" = 4)) %>% 
  kableExtra::column_spec(c(1, 5), latex_column_spec = c("c||", "c|")) %>%
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

