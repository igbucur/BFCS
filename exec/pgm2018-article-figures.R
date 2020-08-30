
# 0. SETUP -------------------------------------------------------------------

library(ggplot2) # plotting the ROC / PRC curves (part of tidyverse)
library(dplyr) # data set processing (part of tidyverse)
library(ggthemes) # provides nice themes for ggplots
library(precrec) # needed for computing the ROC / PRC values
library(DescTools) # for computing the Brier score
library(ggpubr) # needed for plot layout (ggarrange)

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
  ggsave(paste0(figures_dir, "PGM_consistency_", case, ".pdf"), 
         plot_consistency_data(
           get(paste0("Bfs_seq_", case)), 
           priors = list(prior_DMAG), str_name = "DMAG"
         ) + guides(color = FALSE),
         width = 10, height = 5
  )
}

# 2. Figure 4: Comparing the performance of Trigger and BFCS -----------------


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

ggsave(paste0(figures_dir, "PGM_PRC_ROC_curves.pdf"), onefile = FALSE)


# 3. Comparing the calibration using the Brier score ----------------------

#' Plot precision-recall or ROC curve to compare various algorithms.
#'
#' @param df Data frame in long table format containing posterior probabilities 
#' of algorithms under comparison.
#' @param type String specifying type of plot we want to produce (PRC or ROC)
#'
#' @return
#'
#' @keywords internal
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
  tidyr::pivot_wider(names_from = Graph, values_from = trigger:BGe)

# table is in right format, but still needs a bit of work to get to table itself
Hmisc::latex(Brier_scores, file = "", digits = 2)

