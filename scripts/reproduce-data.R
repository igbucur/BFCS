

#' @author Ioan Gabriel Bucur
#' @description This script performs the simulations described in Subsection 4.1
#' of the paper titled "A Bayesian Approach for Inferring Local Causal Structure
#' in Gene Regulatory Networks" and its extension "Large-Scale Local Causal 
#' Inference of Gene Regulatory Relationships". The purpose of this simulation 
#' is to empirically show that BFCS consistently recovers the correct local 
#' causal structure for multivariate Gaussian and conditionally Gaussian data.


# 0. Setup ----------------------------------------------------------------

## Load libraries and scripts
source('R/compute_Bayes_factors.R')
source('R/compute_prior_structures.R')
source('R/compute_BGe_score_vectorized.R')
source('R/simulate_GRN.R')
source('R/infer_GRN.R')
source('R/utils.R')
library(pbmcapply)
library(trigger)


Rcpp::sourceCpp('src/compute_BGe_score.cpp')

## Global variables
prior_DAG <- uniform_prior_GRN_DAG()
prior_DMAG <- uniform_prior_GRN_DMAG()


# 4.1. Consistency of Detecting Local Causal Structures - Figure 3 PGM / 3 IJAR ----

# Number of replications in each simulation
num_reps <- 1000
# Number of observations sequence increases logarithmically for each set of simulations
num_obs_seq <- c(100, 300, 1000, 3000)
# num_obs_seq <- c(100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000)

# Prior on causal structures (DMAG)
prior_DMAG <- uniform_prior_GRN_DMAG()

# List in which we will store the simulation results
Bayes_factors_consistency_reproduction <- list()
Bayes_factors_consistency_reproduction$binomial <- list()
Bayes_factors_consistency_reproduction$gaussian <- list()

if (Sys.info()['sysname'] == 'Windows') {
  num_cores <- 1
} else {
  num_cores <- parallel::detectCores() / 2
}


# a-c) (X_1, X_2, X_3) is multivariate normal data ------------------------


# Figure 3a - Causal model

seed_gaussian_causal <- 1226
set.seed(seed_gaussian_causal)

gamma_gaussian_causal <- rnorm(num_reps)
beta_gaussian_causal <- rnorm(num_reps)

Bayes_factors_consistency_reproduction$gaussian$causal <- lapply(num_obs_seq, function(num_obs) {
  print(paste("Number of observations:", num_obs))
  corr_mats <- do.call('cbind', pbmclapply(1:num_reps, function(i) {
    L <- rnorm(num_obs)
    T_i <- gamma_gaussian_causal[i] * L + rnorm(num_obs)
    T_j <- beta_gaussian_causal[i] * T_i + rnorm(num_obs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }, mc.cores = num_cores))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], num_samples = num_obs)
  
  t(Bfs)
})

# Figure 3b - Independent Model

seed_gaussian_independent <- 1500
set.seed(seed_gaussian_independent)

gamma_gaussian_independent <- rnorm(num_reps)
alpha_gaussian_independent <- rnorm(num_reps)

Bayes_factors_consistency_reproduction$gaussian$independent <- lapply(num_obs_seq, function(num_obs) {
  print(paste("Number of observations:", num_obs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:num_reps, function(i) {
    L <- rnorm(num_obs)
    T_i <- gamma_gaussian_independent[i] * L + rnorm(num_obs)
    T_j <- alpha_gaussian_independent[i] * L + rnorm(num_obs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }, mc.cores = num_cores))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], num_samples = num_obs)
  
  t(Bfs)
})

# Figure 3c - Full Model

seed_gaussian_full <- 1531
set.seed(seed_gaussian_full)

gamma_gaussian_full <- rnorm(num_reps)
alpha_gaussian_full <- rnorm(num_reps)
beta_gaussian_full <- rnorm(num_reps)

Bayes_factors_consistency_reproduction$gaussian$full <- lapply(num_obs_seq, function(num_obs) {
  print(paste("Number of observations:", num_obs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:num_reps, function(i) {
    
    L <- rnorm(num_obs)
    T_i <- gamma_gaussian_full[i] * L + rnorm(num_obs)
    T_j <- alpha_gaussian_full[i] * L + beta_gaussian_full[i] * T_i + rnorm(num_obs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }, mc.cores = num_cores))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], num_samples = num_obs)
  
  t(Bfs)
})


# d-f) X_1 is Bernoulli, (X_2, X_3) | X_1 is Gaussian ---------------------

# Figure 3d - Causal Model

seed_binomial_causal <- 1800
set.seed(seed_binomial_causal)

theta_binomial_causal <- runif(num_reps, 0.1, 0.5)
gamma_binomial_causal <- rnorm(num_reps)
beta_binomial_causal <- rnorm(num_reps)

Bayes_factors_consistency_reproduction$binomial$causal <- lapply(num_obs_seq, function(num_obs) {
  print(paste("Number of observations:", num_obs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:num_reps, function(i) {
    
    L <- rbinom(num_obs, 1, theta_binomial_causal[i]) + 1
    T_i <- gamma_binomial_causal[i] * L + rnorm(num_obs)
    T_j <- beta_binomial_causal[i] * T_i + rnorm(num_obs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }, mc.cores = num_cores))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], num_samples = num_obs)
  
  t(Bfs)
})

## Figure 3e - Independent Model

seed_binomial_independent <- 1611
set.seed(seed_binomial_independent)

theta_binomial_independent <- runif(num_reps, 0.1, 0.5)
gamma_binomial_independent <- rnorm(num_reps)
alpha_binomial_independent <- rnorm(num_reps)

Bayes_factors_consistency_reproduction$binomial$independent <- lapply(num_obs_seq, function(num_obs) {
  print(paste("Number of observations:", num_obs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:num_reps, function(i) {
    
    L <- rbinom(num_obs, 1, theta_binomial_independent[i]) + 1
    T_i <- gamma_binomial_independent[i] * L + rnorm(num_obs)
    T_j <- alpha_binomial_independent[i] * L + rnorm(num_obs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }, mc.cores = num_cores))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], num_samples = num_obs)
  
  t(Bfs)
})

## Figure 3f - Full Model

seed_binomial_full <- 1958
set.seed(seed_binomial_full)

theta_binomial_full <- runif(num_reps, 0.1, 0.5)
gamma_binomial_full <- rnorm(num_reps)
alpha_binomial_full <- rnorm(num_reps)
beta_binomial_full <- rnorm(num_reps)

Bayes_factors_consistency_reproduction$binomial$full <- lapply(num_obs_seq, function(num_obs) {
  print(paste("Number of observations:", num_obs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:num_reps, function(i) {
    
    L <- rbinom(num_obs, 1, theta_binomial_full[i]) + 1
    T_i <- gamma_binomial_full[i] * L + rnorm(num_obs)
    T_j <- alpha_binomial_full[i] * L + beta_binomial_full[i] * T_i + rnorm(num_obs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }, mc.cores = num_cores))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], num_samples = num_obs)
  
  t(Bfs)
})

## Save consistency evaluation results

save(Bayes_factors_consistency_reproduction, file = 'data/Bayes_factors_consistency_reproduction.RData')


# 4.2. Causal Discovery in Gene Regulatory Networks - Figure 4 + Table 2 PGM / Figures 4-7 IJAR ----

compute_simulated_GRN_probabilities <- function(
  GRN, trigger_num_runs = 3, trigger_window_size = 50000, trigger_runs_seed = 1200
) {
  
  # Create trigger object to compute local linkage
  trigger_obj <- trigger.build(marker = GRN$L, exp = GRN$Trait, 
                               marker.pos = GRN$L.pos, exp.pos = GRN$T.pos)
  
  # window.size 50000 is the value used in Chen et al.
  trigger_obj <- trigger.loclink(trigger_obj, window.size = trigger_window_size)
  
  results <- list()
  
  results[['BFCS DAG']] <- compute_BFCS_probabilities(GRN, prior_DAG)
  results[['BFCS DMAG']] <- compute_BFCS_probabilities(GRN, prior_DMAG)
  results[['BFCS loclink']] <- compute_BFCS_probabilities(GRN, prior_DMAG, loc_link = trigger_obj@loc.obj$loc.idx)
  results[['BGe']] <- compute_BGe_probabilities(GRN, prior_DMAG)
  
  set.seed(trigger_runs_seed)
  trigger_seeds <- sample.int(10000L, trigger_num_runs) 
  results[["trigger"]] <- matrix(0, nrow(GRN$Trait), nrow(GRN$Trait))
  
  for(i in 1:trigger_num_runs) {
    
    trigger_run <- paste("trigger", i)
    results[[trigger_run]] <- compute_trigger_probabilities(GRN, trigger_obj, trigger_seeds[i], window.size = trigger_window_size)
    
    # add results of individual run to average, take care that compute_trigger_probabilities may return NULL
    if (is.null(results[[trigger_run]]) || is.null(results[["trigger"]])) {
      results[["trigger"]] <- NULL
    } else {
      results[["trigger"]] <- results[["trigger"]] + results[[trigger_run]]
    }
    
  }
  
  # Divide by the number of runs, perform sanity check
  if (!is.null(results[["trigger"]])) {
    results[["trigger"]] <- results[["trigger"]] / trigger_num_runs
    stopifnot(all(0 <= results[["trigger"]] & results[["trigger"]] <= 1))
  }  
  
  # Set dimension names using the GRN info (e.g. T1-T100 for expression traits)
  results <- lapply(results, function(result) {
    dimnames(result) <- list(rownames(GRN$Trait), rownames(GRN$Trait))
    result
  })
  
  results
}

## Initialize simulation parameters
ngen <- 100 # number of gene markers
nexp <- 100 # number of gene expression levels
simulated_GRN_probabilities_reproduction <- list()

## Simulate GRNs and structure them in a list
simulated_GRN_list <- list(
  sparse_graph = list(
    samples_1e2 = simulate_GRN(ngen, nexp, 100, 1 / nexp, 1000),
    samples_1e3 = simulate_GRN(ngen, nexp, 1000, 1 / nexp, 1000)
  ),
  less_dense_graph = list(
    samples_1e2 = simulate_GRN(ngen, nexp, 100, 0.05, 1000),
    samples_1e3 = simulate_GRN(ngen, nexp, 1000, 0.05, 1000)
  ),
  more_dense_graph = list(
    samples_1e2 = simulate_GRN(ngen, nexp, 100, 0.1, 1000),
    samples_1e3 = simulate_GRN(ngen, nexp, 1000, 0.1, 1000)
  )
)

for (structure_type in c('sparse_graph', 'less_dense_graph', 'more_dense_graph')) {
  for (sample_size in c('samples_1e2', 'samples_1e3')) {
    print(paste("Computing probabilities for", structure_type, "from", sample_size, "observations..."))
    simulated_GRN_probabilities_reproduction[[structure_type]][[sample_size]] <- 
      compute_simulated_GRN_probabilities(simulated_GRN_list[[structure_type]][[sample_size]])
    simulated_GRN_probabilities_reproduction[[structure_type]][['structure']][['direct']] <- 
      simulated_GRN_list[[structure_type]][[sample_size]][['B']]
    simulated_GRN_probabilities_reproduction[[structure_type]][['structure']][['ancestral']] <- 
      as.matrix(Matrix::expm(simulated_GRN_list[[structure_type]][[sample_size]][['B']]) > 0) * 1
  }
}

# Save simulated GNR probabilities
save(simulated_GRN_probabilities_reproduction,
     file = 'data/simulated_GRN_probabilities_reproduction.RData')

# 4.3. Comparing Results from an Experiment on Yeast - Table 3 PGM / 2 IJAR ----

data(yeast)

subset <- 1:100

yeast_GRN <- list(
  L = yeast$marker[subset, ],
  Trait = yeast$exp[subset, ],
  L.pos = yeast$marker.pos[subset, ],
  T.pos = yeast$exp.pos[subset, ]
)

yeast_BFCS_DMAG <- compute_BFCS_probabilities(yeast_GRN)
yeast_trigger_w50k <- get_trigger_probabilities(yeast_GRN)

# 5. Computational and Time Complexity - Figure 8 and 9 IJAR ---------------



#' Function that runs three algorithms (BFCS, BGe, trigger) on a simulated GRN
#' and computes the time it takes for these to compute posterior probabilities.
#' 
#' @details WARNING: For small networks, trigger often does not manage to 
#' properly estimate pi0 and fails, in which case we return NA.
#'
#' @param GRN List containing GRN structure, as produced by simulate_GRN
#'
#' @return numeric vector containing time elapsed for each algorithm, as well
#' as the number of observations, variables (expression traits), and causal
#' structure edges for the simulated GRN
#' @export
#'
#' @examples
#' GRN <- simulate_GRN(ngen = 20, nexp = 50, nobs = 100, prob_edge = 0.05, seed = 1200)
#' simulate_GRN(GRN)
compare_time_on_simulated_GRN <- function(GRN) {
  
  # First call garbage collector
  gc()
  
  num_obs <- ncol(GRN$Trait)
  num_var <- nrow(GRN$Trait)
  
  print(paste("Number of observations:", num_obs, "; Number of variables:", num_var))
  
  time_BFCS <- system.time(compute_BFCS_probabilities(GRN))
  time_BGe <- system.time(compute_BGe_probabilities(GRN))
  
  # # We first need to build the trigger object and look for local linkage
  trigger_obj <- trigger::trigger.build(marker = GRN$L, exp = GRN$Trait,
                                        marker.pos = GRN$L.pos, exp.pos = GRN$T.pos)
  capture.output(trigger_obj <- trigger::trigger.loclink(trigger_obj))
  
  # We only consider the time it takes trigger to compute the probabilities
  time_trigger <- system.time(trigger::trigger_result <- compute_trigger_probabilities(GRN, trigger_obj, seed))
  if (is.null(trigger_result)) time_trigger['elapsed'] <- NA
  
  c(BFCS = time_BFCS['elapsed'],
    BGe = time_BGe['elapsed'],
    trigger = time_trigger['elapsed'],
    num_obs = num_obs, 
    num_var = num_var, 
    num_edges = GRN$num_edges)
}

num_var_seq <- seq(5, 20, 5)

time_versus_network_size_reproduction <- tibble::as_tibble(
  do.call('rbind', pbmclapply(num_var_seq, function(num_var) {
    
    GRN <- simulate_GRN(ngen = num_var, nexp = num_var, nobs = 100, prob_edge = 0.05, seed = 1200)
    
    compare_time_on_simulated_GRN(GRN)
  }, mc.cores = num_cores)))

num_obs_seq <- ceiling(10^(seq(2, 4, 0.25)))

time_versus_sample_size_reproduction <- tibble::as_tibble(
  do.call('rbind', pbmclapply(num_obs_seq, function(num_obs) {
    
    GRN <- simulate_GRN(ngen = 100, nexp = 100, nobs = num_obs, prob_edge = 0.05, seed = 1200)
    
    compare_time_on_simulated_GRN(GRN)
  }, mc.cores = num_cores)))

saveRDS(time_versus_network_size_reproduction, time_versus_sample_size_reproduction,
        file = 'data/computational_and_time_complexity_reproduction.RData')