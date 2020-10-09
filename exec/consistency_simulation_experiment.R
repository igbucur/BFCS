#' @author Ioan Gabriel Bucur
#' @description This script performs the simulations described in Subsection 4.1
#' of the paper titled "A Bayesian Approach for Inferring Local Causal Structure
#' in Gene Regulatory Networks" and its extension "Large-Scale Local Causal 
#' Inference of Gene Regulatory Relationships". The purpose of this simulation 
#' is to empirically show that BFCS consistently recovers the correct local 
#' causal structure for multivariate Gaussian and conditionally Gaussian data.
#' @references 
#' \url{http://proceedings.mlr.press/v72/bucur18a.html}
#' \url{https://www.sciencedirect.com/science/article/abs/pii/S0888613X19301227}


# 0. Setup ----------------------------------------------------------------

library(pbmcapply)

source('R/compute_Bayes_factors.R')
source('R/compute_prior_structures.R')

# Number of replications in each simulation
num_reps <- 1000
# Number of observations sequence increases logarithmically for each set of simulations
num_obs_seq <- c(100, 300, 1000, 3000)
# num_obs_seq <- c(100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000)

# Prior on causal structures (DMAG)
prior_DMAG <- uniform_prior_GRN_DMAG()

# List in which we will store the simulation results
Bayes_factors_consistency <- list()
Bayes_factors_consistency$binomial <- list()
Bayes_factors_consistency$gaussian <- list()

if (Sys.info()['sysname'] == 'Windows') {
  num_cores <- 1
} else {
  num_cores <- parallel::detectCores() / 2
}

# 1. (X_1, X_2, X_3) is multivariate normal data - (Figure 3a, 3b, 3c) ---------

# Figure 3a - Causal model ----

seed_gaussian_causal <- 1226
set.seed(seed_gaussian_causal)

gamma_gaussian_causal <- rnorm(num_reps)
beta_gaussian_causal <- rnorm(num_reps)

Bayes_factors_consistency$gaussian$causal <- lapply(num_obs_seq, function(num_obs) {
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


# Figure 3b - Independent Model ----

seed_gaussian_independent <- 1500
set.seed(seed_gaussian_independent)

gamma_gaussian_independent <- rnorm(num_reps)
alpha_gaussian_independent <- rnorm(num_reps)

Bayes_factors_consistency$gaussian$independent <- lapply(num_obs_seq, function(num_obs) {
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


# Figure 3c - Full Model ----

seed_gaussian_full <- 1531
set.seed(seed_gaussian_full)

gamma_gaussian_full <- rnorm(num_reps)
alpha_gaussian_full <- rnorm(num_reps)
beta_gaussian_full <- rnorm(num_reps)

Bayes_factors_consistency$gaussian$full <- lapply(num_obs_seq, function(num_obs) {
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


## 2. X_1 is Bernoulli, (X_2, X_3) | X_1 is Gaussian (Figure 3d, 3e, 3f) -------

# Figure 3d - Causal Model ----

seed_binomial_causal <- 1800
set.seed(seed_binomial_causal)

theta_binomial_causal <- runif(num_reps, 0.1, 0.5)
gamma_binomial_causal <- rnorm(num_reps)
beta_binomial_causal <- rnorm(num_reps)

Bayes_factors_consistency$binomial$causal <- lapply(num_obs_seq, function(num_obs) {
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

# Figure 3e - Independent Model ----

seed_binomial_independent <- 1611
set.seed(seed_binomial_independent)

theta_binomial_independent <- runif(num_reps, 0.1, 0.5)
gamma_binomial_independent <- rnorm(num_reps)
alpha_binomial_independent <- rnorm(num_reps)

Bayes_factors_consistency$binomial$independent <- lapply(num_obs_seq, function(num_obs) {
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

# Figure 3f - Full Model ----

seed_binomial_full <- 1958
set.seed(seed_binomial_full)

theta_binomial_full <- runif(num_reps, 0.1, 0.5)
gamma_binomial_full <- rnorm(num_reps)
alpha_binomial_full <- rnorm(num_reps)
beta_binomial_full <- rnorm(num_reps)

Bayes_factors_consistency$binomial$full <- lapply(num_obs_seq, function(num_obs) {
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


# 3. Save Consistency Evaluation Results ----------------------------------

save(Bayes_factors_consistency, file = 'data/Bayes_factors_consistency_reproduced.RData')

