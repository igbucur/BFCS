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
nrep <- 1000
# Number of observations increases logarithmically for each set of simulations
nobs_seq <- c(100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000)

# Prior on causal structures (DMAG)
prior_DMAG <- uniform_prior_GRN_DMAG()

# create results subfolder if it is does not already exist
dir.create('results', showWarnings = FALSE)


# 1. (X_1, X_2, X_3) is multivariate normal data - (Figure 3a, 3b, 3c) ---------

# Figure 3a - Causal model ----

seed_T_norm <- 1226
set.seed(seed_T_norm)

gamma_T_norm <- rnorm(nrep)
beta_T_norm <- rnorm(nrep)

Bfs_seq_T_norm <- lapply(nobs_seq, function(nobs) {
  print(paste("No. Obs", nobs))
  corr_mats <- do.call('cbind', pbmclapply(1:nrep, function(i) {
    L <- rnorm(nobs)
    T_i <- gamma_T_norm[i] * L + rnorm(nobs)
    T_j <- beta_T_norm[i] * T_i + rnorm(nobs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], no_samples = nobs)
  
  t(Bfs)
})

print(sapply(Bfs_seq_T_norm, function(Bfs) {
  Bfs <- Bfs * prior_DMAG
  mean(Bfs[6, ] / colSums(Bfs))
}))


# Figure 3b - Independent Model ----

seed_I_norm <- 1500
set.seed(seed_I_norm)

gamma_I_norm <- rnorm(nrep)
alpha_I_norm <- rnorm(nrep)

Bfs_seq_I_norm <- lapply(nobs_seq, function(nobs) {
  print(paste("No. Obs", nobs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:nrep, function(i) {
    L <- rnorm(nobs)
    T_i <- gamma_I_norm[i] * L + rnorm(nobs)
    T_j <- alpha_I_norm[i] * L + rnorm(nobs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], no_samples = nobs)
  
  t(Bfs)
})

print(sapply(Bfs_seq_I_norm, function(Bfs) {
  Bfs <- Bfs * prior_DMAG
  mean(Bfs[6, ] / colSums(Bfs))
}))


# Figure 3c - Full Model ----

seed_F_norm <- 1531
set.seed(seed_F_norm)

gamma_F_norm <- rnorm(nrep)
alpha_F_norm <- rnorm(nrep)
beta_F_norm <- rnorm(nrep)

Bfs_seq_F_norm <- lapply(nobs_seq, function(nobs) {
  print(paste("No. Obs", nobs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:nrep, function(i) {
    
    L <- rnorm(nobs)
    T_i <- gamma_F_norm[i] * L + rnorm(nobs)
    T_j <- alpha_F_norm[i] * L + beta_F_norm[i] * T_i + rnorm(nobs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], no_samples = nobs)
  
  t(Bfs)
})

print(sapply(Bfs_seq_F_norm, function(Bfs) {
  Bfs <- Bfs * prior_DMAG
  mean(Bfs[6, ] / colSums(Bfs))
}))

# Save Consistency Evaluation Results ----

save(seed_T_norm, seed_I_norm, seed_F_norm,
     gamma_T_norm, gamma_I_norm, gamma_F_norm,
     beta_T_norm, beta_F_norm, alpha_I_norm, alpha_F_norm,
     Bfs_seq_T_norm, Bfs_seq_I_norm, Bfs_seq_F_norm,
     nobs_seq, nrep, file = 'results/consistency_norm.RData')


## 2. X_1 is Bernoulli, (X_2, X_3) | X_1 is Gaussian (Figure 3d, 3e, 3f) -------

# Figure 3d - Causal Model ----

seed_T_binom <- 1800
set.seed(seed_T_binom)

theta_T_binom <- runif(nrep, 0.1, 0.5)
gamma_T_binom <- rnorm(nrep)
beta_T_binom <- rnorm(nrep)

Bfs_seq_T_binom <- lapply(nobs_seq, function(nobs) {
  print(paste("No. Obs", nobs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:nrep, function(i) {
    
    L <- rbinom(nobs, 1, theta_T_binom[i]) + 1
    T_i <- gamma_T_binom[i] * L + rnorm(nobs)
    T_j <- beta_T_binom[i] * T_i + rnorm(nobs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], no_samples = nobs)
  
  t(Bfs)
})

print(sapply(Bfs_seq_T_binom, function(Bfs) {
  Bfs <- Bfs * prior_DMAG
  mean(Bfs[6, ] / colSums(Bfs))
}))

# Figure 3e - Independent Model ----

seed_I_binom <- 1611
set.seed(seed_I_binom)

theta_I_binom <- runif(nrep, 0.1, 0.5)
gamma_I_binom <- rnorm(nrep)
alpha_I_binom <- rnorm(nrep)

Bfs_seq_I_binom <- lapply(nobs_seq, function(nobs) {
  print(paste("No. Obs", nobs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:nrep, function(i) {
    
    L <- rbinom(nobs, 1, theta_I_binom[i]) + 1
    T_i <- gamma_I_binom[i] * L + rnorm(nobs)
    T_j <- alpha_I_binom[i] * L + rnorm(nobs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], no_samples = nobs)
  
  t(Bfs)
})

print(sapply(Bfs_seq_I_binom, function(Bfs) {
  Bfs <- Bfs * prior_DMAG
  mean(Bfs[6, ] / colSums(Bfs))
}))

# Figure 3f - Full Model ----

seed_F_binom <- 1958
set.seed(seed_F_binom)

theta_F_binom <- runif(nrep, 0.1, 0.5)
gamma_F_binom <- rnorm(nrep)
alpha_F_binom <- rnorm(nrep)
beta_F_binom <- rnorm(nrep)

Bfs_seq_F_binom <- lapply(nobs_seq, function(nobs) {
  print(paste("No. Obs", nobs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:nrep, function(i) {
    
    L <- rbinom(nobs, 1, theta_F_binom[i]) + 1
    T_i <- gamma_F_binom[i] * L + rnorm(nobs)
    T_j <- alpha_F_binom[i] * L + beta_F_binom[i] * T_i + rnorm(nobs)
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], no_samples = nobs)
  
  t(Bfs)
})

print(sapply(Bfs_seq_F_binom, function(Bfs) {
  Bfs <- Bfs * prior_DMAG
  mean(Bfs[6, ] / colSums(Bfs))
}))

# Save Consistency Evaluation Results ----

save(seed_T_binom, seed_I_binom, seed_F_binom,
     theta_T_binom, theta_I_binom, theta_F_binom,
     gamma_T_binom, gamma_I_binom, gamma_F_binom,
     beta_T_binom, beta_F_binom, alpha_I_binom, alpha_F_binom,
     Bfs_seq_T_binom, Bfs_seq_I_binom, Bfs_seq_F_binom,
     nobs_seq, nrep, file = 'results/consistency_binom.RData')

