#' @author Ioan Gabriel Bucur
#' @description This script performs the simulations described in Subsection 4.1
#' of the paper titled "A Bayesian Approach for Inferring Local Causal
#' Structure in Gene Regulatory Networks". The purpose of this simulation is
#' to empirically show that BFCS consistently recovers the correct local
#' causal structure for multivariate Gaussian and conditionally Gaussian data.
#' @references \url{http://proceedings.mlr.press/v72/bucur18a.html}


# 0. SETUP ----------------------------------------------------------------
library(pbmcapply)
source('R/compute_Bayes_factors.R')
source('R/compute_prior_structures.R')

nobs_seq <- c(100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000)
nrep <- 1000

prior_DAG <- uniform_prior_GRN_DAG()
prior_DMAG <- uniform_prior_GRN_DMAG()


# 1a.  Multivariate normal data, true (causal) model -----------------------------

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
  Bfs <- Bfs * prior_DAG
  mean(Bfs[6, ] / colSums(Bfs))
}))

# 1b. Multivariate normal data, independent model --------------------------------

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
  Bfs <- Bfs * prior_DAG
  mean(Bfs[6, ] / colSums(Bfs))
}))

# 1c. Multivariate normal data, full model --------------------------------

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
  Bfs <- Bfs * prior_DAG
  mean(Bfs[6, ] / colSums(Bfs))
}))


# 1d. Save data, multivariate normal --------------------------------------

save(seed_T_norm, seed_I_norm, seed_F_norm,
     gamma_T_norm, gamma_I_norm, gamma_F_norm,
     beta_T_norm, beta_F_norm, alpha_I_norm, alpha_F_norm,
     Bfs_seq_T_norm, Bfs_seq_I_norm, Bfs_seq_F_norm,
     nobs_seq, nrep, file = '~/CHiLL/results/PGM/BFCS_consistency_norm_fixparam.RData')


# 2a. X_1 binomial distributed, causal model ------------------------------

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
  Bfs <- Bfs * prior_DAG
  mean(Bfs[6, ] / colSums(Bfs))
}))



# 2b. Binomial X_1, independent model -------------------------------------

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
  Bfs <- Bfs * prior_DAG
  mean(Bfs[6, ] / colSums(Bfs))
}))

# 2c. Binomial X_1, full model --------------------------------------------

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
  Bfs <- Bfs * prior_DAG
  mean(Bfs[6, ] / colSums(Bfs))
}))


# 2d. Save data, Binomial L -----------------------------------------------

save(seed_T_binom, seed_I_binom, seed_F_binom,
     theta_T_binom, theta_I_binom, theta_F_binom,
     gamma_T_binom, gamma_I_binom, gamma_F_binom,
     beta_T_binom, beta_F_binom, alpha_I_binom, alpha_F_binom,
     Bfs_seq_T_binom, Bfs_seq_I_binom, Bfs_seq_F_binom,
     nobs_seq, nrep, file = '~/CHiLL/results/PGM/BFCS_consistency_binom_fixparam.RData')


# 3a. Binomial X_1, causal model, inverse rank transformation for T -------

seed_T_irank <- 1842
set.seed(seed_T_irank)

theta_T_irank <- runif(nrep, 0.1, 0.5)
gamma_T_irank <- rnorm(nrep)
beta_T_irank <- rnorm(nrep)


Bfs_seq_T_irank <- lapply(nobs_seq, function(nobs) {
  print(paste("No. Obs", nobs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:nrep, function(i) {
    
    L <- rbinom(nobs, 1, theta_T_irank[i]) + 1
    T_i <- gamma_T_irank[i] * L + rnorm(nobs)
    T_j <- beta_T_irank[i] * T_i + rnorm(nobs)
    
    T_i <- qnorm(rank(T_i) / (nobs + 1)) # re-normalize data
    T_j <- qnorm(rank(T_j) / (nobs + 1)) # re-normalize data
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], no_samples = nobs)
  
  t(Bfs)
})


print(sapply(Bfs_seq_T_irank, function(Bfs) {
  Bfs <- Bfs * prior_DAG
  mean(Bfs[6, ] / colSums(Bfs))
}))


# 3b. Binomial X_1, indep model, inverse rank transformation for T --------

seed_I_irank <- 1924
set.seed(seed_I_irank)

theta_I_irank <- runif(nrep, 0.1, 0.5)
gamma_I_irank <- rnorm(nrep)
alpha_I_irank <- rnorm(nrep)

Bfs_seq_I_irank <- lapply(nobs_seq, function(nobs) {
  print(paste("No. Obs", nobs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:nrep, function(i) {
    
    L <- rbinom(nobs, 1, theta_I_irank[i]) + 1
    T_i <- gamma_I_irank[i] * L + rnorm(nobs)
    T_j <- alpha_I_irank[i] * L + rnorm(nobs)
    
    T_i <- qnorm(rank(T_i) / (nobs + 1)) # re-normalize data
    T_j <- qnorm(rank(T_j) / (nobs + 1)) # re-normalize data
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], no_samples = nobs)
  
  t(Bfs)
})

print(sapply(Bfs_seq_I_irank, function(Bfs) {
  Bfs <- Bfs * prior_DAG
  mean(Bfs[6, ] / colSums(Bfs))
}))


# 3c. Binomial X_1, full model, inverse rank transformation for T ---------

seed_F_irank <- 1453
set.seed(seed_F_irank)

theta_F_irank <- runif(nrep, 0.1, 0.5)
gamma_F_irank <- rnorm(nrep)
alpha_F_irank <- rnorm(nrep)
beta_F_irank <- rnorm(nrep)

Bfs_seq_F_irank <- lapply(nobs_seq, function(nobs) {
  print(paste("No. Obs", nobs))
  
  corr_mats <- do.call('cbind', pbmclapply(1:nrep, function(i) {
    
    L <- rbinom(nobs, 1, theta_F_irank[i]) + 1
    T_i <- gamma_F_irank[i] * L + rnorm(nobs)
    T_j <- alpha_F_irank[i] * L + beta_F_irank[i] * T_i + rnorm(nobs)
    
    T_i <- qnorm(rank(T_i) / (nobs + 1)) # re-normalize data
    T_j <- qnorm(rank(T_j) / (nobs + 1)) # re-normalize data
    
    cor <- cor(cbind(L, T_i, T_j))
    cor[upper.tri(cor)]
  }))
  
  Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], no_samples = nobs)
  
  t(Bfs)
})

print(sapply(Bfs_seq_F_irank, function(Bfs) {
  Bfs <- Bfs * prior_DAG
  mean(Bfs[6, ] / colSums(Bfs))
}))


# 3d. Save data, binomial L and inverse rank transformation on T ----------

save(seed_T_irank, seed_I_irank, seed_F_irank,
     theta_T_irank, theta_I_irank, theta_F_irank,
     gamma_T_irank, gamma_I_irank, gamma_F_irank,
     beta_T_irank, beta_F_irank, alpha_I_irank, alpha_F_irank,
     Bfs_seq_T_irank, Bfs_seq_I_irank, Bfs_seq_F_irank,
     nobs_seq, nrep, file = '~/CHiLL/results/PGM/BFCS_consistency_irank_fixparam.RData')


# 4. Use Kendall / Spearman correlation -----------------------------------
# NOTE: does not help

# set.seed(1401)
# 
# corr_mats <- replicate(nrep, {
#   theta <- runif(1, 0.1, 0.5)
#   
#   gamma <- rnorm(1)
#   beta <- rnorm(1)
#   L <- rbinom(nobs, 1, theta) + 1
#   T_i <- gamma * L + rnorm(nobs)
#   T_j <- beta * T_i + rnorm(nobs)
#   
#   cor <- cor(cbind(L, T_i, T_j), method = "kendall")
#   cor[upper.tri(cor)]
# })
# 
# 
# Bfs <- compute_Bayes_factors_vectorized(corr_mats[1, ], corr_mats[2, ], corr_mats[3, ], no_samples = nobs)
# 
# # Bfs_check <- apply(data_sets, 3, function(data) {
# #   compute_Bayes_factors(cor(data), no_samples = nobs)
# # })
# 
# Bfs_DAG <- prior_DAG * t(Bfs)
# Bfs_DMAG <- prior_DMAG * t(Bfs)
# 
# print(mean(Bfs_DAG[6, ] / colSums(Bfs_DAG)))
# print(mean(Bfs_DMAG[6, ] / colSums(Bfs_DMAG)))
# 
# 
# 
# pl_pval <- function(L, Tr, B = 50, seed = 123) {
#   
#   stat <- link.stat.c(Tr, L)
#   
#   stat0 <- sapply(1:B, function(i) {
#     link.stat.c(sample(Tr), L)
#   })
#   
#   #microbenchmark(link.stat.xx.c(t(replicate(B, sample(Tr), simplify = TRUE)), rbind(L, L)))
#   #prob <- 1 - edge.lfdr(edge.pvalue(stat, stat0), pi0 = 1)
#   
#   edge.pvalue(stat, stat0)
# }

# results <- apply(results, 2, function(p) p / sum(p))
# table(apply(results, 2, which.max))
# apply(results, 1, mean)

# 1. Generate from L->T_i->T_j --------------------------------------------

# set.seed(1800)
# results <- replicate(nrep, {
#   theta <- runif(1, 0.1, 0.5)
#   gamma <- rnorm(1)
#   beta <- rnorm(1)
#   L <- rbinom(nobs, 1, theta)
#   T_i <- gamma * L + rnorm(nobs)
#   T_j <- beta * T_i + rnorm(nobs)
#   
#   compute_Bayes_factors(cor(cbind(L, T_i, T_j)), no_samples = nobs) * prior_GRN
# })
# results <- apply(results, 2, function(p) p / sum(p))
# table(apply(results, 2, which.max))
# apply(results, 1, mean)


# NOTE: needs too much memory
# data_sets <- replicate(nrep, {
#   # theta <- runif(1, 0.1, 0.5)
#   gamma <- rnorm(1)
#   beta <- rnorm(1)
#   L <- rbinom(nobs, 1, 0.3) + 1
#   T_i <- gamma * L + rnorm(nobs)
#   T_j <- beta * T_i + rnorm(nobs)
#   
#   cbind(L, T_i, T_j)
# 
#   # compute_Bayes_factors(cor(cbind(L, T_i, T_j)), no_samples = nobs) * prior_GRN
# })

# correlations <- apply(data_sets, 3, cor)
# Bfs <- compute_Bayes_factors_vectorized(correlations[2, ], correlations[3, ], correlations[6, ], no_samples = nobs)