#' Experiment on simulated data for PGM paper
#' @author Ioan Gabriel Bucur
# 0. SETUP ----------------------------------------------------------------

## Load libraries and scripts
source('R/compute_Bayes_factors.R')
source('R/compute_prior_structures.R')
library(pcalg)
library(pbmcapply)
library(trigger)

## Global variables
prior_DAG <- uniform_prior_GRN_DAG()
prior_DMAG <- uniform_prior_GRN_DMAG()


# Simulate GRNs -----------------------------------------------------------

## Initialize simulation parameters
set.seed(1020)
ngen <- 100 # number of gene markers
nexp <- ngen # number of gene expression levels

simulate_GRN <- function(ngen, nexp, nobs, prob_edge, seed = NULL) {
  
  set.seed(seed)
  
  # dense graph
  G <- pcalg:::randomDAG(nexp, prob_edge, lB = 1, uB = 1, V = paste0("T", 1:nexp))
  
  stopifnot(nexp >= ngen)
  
  ## Generate data
  
  theta <- runif(ngen, 0.1, 0.5) 
  
  L <- t(sapply(theta, rbinom, n = nobs, size = 1)) + 1
  rownames(L) <- paste0("L", 1:ngen)
  L.pos <- matrix(1, ngen, 2); L.pos[, 2] <- 1:ngen # one chromosome, each in the same position
  
  
  B <- t(as(G, "matrix")); print(paste("Number of edges:", sum(B)))
  C <- matrix(0, nexp, ngen); C[1:ngen, 1:ngen] <- diag(ngen)
  
  Tr <- solve(diag(nexp) - B) %*% (C %*% L + matrix(rnorm(nexp * nobs), nexp, nobs))
  # Tr <- t(apply(Tr, 1, function(row) row / sd(row))) # normalize to Var[X] = 1
  
  Trait <- Tr
  Trait <- t(apply(Tr, 1, function(x) qnorm( rank(x) / (length(x) + 1) ) ) ) # re-normalize data
  T.pos <- matrix(1, nexp, 3); T.pos[, 2:3] <- 1:nexp
  
  data <- t(rbind(L, Tr))
  cor.matrix <- cor(data) # unnormalized
  
  list(L = L, Trait = Trait, L.pos = L.pos, T.pos = T.pos, 
       cor.matrix = cor.matrix, B = B)
}

sparse_GRN_n1e2 <- simulate_GRN(ngen, nexp, 100, 1 / nexp, 1020)
sparse_GRN_n1e3 <- simulate_GRN(ngen, nexp, 1000, 1 / nexp, 1020)
dense_GRN_n1e2 <- simulate_GRN(ngen, nexp, 100, 0.1, 1020)
dense_GRN_n1e3 <- simulate_GRN(ngen, nexp, 1000, 0.1, 1020)


# 1. Run trigger on simulated GRNs ----------------------------------------

set.seed(1020)
trigger_seeds <- sample.int(10000L, 3) 

run_trigger_simulated_GRN <- function(GRN, sim_name, n_runs = 3, seeds = NULL) {
  
  trigger_obj <- trigger.build(marker = GRN$L, exp = GRN$Trait, 
                               marker.pos = GRN$L.pos, exp.pos = GRN$T.pos)
  
  # window.size 50kb is the value used in Chen et al.
  trigger_obj <- trigger.loclink(trigger_obj, window.size = 50000)
  
  for(i in 1:n_runs) {
    tryCatch({
      oldwd <- setwd("~/")
      trig_prob <- t(trigger.net(trigger_obj, seed = seeds[i]))
      write.table(trig_prob, file = paste0(oldwd, '/results/', sim_name, '_run', i, '.txt'))
    }, finally = {
      file.remove('~/net_trigg_prob.txt')
      setwd(oldwd)
    })
  }
}

run_trigger_simulated_GRN(sparse_GRN_n1e2, 'trigger_sparse_n1e2', seeds = trigger_seeds)
run_trigger_simulated_GRN(sparse_GRN_n1e3, 'trigger_sparse_n1e3', seeds = trigger_seeds)
run_trigger_simulated_GRN(dense_GRN_n1e2, 'trigger_dense_n1e2', seeds = trigger_seeds)
run_trigger_simulated_GRN(dense_GRN_n1e3, 'trigger_dense_n1e3', seeds = trigger_seeds)


# 2. Run BFCS on simulated GRNs -------------------------------------------

run_BFCS_simulated_GRN <- function(GRN, filename_root) {
  
  ngen <- nrow(GRN$L)
  nexp <- nrow(GRN$Trait)
  # eobs <- log10(ncol(GRN$L))
  
  cLT <- GRN$cor.matrix[1:ngen, (ngen+1):(ngen+nexp)]
  cTT <- GRN$cor.matrix[(ngen+1):(ngen+nexp), (ngen+1):(ngen+nexp)]
  
  BFCS_sim <- lapply(1:nexp, function(i) {
    compute_Bayes_factors_vectorized(rep(cLT[, i], nexp), c(cLT), rep(cTT[i, ], each = ngen),
                                     no_samples = ncol(GRN$L))
  })
  
  BFCS_prob_DAG <- sapply(1:nexp, function(i) {
    Bf <- prior_DAG * t(BFCS_sim[[i]])
    prob <- matrix(Bf[6, ] / colSums(Bf), ngen, nexp)
    prob <- apply(prob, 2, max)
    prob[i] <- 0
    prob
  })
  write.table(BFCS_prob_DAG, file = paste0('results/', filename_root, '_DAG.txt'))
  
  # print(paste("TPR:", sum(BFCS_prob_DAG * GRN$B) / sum(GRN$B)))
  # nB <- 1 - GRN$B; diag(nB) <- 0
  # print(paste("FPR:", sum(BFCS_prob_DAG * nB) / sum(nB)))
  
  BFCS_prob_DMAG <- sapply(1:nexp, function(i) {
    Bf <- prior_DMAG * t(BFCS_sim[[i]])
    prob <- matrix(Bf[6, ] / colSums(Bf), ngen, nexp)
    prob <- apply(prob, 2, max)
    prob[i] <- 0
    prob
  })
  write.table(BFCS_prob_DMAG, file = paste0('results/', filename_root, '_DMAG.txt'))
  
  # print(paste("TPR:", sum(BFCS_prob_DMAG * GRN$B) / sum(GRN$B)))
  # nB <- 1 - GRN$B; diag(nB) <- 0
  # print(paste("FPR:", sum(BFCS_prob_DMAG * nB) / sum(nB)))
}

run_BFCS_simulated_GRN(sparse_GRN_n1e2, 'BFCS_sparse_n1e2')
run_BFCS_simulated_GRN(sparse_GRN_n1e3, 'BFCS_sparse_n1e3')
run_BFCS_simulated_GRN(dense_GRN_n1e2, 'BFCS_dense_n1e2')
run_BFCS_simulated_GRN(dense_GRN_n1e3, 'BFCS_dense_n1e3')


# 3. Score local structures with BGe --------------------------------------

library(BiDAG)
# library(bnlearn)

# DAGs on triples (amat format, upper diagonal)

# 1. Empty (1 DAG w/ BK)
empty <- matrix(0, 3, 3)

# 2. Independent (4 DAG w/ BK)
indep_1 <- empty; indep_1[2, 3] <- 1 # alternatively indep_1[3, 2] <- 1
indep_2 <- empty; indep_2[1, 3] <- 1
indep_3 <- empty; indep_3[1, 2] <- 1

# 3. Causal (3 DAG w/ BK)
causal_123 <- indep_2; causal_123[1, 2] <- 1
causal_231 <- indep_3; causal_231[2, 3] <- 1
causal_312 <- indep_2; causal_312[3, 2] <- 1

# 4. Acausal (2 DAG w / BK)
acausal_12 <- indep_2; acausal_12[2, 3] <- 1
acausal_23 <- empty; acausal_23[2, 1] <- acausal_23[3, 1] <- 1
acausal_31 <- indep_3; acausal_31[3, 2] <- 1

# 5. Full (2 DAG w / BK)
full <- acausal_12; full[1, 2] <- 1
# alternatively full <- acausal_31; full[1, 3] <- 1 

# Markov equivalence classes with BK
mec_list <- list(
  empty, 
  indep_1,
  indep_2,
  indep_3,
  causal_123,
  causal_231,
  causal_312,
  acausal_12,
  acausal_23,
  acausal_31,
  full
)

run_BiDAG_simulated_GRN <- function(GRN, mec_list, results_root) {
  
  library(pbmcapply)
  
  ngen <- nrow(GRN$L)
  nexp <- nrow(GRN$Trait)
  eobs <- log10(ncol(GRN$L))

    
  result_BGe <- pbmclapply(1:nexp, function(i) {
    
    prob_mat <- matrix(0, nexp, 2)
    
    T_i <- GRN$Trait[i, ]
    
    for (j in 1:nexp) {
      
      if (i == j) next
      T_j <- GRN$Trait[j, ]
      
      probs <- sapply(1:ngen, function(k) {

        L_k <- GRN$L[k, ]
        
        triple <- data.frame(L_k, T_i, T_j)
        triple <- apply(triple, 2, scale) # normalize
        
        # Compute posterior probabilities
        lev <- sapply(mec_list, function(mec) {
          DAGscore(3, scoreparameters(3, "bge", triple), mec)
        })

        ev <- exp(lev - max(lev)) * prior_DAG
        ev <- ev / sum(ev) # normalize probabilities
        ev[6] # probability of L_k -> T_i -> T_j
      })
      
      prob_mat[j, 1] <- max(probs)
      prob_mat[j, 2] <- which.max(probs)
    }
    
    prob_mat
  })
  
  probs_BGe <- simplify2array(lapply(result_BGe, '[', , 1))
  ids_BGe <- simplify2array(lapply(result_BGe, '[', , 2))


  write.table(probs_BGe, file = paste0(results_root, '_n1e', eobs, '.txt'))
  
  print(paste("TPR:", sum(probs_BGe * GRN$B) / sum(GRN$B)))
  nB <- 1 - GRN$B; diag(nB) <- 0
  print(paste("FPR:", sum(probs_BGe * nB) / sum(nB)))
  
  list(probs = probs_BGe, ids = ids_BGe)
}

run_BiDAG_simulated_GRN(sparse_GRN_n1e2, mec_list, 'results/BGe_BiDAG_sG')
run_BiDAG_simulated_GRN(sparse_GRN_n1e3, mec_list, 'results/PGM/simulated_GRN/BiDAG_BGe_prob_sG')
run_BiDAG_simulated_GRN(dense_GRN_n1e2, mec_list, 'results/PGM/simulated_GRN/BiDAG_BGe_prob_dG')
run_BiDAG_simulated_GRN(dense_GRN_n1e3, mec_list, 'results/PGM/simulated_GRN/BiDAG_BGe_prob_dG')

run_BiDAG_simulated_GRN(sparse_GRN_n1e2, mec_list, 'results/PGM/simulated_GRN/BiDAG_BGe_prob_scaled_sG')
run_BiDAG_simulated_GRN(sparse_GRN_n1e3, mec_list, 'results/PGM/simulated_GRN/BiDAG_BGe_prob_scaled_sG')
run_BiDAG_simulated_GRN(dense_GRN_n1e2, mec_list, 'results/PGM/simulated_GRN/BiDAG_BGe_prob_scaled_dG')
run_BiDAG_simulated_GRN(dense_GRN_n1e3, mec_list, 'results/PGM/simulated_GRN/BiDAG_BGe_prob_scaled_dG')


# DAGs on triples (bnlearn format)

# 1. Empty (1 DAG w/ BK)
empty <- empty.graph(c("L_k", "T_i", "T_j"))
#empty <- empty.graph(names(data_triple))

# 2. Independent (4 DAG w/ BK)
indep_1_ij <- set.arc(empty, "T_i", "T_j")
indep_1_ji <- set.arc(empty, "T_j", "T_i")
indep_2 <- set.arc(empty, "L_k", "T_j")
indep_3 <- set.arc(empty, "L_k", "T_i")

# 3. Causal (3 DAG w/ BK)
causal_123 <- set.arc(indep_2, "L_k", "T_i")
causal_231 <- set.arc(indep_3, "T_i", "T_j")
causal_312 <- set.arc(indep_2, "T_j", "T_i")

# 4. Acausal (2 DAG w / BK)
acausal_12 <- set.arc(indep_2, "T_i", "T_j")
acausal_23 <- set.arc(set.arc(empty, "T_i", "L_k"), "T_j", "L_k")
acausal_31 <- set.arc(indep_3, "T_j", "T_i")

# 5. Full (2 DAG w / BK)
full_ij <- set.arc(acausal_12, "L_k", "T_i")
full_ji <- set.arc(acausal_31, "L_k", "T_j")

dag_list <- list(
  # 1. Empty (1 DAG w/ BK)
  empty,
  indep_1_ij,
  indep_1_ji,
  indep_2,
  indep_3,
  causal_123,
  causal_231,
  causal_312,
  acausal_12,
  acausal_31,
  full_ij,
  full_ji
)

mec_list <- list(
  empty,
  indep_1_ij,
  indep_2,
  indep_3,
  causal_123,
  causal_231,
  causal_312,
  acausal_12,
  acausal_23,
  acausal_31,
  full_ij
)

lev_bnlearn <- sapply(mec_list, score, data = data_triple, type = "bge", iss = 5)
lev_bidag <- sapply(mec_list, function(g) DAGscore(3, scoreparameters(3, "bge", data_triple, bgepar = list(am = 1, aw = 5)), amat(g)))

lev_bge <- sapply(mec_list, my_bge, data_triple = data_triple)

# lev_to_prob <- function(lev) {
#   ev <- exp(lev - max(lev)) * prior_DAG
#   ev / sum(ev)
# }


# 4. Run Trigger Reimplementation -----------------------------------------

# trigger.pl <- compute_primary_linkage(L, Trait)
# trigger.pl.loci <- apply(trigger.pl, 1, which.max)
# trigger.pl.prob <- apply(trigger.pl, 1, max)
# 
# 
# # trigger.pl <- compute_primary_linkage_stat(L, Trait, L.pos, T.pos, B = 5, seed = 1753)
# trigger.sl.ci <- pbmclapply(1:nexp, compute_sl_ci_combined, L = L, Trait = Trait, cis.locus = trigger.pl.loci, B = 50)
# 
# trigger.sl.ci.1 <- trigger.sl.ci[[1]]
# trigger.sl.ci.1.redo <- compute_sl_ci_combined(L, Trait, cis.locus = trigger.pl.loci, B = 50, seed = NULL)
# 
# trigger.sl.ci.alt <- pbmclapply(1:nexp, compute_sl_ci_combined, L = L, Trait = Trait, mc.cores = 1, cis.locus = prob.loc, B = 50)
# # trigger.sl <- sapply(1:nexp, function(i) compute_secondary_linkage_stat(L, Trait, trigger.pl.loci, i = i, B = 50))
# # trigger.ci <- sapply(1:nexp, function(i) compute_indep_test_stat(L, Trait, trigger.pl.loci, i = i, B = 50))
# timp_prob <- sapply(1:nexp, function(i) {
#   trigger.pl.prob[i] * trigger.sl.ci[[i]][,1] * trigger.sl.ci[[i]][,2]
# })
# 
# print(paste("TPR:", sum(timp_prob * B) / sum(B)))
# nB <- 1 - B; diag(nB) <- 0
# print(paste("FPR:", sum(timp_prob * nB) / sum(nB)))



# 5. Compare Results ------------------------------------------------------

# FDR, FPR, TPR
# thresholds <- seq(0, 1, 0.001)
# 
# BFCS_TPR <- sapply(thresholds, function(t) {
#   discoveries <- BFCS_prob > t
#   sum(discoveries * B) / sum(B)
# })
# trig_TPR <- sapply(thresholds, function(t) {
#   discoveries <- timp_prob > t
#   sum(discoveries * B) / sum(B)
# })
# plot(BFCS_TPR ~ thresholds, type = 'l', ylim = c(0, 1), col = 'blue', ylab = 'TPR', xlab = 'threshold')
# lines(trig_TPR ~ thresholds, col = 'red')
# legend("topright", col = c('blue', 'red'), lty = c(1, 1), legend = c("BFCS", "Trigger"))
# 
# BFCS_FPR <- sapply(thresholds, function(t) {
#   discoveries <- BFCS_prob > t
#   sum(discoveries * nB) / sum(nB)
# })
# trig_FPR <- sapply(thresholds, function(t) {
#   discoveries <- timp_prob > t
#   sum(discoveries * nB) / sum(nB)
# })
# plot(BFCS_FPR ~ thresholds, type = 'l', ylim = c(0, 1), col = 'blue', ylab = 'FPR', xlab = 'threshold')
# lines(trig_FPR ~ thresholds, col = 'red')
# legend("topright", col = c('blue', 'red'), lty = c(1, 1), legend = c("BFCS", "Trigger"))
# 
# 
# BFCS_FDR <- sapply(thresholds, function(t) {
#   discoveries <- BFCS_prob > t
#   ifelse(sum(discoveries) == 0, 0, sum(discoveries * nB) / sum(discoveries))
# })
# trig_FDR <- sapply(thresholds, function(t) {
#   discoveries <- timp_prob > t
#   ifelse(sum(discoveries) == 0, 0, sum(discoveries * nB) / sum(discoveries))
# })
# 
# plot(BFCS_FDR ~ thresholds, type = 'l', ylim = c(0, 1), col = 'blue', ylab = 'FDR', xlab = 'threshold')
# lines(trig_FDR ~ thresholds, col = 'red')
# legend("bottomleft", col = c('blue', 'red'), lty = c(1, 1), legend = c("BFCS", "Trigger"))
# 
# library(dplyr)
# BFCS_sorted_probs <- cbind(prob = c(BFCS_prob), true_reg = c(B), expand.grid(row = 1:nexp, column = 1:nexp)) %>%
#   filter(row != column) %>%
#   arrange(desc(prob)) %>%
#   mutate(rank = row_number()) %>%
#   mutate(fdr = cumsum(1 - true_reg) / rank)
# 
# trig_sorted_probs <- cbind(prob = c(trig_prob), true_reg = c(B), expand.grid(row = 1:nexp, column = 1:nexp)) %>%
#   filter(row != column) %>%
#   arrange(desc(prob)) %>%
#   mutate(rank = row_number()) %>%
#   mutate(fdr = cumsum(1 - true_reg) / rank)


# Slightly modified scoreparameters for BGe score, made compatible with BCCD
score_bge <- function (n,  data, weightvector = NULL, 
          bgepar = list(am = 1, aw = NULL), nodeslabels = NULL) 
{
  if (anyNA(data)) {
    stop("Dataset contains missing data")
  }
  if (ncol(data) != n) {
    stop("n and number of columns in data do not match")
  }
  if (!is.null(weightvector)) {
    if (length(weightvector) != nrow(data)) {
      stop("Length of the weightvector does not match the number of columns (observations) in data")
    }
  }
  if (is.null(nodeslabels)) {
    if (all(is.character(colnames(data)))) {
      nodeslabels <- colnames(data)
    }
    else {
      nodeslabels <- sapply(c(1:n), function(x) paste("v", 
                                                      x, sep = ""))
    }
  }
  colnames(data) <- nodeslabels
  initparam <- list()
  initparam$labels <- nodeslabels
  initparam$type <- "bge"
  initparam$weightvector <- weightvector
  

  initparam$data <- data
  if (is.null(weightvector)) {
    N <- nrow(data)
    covmat <- cov(data)
    means <- colMeans(data)
  }
  else {
    N = sum(weightvector)
    forcov <- cov.wt(data, wt = weightvector, cor = TRUE)
    covmat <- forcov$cov
    means <- forcov$center
  }
  if (is.null(bgepar$aw)) {
    bgepar$aw <- n + bgepar$am + 1
  }
  mu0 <- numeric(n)
  T0scale <- bgepar$am * (bgepar$aw - n - 1)/(bgepar$am + 
                                                1)
  T0 <- diag(T0scale, n, n)
  initparam$TN <- T0 + N * covmat + ((bgepar$am * 
                                              N)/(bgepar$am + N)) * (mu0 - means) %*% t(mu0 - means)
  initparam$awpN <- bgepar$aw + N
  constscorefact <- -(N/2) * log(pi) + (1/2) * log(bgepar$am/(bgepar$am + 
                                                                N))
  initparam$scoreconstvec <- numeric(n)
  for (j in 1:n) {
    awp <- bgepar$aw - n + j
    initparam$scoreconstvec[j] <- constscorefact - lgamma(awp/2) + 
      lgamma((awp + N)/2) + ((awp + j - 1)/2) * log(T0scale)
  }
 
  attr(initparam, "class") <- "scoreparameters"
  initparam
}


# Test the PC algorithm on the data
# pc.res <- pc(suffStat = list(C = cor.matrix, n = nobs), indepTest = gaussCItest, alpha = 0.01,
#              labels = c(paste0("L", 1:ngen), paste0("T", 1:nexp)))
