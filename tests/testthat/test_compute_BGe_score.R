context("BGe score")

test_that("The BGe score computed by compute_BGe_score is the same as the one
          obtained by using the BiDAG package", {
  
  num_samples <- 100
  data_triple <- data.frame(L_k = rnorm(num_samples), 
                            T_i = rnorm(num_samples), 
                            T_j = rnorm(num_samples))
  
  mec_bn <- get_Markov_equivalence_class_list(format = "bnlearn")
  mec_amat <- get_Markov_equivalence_class_list(format = "amat")
  
  log_evidence <- sapply(mec_bn, function(mec) {
    compute_BGe_score(num_samples, colMeans(data_triple), cov(data_triple),
                      get_parents_list(mec), 5, rep(0, 3), 1)
  })
  
  log_evidence_BiDAG <- sapply(mec_amat, function(mec) {
    BiDAG::DAGscore(BiDAG::scoreparameters(3, "bge", data_triple), mec)
  })
  
  expect_equal(log_evidence, log_evidence_BiDAG)
})
