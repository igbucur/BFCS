context("BFCS")

test_that("compute_BFCS_vectorized is consistent with compute_BFCS", {
  
  corrs <- apply(rWishart(100, 4, diag(3)), 3, cov2cor)
  
  for (num_samples in c(100, 300, 1000, 3000, 10000, 30000, 10000, 30000, 100000)) {
    Bfs_vec <- compute_BFCS_vectorized(corrs[2, ], corrs[3, ], corrs[6, ], num_samples)
    Bfs_iter <- t(apply(corrs, 2, function(corr) {
      compute_BFCS(matrix(corr, 3, 3), num_samples)
    }))
    
    expect_equivalent(Bfs_iter, Bfs_vec)
  }
  
})