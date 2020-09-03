context("basic_correlations")

data("grp_cor_data")
exp_data <- grp_cor_data$data
rownames(exp_data) <- paste0("f", seq(1, nrow(exp_data)))
colnames(exp_data) <- paste0("s", seq(1, ncol(exp_data)))

sample_info <- data.frame(id = colnames(exp_data), class = grp_cor_data$class)

sample_classes <- sample_info$class

test_that("correlations don't change", {
  data_cor <- pairwise_correlation(t(exp_data))
  expect_equal_to_reference(data_cor, "data_cor1")
})

test_that("adding missing data works", {
  exp_data2 <- exp_data
  set.seed(1234)
  make_na <- rep(FALSE, nrow(exp_data2))
  s1_missing <- make_na
  s1_missing[sample(length(make_na), 20)] <- TRUE
  s2_missing <- make_na
  s2_missing[sample(which(!s1_missing), 20)] <- TRUE
  
  exp_data2[s1_missing, 1] <- NA
  exp_data2[s2_missing, 2] <- NA
  
  data_cor <- pairwise_correlation(t(exp_data))
  data_cor2 <- pairwise_correlation(t(exp_data2))
  
  # because we added missing values in samples 1 & 2, we expect that their correlations
  # should be less than they were before.
  expect_lt(data_cor2$cor[1,2], data_cor$cor[1,2])
  expect_equal_to_reference(data_cor2, "data_cor2")
  expect_equal(data_cor2$cor[1,1], 1)
  expect_equal(data_cor2$keep[1,1], FALSE)
  # we didn't add any missing data to the other samples, so their correlation
  # should remain unchanged.
  expect_equal(data_cor2$cor[3:ncol(exp_data),3:ncol(exp_data)], data_cor$cor[3:ncol(exp_data),3:ncol(exp_data)])
})
