data(grp_cor_data)
data_matrix = grp_cor_data$data
rownames(data_matrix) = paste0("f", seq(1, nrow(data_matrix)))
colnames(data_matrix) = paste0("s", seq(1, ncol(data_matrix)))

sample_info = data.frame(id = colnames(data_matrix), class = grp_cor_data$class)

log_pca = visqc_run_pca(t(log1p(data_matrix)))
pr_pca = prcomp(t(log1p(data_matrix)), scale. = FALSE)

test_that("the outputs are equivalent", {
  expect_equal(log_pca$scores |> as.matrix(), pr_pca$x)
  expect_equal(log_pca$loadings |> as.matrix(), pr_pca$rotation)
})

test_that("functions work with both matrix and data.frames", {
  test_info = sample_info[, "class", drop = FALSE]
  anova_matrix = visqc_test_pca_scores(
    pr_pca$x,
    test_info
  )
  anova_df = visqc_test_pca_scores(
    log_pca$scores,
    test_info
  )
  expect_equal(anova_matrix, anova_df)

  cor_matrix = visqc_cor_pca_scores(
    pr_pca$x,
    test_info
  )
  cor_df = visqc_cor_pca_scores(
    log_pca$scores,
    test_info
  )
  expect_equal(cor_matrix, cor_df)

  loading_matrix = visqc_test_pca_loadings(pr_pca$rotation, c("PC1", "PC2"))
  loading_df = visqc_test_pca_loadings(log_pca$loadings, c("PC1", "PC2"))
  expect_equal(loading_matrix$PC1$p.value, loading_df$PC1$p.value)

  variance_matrix = visqc_score_contributions(pr_pca$x)
  variance_df = log_pca$variances
  expect_equal(variance_matrix, variance_df)
})
