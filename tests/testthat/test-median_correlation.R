data(grp_cor_data)

grp1 <- grp_cor_data$data[, 1:10]
grp2 <- grp_cor_data$data[, 11:20]

grp1_class <- grp_cor_data$class[1:10]
grp_class <- grp_cor_data$class

test_that("median single works", {
  grp1_cor <- cor(grp1)
  expect_snapshot_value(median_correlations(grp1_cor), style = "serialize")
  grp1_class <- rep("grp1", 10)
  expect_snapshot_value(median_correlations(grp1_cor, grp1_class), style = "serialize")
  grp1_class <- factor(grp1_class)
  expect_snapshot_value(median_correlations(grp1_cor, grp1_class), style = "serialize")
})

test_that("median double works", {
  grp1_grp2 <- cor(cbind(grp1, grp2))
  
  
  expect_snapshot_value(median_correlations(grp1_grp2, grp_class), style = "serialize")
  grp_class <- factor(grp_class)
  expect_snapshot_value(median_correlations(grp1_grp2, grp_class), style = "serialize")

  expect_snapshot_value(median_correlations(grp1_grp2, grp_class, between_classes = TRUE), style = "serialize")
})

test_that("median swapped works", {
  all_samples <- cbind(grp1, grp2)
  all_samples[, 3] <- grp2[, 3]
  all_samples[, 15] <- grp1[, 3]
  
  all_cor <- cor(all_samples)
  
  expect_snapshot_value(median_correlations(all_cor, grp_class), style = "serialize")
})


test_that("median with rownames works", {
  all_samples <- cbind(grp1, grp2)
  all_samples[, 3] <- grp2[, 3]
  all_samples[, 15] <- grp1[, 3]
  
  all_cor <- cor(all_samples)
  
  rownames(all_cor) <- colnames(all_cor) <- paste0("s", seq(1, 20))
  
  expect_snapshot_value(median_correlations(all_cor, grp_class), style = "serialize")
})
