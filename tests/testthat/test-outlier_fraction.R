data("grp_cor_data")
in_data <- t(grp_cor_data$data)
sample_classes <- grp_cor_data$class

test_that("no outlier fraction", {
  no_outlier <- vapply(seq(1, 10), function(x){rnorm(200, 0, 1)}, numeric(200))
  out_frac <- outlier_fraction(no_outlier)
  expect_equal(out_frac$frac, rep(0, 200))
  
})

test_that("inserted outlier fraction", {
  test_data <- vapply(seq(1, 2), function(x){rnorm(200, 0, 1)}, numeric(200))
  test_data[1:6,1] <- c(-8, -7, -6, 6, 7, 8)
  out_frac <- outlier_fraction(test_data)
  
  expect_equal(out_frac$frac, rep(c(0.5, 0), times = c(6, 194)))
})

test_that("zero argument works", {
  set.seed(456689)
  test_data <- vapply(seq(1, 2), function(x){rnorm(200, 5, 0.5)}, numeric(200))
  test_data[1:7,1] <- c(0, 0, 0, 9, 10, 11, 12)
  
  test_data[sample(seq(8, 200), 5), 1] <- 0
  zero_1 = outlier_fraction(test_data)
  expect_snapshot_value(zero_1, style = "serialize")
  
  zero_2 <- outlier_fraction(test_data, remove_0 = TRUE)
  expect_snapshot_value(zero_2, style = "serialize")
  
  out_1_2 <- identical(zero_1, zero_2)
  expect_false(out_1_2)
})

test_that("other arguments work", {
  set.seed(456689)
  test_data <- vapply(seq(1, 2), function(x){rnorm(200, 5, 1)}, numeric(200))
  test_data[1:6,1] <- c(0, 0, 0, 10, 11, 12)
  
  out_frac <- outlier_fraction(test_data)
  
  expect_equal(out_frac$frac, rep(c(0.5, 0), times = c(6, 194)))
  
  out_frac2 <- outlier_fraction(test_data, n_sd = 6)
  expect_equal(out_frac2$frac, rep(c(0, 0.5, 0), times = c(4, 2, 194)))
  
  out_frac3 <- outlier_fraction(test_data, n_trim = 10)
  expect_equal(out_frac3$frac, rep(c(0.5, 0), times = c(6, 194)))
})

test_that("double classes work", {
  data("grp_cor_data")
  in_data <- t(grp_cor_data$data)
  sample_classes <- grp_cor_data$class
  
  out_frac <- outlier_fraction(in_data)
  expect_snapshot_value(out_frac, style = "serialize")
  
  out_frac2 <- outlier_fraction(in_data, sample_classes)
  expect_snapshot_value(out_frac2, style = "serialize")
  
  out_frac3 <- outlier_fraction(in_data, as.factor(sample_classes))
  expect_equal(out_frac3, out_frac2)
  
  in_data2 <- in_data
  in_data2[3, ] <- in_data[12, ]
  in_data2[12, ] <- in_data[3, ]
  out_frac4 <- outlier_fraction(in_data2, sample_classes)
  expect_snapshot_value(out_frac4, style = "serialize")
})

test_that("too few entries handled properly", {
  test_data <- matrix(c(rep(0, 3), c(-8, -7, -6, 6)), nrow = 7, ncol = 1)
  out_frac <- outlier_fraction(test_data, remove_0 = TRUE)
  expect_equal(out_frac$frac, rep(0, 7))
})
