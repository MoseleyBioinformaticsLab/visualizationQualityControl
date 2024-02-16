# set up the data
test_data <- matrix(rnorm(80, mean = 20), nrow = 10, ncol = 8)
rownames(test_data) <- as.character(seq(1, nrow(test_data)))
test_data[1, c(1, 3, 6)] <- 0
test_data[3, c(1,2,3)] <- 0
test_data[4, c(1,2,3,4)] <- 0

test_data[7, 7] <- 0
test_data[8, c(1, 2, 6, 8)] <- 0

sample_classes <- rep(c("A", "B"), each = 4)


test_that("integer filtering works without classes", {
          f_1 <- keep_non_zero_percentage(test_data, keep_num = 5)
          expect_equal(c("1", "2", "3", "5", "6", "7", "9", "10"), names(f_1)[f_1])
          })
  

test_that("integer filtering works with classes", {
          f_2 <- keep_non_zero_percentage(test_data, sample_classes, 3)
          expect_equal(c("1", "2", "3", "4", "5", "6", "7", "9", "10"), names(f_2)[f_2])
          })

test_that("percentage filtering works without classes", {
          f_3 <- keep_non_zero_percentage(test_data, keep_num = 0.6)
          expect_equal(c("1", "2", "3", "5", "6", "7", "9", "10"), names(f_3)[f_3])
          })

test_that("percentage filtering works with classes", {
          f_4 <- keep_non_zero_percentage(test_data, sample_classes, 0.7)
          expect_equal(c("1", "2", "3", "4", "5", "6", "7", "9", "10"), names(f_4)[f_4])
          })

test_that("non-zero error handling works", {
  expect_error(keep_non_zero_percentage(t(test_data), sample_classes, -1))
  expect_error(keep_non_zero_percentage(t(test_data), sample_classes, 100))
})

test_missing = test_data
test_missing[test_missing == 0] = NA

test_that("missing: integer filtering works without classes", {
  f_1 <- keep_non_missing_percentage(test_missing, keep_num = 5)
  expect_equal(c("1", "2", "3", "5", "6", "7", "9", "10"), names(f_1)[f_1])
})


test_that("missing: integer filtering works with classes", {
  f_2 <- keep_non_missing_percentage(test_missing, sample_classes, 3)
  expect_equal(c("1", "2", "3", "4", "5", "6", "7", "9", "10"), names(f_2)[f_2])
})

test_that("missing: percentage filtering works without classes", {
  f_3 <- keep_non_missing_percentage(test_missing, keep_num = 0.6)
  expect_equal(c("1", "2", "3", "5", "6", "7", "9", "10"), names(f_3)[f_3])
})

test_that("missing: percentage filtering works with classes", {
  f_4 <- keep_non_missing_percentage(test_missing, sample_classes, 0.7)
  expect_equal(c("1", "2", "3", "4", "5", "6", "7", "9", "10"), names(f_4)[f_4])
})

test_norows = test_missing
rownames(test_norows) = NULL
test_that("missing: null rownames still works", {
  f_5 <- keep_non_missing_percentage(test_missing, sample_classes, 0.7)
  expect_equal(c("1", "2", "3", "4", "5", "6", "7", "9", "10"), names(f_5)[f_5])
})

test_na0 = test_missing
test_na0[8, 1] = 0

test_that("missing: mix of NA and 0 works", {
  f_6 <- keep_non_missing_percentage(test_missing, sample_classes, 0.7, missing_value = c(NA, 0))
  expect_equal(c("1", "2", "3", "4", "5", "6", "7", "9", "10"), names(f_6)[f_6])
})

test_that("missing: error handling works", {
  expect_error(keep_non_missing_percentage(test_missing, sample_classes, -1))
  expect_error(keep_non_zero_percentage(test_missing, sample_classes, 100))
})
