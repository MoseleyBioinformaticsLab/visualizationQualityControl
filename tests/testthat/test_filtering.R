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
          f_1 <- t(keep_non_zero_percentage(t(test_data), keep_num = 5))
          expect_equal(c("1", "2", "3", "5", "6", "7", "9", "10"), rownames(f_1))
          })
  

test_that("integer filtering works with classes", {
          f_2 <- t(keep_non_zero_percentage(t(test_data), sample_classes, 3))
          expect_equal(c("1", "2", "3", "4", "5", "6", "7", "9", "10"), rownames(f_2))
          })

test_that("percentage filtering works without classes", {
          f_3 <- t(keep_non_zero_percentage(t(test_data), keep_num = 0.6))
          expect_equal(c("1", "2", "3", "5", "6", "7", "9", "10"), rownames(f_3))
          })

test_that("percentage filtering works with classes", {
          f_4 <- t(keep_non_zero_percentage(t(test_data), sample_classes, 0.7))
          expect_equal(c("1", "2", "3", "4", "5", "6", "7", "9", "10"), rownames(f_4))
          })

test_that("error handling works", {
  expect_error(keep_non_zero_percentage(t(test_data), sample_classes, -1))
  expect_error(keep_non_zero_percentage(t(test_data), sample_classes, 100))
})
