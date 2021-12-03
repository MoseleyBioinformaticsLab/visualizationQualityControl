test_that("defaults work", {
  col_1 <- generate_group_colors(3)
  col_2 <- generate_group_colors(3)
  
  expect_equal(col_1, col_2)
  
  col_3 <- generate_group_colors(5)
  col_4 <- generate_group_colors(5)
  
  expect_true(!(identical(col_3, col_4)))
})

test_that("overrides work", {
  col_1 <- generate_group_colors(3)
  col_2 <- generate_group_colors(3, TRUE)
  expect_true(!(identical(col_1, col_2)))
  
  col_3 <- generate_group_colors(5, FALSE)
  col_4 <- generate_group_colors(5, FALSE)
  expect_equal(col_3, col_4)
  
  col_5 <- generate_group_colors(5, TRUE)
  expect_true(!(identical(col_3, col_5)))
})

test_that("setting set.seed works", {
  set.seed(1234)
  col_1 <- generate_group_colors(5)
  
  set.seed(1234)
  col_2 <- generate_group_colors(5, TRUE)
  expect_equal(col_1, col_2)
})
