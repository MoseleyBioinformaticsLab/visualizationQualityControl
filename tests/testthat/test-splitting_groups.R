context("test-splitting_groups")

data("grp_cor_data")

matrix_data = grp_cor_data$data
class_data = grp_cor_data$class


test_that("character and factor works", {
  
  split_character = split_groups(matrix_data, class_data)
  expect_equal(c("grp1" = 10, "grp2" = 10), sapply(split_character$groups, length))
  
  class_data2 = as.factor(class_data)
  split_factor = split_groups(matrix_data, class_data2)
  expect_equal(c("grp1" = 10, "grp2" = 10), sapply(split_factor$groups, length))
})

test_that("data.frame works", {
  class_df = data.frame(group = class_data, stringsAsFactors = FALSE)
  
  split_df = split_groups(matrix_data, class_df)
  expect_equal(c("grp1" = 10, "grp2" = 10), sapply(split_df$groups, length))
  
  class_df2 = data.frame(group = class_data, group2 = 1)
  split_df2 = split_groups(matrix_data, class_df2)
  expect_equal(c("grp1.1" = 10, "grp2.1" = 10), sapply(split_df2$groups, length))
})

test_that("rotating matrix works", {
  matrix2 = t(matrix_data)
  split_char2 = split_groups(matrix2, class_data)
  expect_equal(c("grp1" = 10, "grp2" = 10), sapply(split_char2$groups, length))
})
