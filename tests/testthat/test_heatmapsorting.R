test_that("reordering works without rownames", {
  data(grp_cor)
  data(grp_info)
  
  rownames(grp_info) <- NULL
  rownames(grp_cor) <- colnames(grp_cor) <- NULL
  
  
  expected_order <- c(1,3,2,4,5,7,9,6,8,10)
  new_order <- similarity_reorderbyclass(grp_cor, grp_info, "sub_1")
  expect_equal(expected_order, new_order$indices)
  expect_null(new_order$names)
  expect_equal(as.character(expected_order), labels(new_order$dendrogram))
})

test_that("reordering works with rownames", {
  data(grp_cor)
  data(grp_info)
  expected_order <- c(1,3,2,4,5,7,9,6,8,10)
  new_order <- similarity_reorderbyclass(grp_cor, grp_info, "sub_1")
  expect_equal(expected_order, new_order$indices)
  expect_equal(new_order$names, rownames(grp_cor)[new_order$indices])
  expect_equal(new_order$names, labels(new_order$dendrogram))
})

test_that("changing rownames throws error", {
  data(grp_cor)
  data(grp_info)
  
  rownames(grp_cor) <- paste("c", seq(1,10), sep = "")
  expect_error(similarity_reorderbyclass(grp_cor, grp_info))
})

test_that("changing matrix changes sorting", {
  data(grp_cor)
  data(grp_info)
  
  tmp_order <- c(1,9,2,3,4,5,6,7,8,10)
  grp_cor <- grp_cor[tmp_order, tmp_order]
  grp_info <- grp_info[tmp_order,]
  
  expected_order <- c(1,4,3,5,6,2,8,7,9,10)
  new_order <- similarity_reorderbyclass(grp_cor, grp_info, "sub_1")
  expect_equal(expected_order, new_order$indices)
})

test_that("null sample_classes works", {
  data(grp_cor)
  expected_order <- c(7,9,6,8,10,1,5,4,2,3)
  new_order <- similarity_reorderbyclass(grp_cor, transform = "sub_1")
  expect_equal(expected_order, new_order$indices)
})

test_that("factor works as well as data.frame", {
  data(grp_cor)
  sample_classes <- as.factor(rep("1", 10))
  expected_order <- c(7,9,6,8,10,1,5,4,2,3)
  new_order <- similarity_reorderbyclass(grp_cor, sample_classes, transform = "sub_1")
  expect_equal(expected_order, new_order$indices)
})
