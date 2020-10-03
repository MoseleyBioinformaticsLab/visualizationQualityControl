run_long_kendallt <- as.logical(Sys.getenv("run_long_kendallt"))
if (is.na(run_long_kendallt)) {
  run_long_kendallt <- FALSE
}


test_that("basic kendall-tau matches base R", {
  x = seq(1, 10)
  y = seq(1, 10)
  expect_equal(ici_kendallt(x, y), cor(x, y, method = "kendall"))
  
  y[2] = 15
  expect_equal(ici_kendallt(x, y), cor(x, y, method = "kendall"))
  
  y = seq(10, 1) # should give -1
  expect_equal(ici_kendallt(x, y), cor(x, y, method = "kendall"))
  
  y[2] = 15
  expect_equal(ici_kendallt(x, y), cor(x, y, method = "kendall"))
})

test_that("difference and reference match - short", {
  x = sort(rnorm(100))
  y = x + 1
  y[1:20] = NA
  
  expect_equal(ici_kendallt(x, y, "global"), visualizationQualityControl:::ici_ref_kendallt(x, y, "global"))
})

test_that("matrix kendall works", {
  x = sort(rnorm(100))
  y = x + 1
  y[1:20] = NA
  
  test_mat = rbind(x, y)
  matrix_cor = visqc_ici_kendallt(test_mat, exclude_na = TRUE, exclude_0 = FALSE,
                                  perspective = "global", scale_max = FALSE)
  expect_equal(ici_kendallt(x, y, "global"), matrix_cor$cor[2, 1])
})

test_that("difference and reference match - long", {
  skip_if_not(run_long_kendallt)
  
  x = seq(1, 10)
  y = seq(1, 10)
  
  n_na = seq(1, 20)
  
  where_na = purrr::map(n_na, function(in_na){
    na_comb = combn(20, in_na)
    asplit(na_comb, 2)
  })
  where_na = unlist(where_na, recursive = FALSE)
  
  forward_na_diff = purrr::map_df(where_na, function(use_na){
    #message(.y)
    tmp_x = x
    tmp_y = y
    y_na = use_na[use_na > 10] - 10
    x_na = use_na[use_na <= 10]
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    data.frame(g = ici_kendallt(tmp_x, tmp_y, "global"), l = ici_kendallt(tmp_x, tmp_y, "local"))
  })
  
  y2 = seq(10, 1)
  reverse_na_diff = purrr::map_df(where_na, function(use_na){
    tmp_x = x
    tmp_y = y2
    y_na = use_na[use_na > 10] - 10
    y_na = 10 - y_na + 1
    x_na = use_na[use_na <= 10]
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    data.frame(g = ici_kendallt(tmp_x, tmp_y, "global"), l = ici_kendallt(tmp_x, tmp_y, "local"))
  })
  
  
  forward_na_old = purrr::map_df(where_na, function(use_na){
    tmp_x = x
    tmp_y = y
    y_na = use_na[use_na > 10] - 10
    x_na = use_na[use_na <= 10]
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    data.frame(g = visualizationQualityControl:::ici_ref_kendallt(tmp_x, tmp_y, "global"),
               l = visualizationQualityControl:::ici_ref_kendallt(tmp_x, tmp_y, "local"))
  })
  
  reverse_na_old = purrr::map_df(where_na, function(use_na){
    tmp_x = x
    tmp_y = y2
    y_na = use_na[use_na > 10] - 10
    y_na = 10 - y_na + 1
    x_na = use_na[use_na <= 10]
    tmp_y[y_na] = NA
    tmp_x[x_na] = NA
    data.frame(g = visualizationQualityControl:::ici_ref_kendallt(tmp_x, tmp_y, "global"),
               l = visualizationQualityControl:::ici_ref_kendallt(tmp_x, tmp_y, "local"))
  })
  
  expect_equal(forward_na_diff$g, forward_na_old$g)
  expect_equal(forward_na_diff$l, forward_na_old$l)
  
  expect_equal(reverse_na_diff$g, reverse_na_old$g)
  expect_equal(reverse_na_diff$l, reverse_na_old$l)
  
})
