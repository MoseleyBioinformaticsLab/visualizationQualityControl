library(visualizationQualityControl)

x = seq(1, 10)
y = seq(1, 10)

n_na = seq(1, 8)

set.seed(1234)
select_na = purrr::map(n_na, function(use_na){
  sample(x, use_na)
})

forward_na = purrr::map_dbl(select_na, function(use_na){
  tmp_y = y
  tmp_y[use_na] = NA
  kendallt(x, tmp_y, "global")
})

y2 = seq(10, 1)
reverse_na = purrr::map_dbl(select_na, function(use_na){
  tmp_y = y2
  tmp_y[10 - use_na + 1] = NA
  kendallt(x, tmp_y, "global")
})

cbind(forward_na, reverse_na)
