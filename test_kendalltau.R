library(visualizationQualityControl)
library(ggplot2)
library(dplyr)

x = seq(1, 10)
y = seq(1, 10)

n_na = seq(1, 20)

where_na = purrr::map(n_na, function(in_na){
  na_comb = combn(20, in_na)
  asplit(na_comb, 2)
})
where_na = unlist(where_na, recursive = FALSE)

prog_bar = knitrProgressBar::progress_estimated(length(where_na))

forward_na = purrr::map_dbl(where_na, function(use_na){
  #message(.y)
  knitrProgressBar::update_progress(prog_bar)
  tmp_x = x
  tmp_y = y
  y_na = use_na[use_na > 10] - 10
  x_na = use_na[use_na <= 10]
  tmp_y[y_na] = NA
  tmp_x[x_na] = NA
  ici_kendallt(tmp_x, tmp_y, "global")
})

prog_bar = knitrProgressBar::progress_estimated(length(where_na))
y2 = seq(10, 1)
reverse_na = purrr::map_dbl(where_na, function(use_na){
  knitrProgressBar::update_progress(prog_bar)
  tmp_x = x
  tmp_y = y2
  y_na = use_na[use_na > 10] - 10
  y_na = 10 - y_na + 1
  x_na = use_na[use_na <= 10]
  tmp_y[y_na] = NA
  tmp_x[x_na] = NA
  ici_kendallt(tmp_x, tmp_y, "global")
})

all_na = data.frame(positive = forward_na, negative = reverse_na)
all_na = dplyr::mutate(all_na, diff = -1 * negative - positive)

ggplot(dplyr::slice_sample(all_na, n = 1000), aes(x = positive, y = negative)) + geom_point()

zero_diff = dplyr::filter(all_na, diff == 0)

long_na = tidyr::pivot_longer(all_na, -diff, names_to = "type", values_to = "correlation")

ggplot(long_na, aes(x = correlation)) + geom_histogram(bins = 100) + 
  facet_wrap(~ type, ncol = 1)

dplyr::group_by(long_na, type) %>%
  dplyr::summarize(mean = abs(mean(correlation)))

# just for kicks, lets now test the pearson correlation
prog_bar = knitrProgressBar::progress_estimated(length(where_na))

positive_pearson = purrr::map_dbl(where_na, function(use_na){
  #message(.y)
  knitrProgressBar::update_progress(prog_bar)
  tmp_x = x
  tmp_y = y
  y_na = use_na[use_na > 10] - 10
  x_na = use_na[use_na <= 10]
  tmp_y[y_na] = NA
  tmp_x[x_na] = NA
  in_matrix = rbind(tmp_x, tmp_y)
  out_res = locally_it_weighted_pairwise_correlation(in_matrix)
  out_res$cor[1, 2]
})

prog_bar = knitrProgressBar::progress_estimated(length(where_na))
y2 = seq(10, 1)
negative_pearson = purrr::map_dbl(where_na, function(use_na){
  knitrProgressBar::update_progress(prog_bar)
  tmp_x = x
  tmp_y = y2
  y_na = use_na[use_na > 10] - 10
  y_na = 10 - y_na + 1
  x_na = use_na[use_na <= 10]
  tmp_y[y_na] = NA
  tmp_x[x_na] = NA
  in_matrix = rbind(tmp_x, tmp_y)
  out_res = locally_it_weighted_pairwise_correlation(in_matrix)
  out_res$cor[1, 2]
})


pearson_cor = data.frame(correlation = c(positive_pearson,
                                         negative_pearson),
                         type = rep(c("positive", "negative"), each = length(positive_pearson)))
ggplot(pearson_cor, aes(x = correlation)) + geom_histogram(bins = 100) + 
  facet_wrap(~ type, ncol = 1)

dplyr::group_by(pearson_cor, type) %>%
  dplyr::summarise(mean = abs(mean(correlation, na.rm = TRUE)))

pearson_wide = data.frame(positive = positive_pearson, negative = negative_pearson)
pearson_wide = dplyr::mutate(pearson_wide, diff = positive - -1 * negative)
pearson_wide = dplyr::filter(pearson_wide, !is.na(negative), !is.na(positive))
ggplot(dplyr::slice_sample(pearson_wide, n = 10000), aes(x = positive, y = negative)) + geom_point()
