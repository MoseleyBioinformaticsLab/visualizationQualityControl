#' compute kendall tau
#' 
#' Given two vectors of data, computes the Kendall Tau correlation between them.
#' This version has logic for handling missing data in X and Y.
#' 
#' @param x vector of x data
#' @param y vector of y data
#' @param perspective how to treat missing data, see details
#' 
#' @return numeric
#' @export
#' 
#' @examples 
#' data("grp_cor_data")
#' exp_data = grp_cor_data$data
#' x = exp_data[, 1]
#' y = exp_data[, 2]
#' kendallt(x, y)
#' cor(x, y, method = "kendall") 
#' 
#' x = sort(rnorm(100))
#' y = x + 1
#' y2 = y
#' y2[1:10] = NA
#' kendallt(x, y)
#' kendallt(x, y2, "global")
#' kendallt(x, y2)
kendallt = function(x, y, perspective = "local"){
  if (length(x) != length(y)) {
    stop("x and y vector lengths are not the same!")
  }
  #pairpoints = combn(length(x), 2)
  
  # for local perspective
  # number of comparisons should be changed to (n * (n - 1)), this lets us modify n
  # when we have matching NA's in both x and y
  n = length(x)
  
  # if we don't do this, then they will get counted in the concordant pairs when they shouldn't
  # in the local version.
  # Note, we actually want to see these for the "global" version
  if (perspective %in% "local") {
    matching_na = (is.na(x) & is.na(y))
    n_matching_na = sum(matching_na)
    x = x[!matching_na]
    y = y[!matching_na]
  }
  
  # creates two matrices to hold the pairwise data in columnar format
  # x_i in column 1, x_j in column 2, and same for y
  x_index = t(combn(length(x), 2))
  y_index = t(combn(length(y), 2))
  
  x_pairs = matrix(x[c(x_index[, 1], x_index[, 2])], ncol = 2, byrow = FALSE)
  y_pairs = matrix(y[c(y_index[, 1], y_index[, 2])], ncol = 2, byrow = FALSE)
  
  # xi > xj and yi > yj                ## 1
  # xi < xj and yi < yj                ## 2
  # xi > xj and yi and not yj          ## 3
  # xi < xj and not yi and yj          ## 4
  # xi and not xj and yi > yj          ## 5
  # not xi and xj and yi < yj          ## 6
  # xi and not xj and yi and not yj    ## 7
  # not xi and xj and not yi and yj    ## 8
  
  
  concordant_pairs = 
    (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2])) |     #3 1
    (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) |   ## 2
    (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & !is.na(y_pairs[, 1]) & is.na(y_pairs[, 2])) |                              ## 3
    (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2])) |                              ## 4
    (!is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2])) |                              ## 5
    (is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) |                              ## 6
    (!is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & !is.na(y_pairs[, 1]) & is.na(y_pairs[, 2])) |                                                         ## 7
    (is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]))
  
  # xi > xj and yi < yj                 ## 1
  # xi < xj and yi > yj                 ## 2
  # xi > xj and not yi and yj           ## 3
  # xi < xj and yi and not yj           ## 4
  # xi and not xj and yi < yj           ## 5
  # not xi and xj and yi > yj           ## 6
  # xi and not xj and not yi and yj     ## 7
  # not xi and xj and yi and not yj     ## 8
  
  if (perspecitive == "global") {
    discordant_pairs = 
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) |
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2])) |
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2])) |
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & !is.na(y_pairs[, 1]) & is.na(y_pairs[, 2])) |
      (!is.na(x_pairs[, 1]) & is.na(x_pairs[, 2])  & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) |
      (is.na(x_pairs[, 1])  & !is.na(x_pairs[, 2])  & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2])) |
      (!is.na(x_pairs[, 1]) & is.na(x_pairs[, 2])  & is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2])) |
      (is.na(x_pairs[, 1])  & !is.na(x_pairs[, 2])  & !is.na(y_pairs[, 1]) & is.na(y_pairs[, 2]))
  } else {
    discordant_pairs = 
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) |
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2])) |
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2])) |
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & !is.na(y_pairs[, 1]) & is.na(y_pairs[, 2])) |
      (!is.na(x_pairs[, 1]) & is.na(x_pairs[, 2])  & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) |
      (is.na(x_pairs[, 1])  & !is.na(x_pairs[, 2])  & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2]))
  }
  
  
  x_ties = (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] == x_pairs[, 2]) & (!is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] != y_pairs[, 2]))) |
    (is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] != y_pairs[, 2]))
  
  y_ties = (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] != x_pairs[, 2]) & (!is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] == y_pairs[, 2]))) |
    (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] != x_pairs[, 2]) & is.na(y_pairs[, 1]) & is.na(y_pairs[, 2]))
  
  both_ties = (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] == x_pairs[, 2]) & (!is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] == y_pairs[, 2]))) |
    (is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & is.na(y_pairs[, 1]) & is.na(y_pairs[, 2]))
  
  x_na_ties = is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & (!is.na(y_pairs[, 1]) | !is.na(y_pairs[, 2]))
  sum_x_na_ties = sum(x_na_ties)
  y_na_ties = (!is.na(x_pairs[, 1]) | !is.na(x_pairs[, 2])) & is.na(y_pairs[, 1]) & is.na(y_pairs[, 2]) 
  sum_y_na_ties = sum(y_na_ties)
  
  all_na = is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & is.na(y_pairs[, 1]) & is.na(y_pairs[, 2])
  sum_all_na = sum(all_na) / 2

  
  sum_concordant = sum(concordant_pairs)
  sum_discordant = sum(discordant_pairs)
  
  if (perspective == "global") {
    sum_x_ties = sum(x_ties) + sum_x_na_ties + sum_all_na
    sum_y_ties = sum(y_ties) + sum_y_na_ties + sum_all_na
  } else {
    sum_x_ties = sum(x_ties)
    sum_y_ties = sum(y_ties)
  }
  
  log_multiplier = log(sum_concordant + sum_discordant + sum_x_ties) + log(sum_concordant + sum_discordant + sum_y_ties)
  log_tau = log((sum_concordant - sum_discordant)^2) - log_multiplier
  k_tau = exp(log_tau)
  return(k_tau)
}

#' information-theoretic kendall tau
#' 
#' Given a data-matrix, computes the information-theoretic Kendall-tau-b between
#' all samples.
#' 
#' @param data_matrix samples are rows, features are columns
#' @param exclude_na should NA values be treated as NA?
#' @param exclude_inf should Inf values be treated as NA?
#' @param exclude_0 should zero values be treated as NA?
#' @param zero_value what is the actual zero value?
#' @param perspective how to treat missing data in denominator and ties, see details
#' 
#' @return numeric
#' @export
#' 
visqc_it_kendallt = function(data_matrix, exclude_na = TRUE, exclude_inf = TRUE, 
                                    exclude_0 = TRUE, zero_value = 0, perspective = "local"){
  
  # assume row-wise (because that is what the description states), so need to transpose
  # because `cor` actually does things columnwise.
  data_matrix <- t(data_matrix)
  na_loc <- matrix(FALSE, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  inf_loc <- na_loc
  zero_loc <- na_loc
  
  if (exclude_na) {
    na_loc <- is.na(data_matrix)
  }
  
  if (exclude_inf) {
    inf_loc <- is.infinite(data_matrix)
  }
  
  if (exclude_0) {
    zero_loc <- data_matrix == zero_value
  }
  
  exclude_loc <- na_loc | zero_loc | inf_loc
  
  # set everything to NA and let R take care of it
  
  cor_matrix = matrix(NA, nrow = ncol(data_matrix), ncol = ncol(data_matrix))
  for (icol in seq(1, ncol(data_matrix) - 1)) {
    for (jcol in seq(2, ncol(data_matrix))) {
      cor_matrix[icol, jcol] = cor_matrix[jcol, icol] = kendallt(data_matrix[, icol], data_matrix[, jcol], perspective = perspective)
    }
  }
  calc_cor <- cor(data_matrix, use = use, method = method)
  
  return(list(cor = calc_cor, keep = t(!exclude_loc)))
}
