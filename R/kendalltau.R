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
  
  matching_na = (is.na(x) & is.na(y))
  n_matching_na = sum(matching_na)
  x = x[!matching_na]
  y = y[!matching_na]
  
  if (perspective %in% "local") {
    n = length(x) - (n_matching_na)
  }
  
  n_pairs = (n * (n - 1)) / 2
  
  # creates two matrices to hold the pairwise data in columnar format
  # x_i in column 1, x_j in column 2, and same for y
  # FIX THIS!!! It turns out NA's are bad for combn, need to use indices first and then make it
  x_pairs = t(combn(x, 2))
  y_pairs = t(combn(y, 2))
  
  # this does not currently test for ties, those should be added using the 
  # Tau-B 
  # https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient#Tau-b
  

  # instead of row-wise, can we just do a great big | like this
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
  
  discordant_pairs = 
    (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) |
    (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2])) |
    (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2])) |
    (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & !is.na(y_pairs[, 1]) & is.na(y_pairs[, 2])) |
    (!is.na(x_pairs[, 1]) & is.na(x_pairs[, 2])  & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) |
    (is.na(x_pairs[, 1])  & !is.na(x_pairs[, 2])  & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2])) |
    (!is.na(x_pairs[, 1]) & is.na(x_pairs[, 2])  & is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2])) |
    (is.na(x_pairs[, 1])  & !is.na(x_pairs[, 2])  & !is.na(y_pairs[, 1]) & is.na(y_pairs[, 2]))
  
  x_ties = (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] == x_pairs[, 2]) & (!is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] != y_pairs[, 2]))) |
    (is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] != y_pairs[, 2]))
  
  y_ties = (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] != x_pairs[, 2]) & (!is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] == y_pairs[, 2]))) |
    (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] != x_pairs[, 2]) & is.na(y_pairs[, 1]) & is.na(y_pairs[, 2]))
  
  both_ties = (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] == x_pairs[, 2]) & (!is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] == y_pairs[, 2]))) |
    (is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & is.na(y_pairs[, 1]) & is.na(y_pairs[, 2]))

  
  sum_concordant = sum(concordant_pairs)
  sum_discordant = sum(discordant_pairs)
  sum_x_ties = sum(x_ties)
  sum_y_ties = sum(y_ties)
  log_multiplier = log(sum_concordant + sum_discordant + sum_x_ties) + log(sum_concordant + sum_discordant + sum_y_ties)
  log_tau = log((sum_concordant - sum_discordant)^2) - log_multiplier
  k_tau = exp(log_tau)
  return(k_tau)
}
