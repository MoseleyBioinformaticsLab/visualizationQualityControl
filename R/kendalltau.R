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
  pairpoints = combn(length(x), 2)
  
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
  x_pairs = cbind(x[pairpoints[1, ]], x[pairpoints[2, ]])
  y_pairs = cbind(y[pairpoints[1, ]], y[pairpoints[2, ]])
  
  na_pairs = rowSums(is.na(x_pairs) > 0) | rowSums(is.na(y_pairs) > 0)
  
  basic_x_sign = sign(x_pairs[!na_pairs, 1] - x_pairs[!na_pairs, 2])
  basic_y_sign = sign(y_pairs[!na_pairs, 1] - y_pairs[!na_pairs, 2])
  
  basic_concordant_pair = sum((basic_x_sign * basic_y_sign) > 0)
  basic_discordant_pair = sum((basic_x_sign * basic_y_sign) < 0)
  
  na_concordant_pair = 0
  na_discordant_pair = 0
  
  
  if (sum(na_pairs) > 0) {
    x_na = x_pairs[na_pairs, ]
    y_na = y_pairs[na_pairs, ]
    
    both_x_na = (rowSums(is.na(x_na)) == 2)
    n_both_x_na = sum(both_x_na)
    if (n_both_x_na == 0) {
      n_both_x_na = 1
    }
    both_y_na = (rowSums(is.na(y_na)) == 2)
    n_both_y_na = sum(both_y_na)
    if (n_both_y_na == 0) {
      n_both_y_na = 1
    }
    
    both_na = both_x_na | both_y_na
    n_both_na = sum(both_na)
    
    if (perspective %in% "local") {
      n_pairs = n_pairs - n_both_na
    }
    
    x_na = x_na[!both_na, ]
    y_na = y_na[!both_na, ]
    
    if (nrow(x_na) > 0) {
      na_sign = rep(NA, nrow(x_na))
      
      # instead of row-wise, can we just do a great big | like this
      # (! is_na(x_i) and ! is_na(x_j) and x_i > x_j and ! is_na(y_i) and ! is_na(y_j) and y_i > y_j) | 
      
      for (irow in seq(1, nrow(x_na))) {
        #message(irow)
        if (is.na(x_na[irow, 1]) && is.na(y_na[irow, 1])) {
          na_sign[irow] = 1
        } else if (is.na(x_na[irow, 1]) && is.na(y_na[irow, 2])) {
          na_sign[irow] = -1
        } else if (is.na(x_na[irow, 2]) && is.na(y_na[irow, 1])) {
          na_sign[irow] = -1
        } else if (is.na(x_na[irow, 2]) && is.na(y_na[irow, 2])) {
          na_sign[irow] = 1
        } else if ((x_na[irow, 1] < x_na[irow, 2]) && (is.na(y_na[irow, 1]))) {
          na_sign[irow] = 1
        } else if ((x_na[irow, 1] < x_na[irow, 2]) && (is.na(y_na[irow, 2]))) {
          na_sign[irow] = -1
        } else if ((x_na[irow, 2] < x_na[irow, 1]) && (is.na(y_na[irow, 1]))) {
          na_sign[irow] = -1
        } else if ((x_na[irow, 2] < x_na[irow, 1]) && (is.na(y_na[irow, 2]))) {
          na_sign[irow] = 1
        } else if ((y_na[irow, 1] < y_na[irow, 2]) && (is.na(x_na[irow, 1]))) {
          na_sign[irow] = 1
        } else if ((y_na[irow, 1] < y_na[irow, 2]) && (is.na(x_na[irow, 2]))) {
          na_sign[irow] = -1
        } else if ((y_na[irow, 2] < y_na[irow, 1]) && (is.na(x_na[irow, 1]))) {
          na_sign[irow] = -1
        } else if ((y_na[irow, 2] < y_na[irow, 1]) && (is.na(x_na[irow, 2]))) {
          na_sign[irow] = 1
        }
      }
      na_concordant_pair = sum(na_sign == 1)
      na_discordant_pair = sum(na_sign == -1)
    } 
    
  } 
  
  
  sum_concordant = basic_concordant_pair + na_concordant_pair
  sum_discordant = basic_discordant_pair + na_discordant_pair
  
  k_tau = (sum_concordant - sum_discordant) / n_pairs
  return(k_tau)
}
