#' compute kendall tau
#' 
#' Given two vectors of data, computes the Kendall Tau correlation between them.
#' This version has logic for handling missing data in X and Y.
#' 
#' @param x vector of x data
#' @param y vector of y data
#' @param remove_both_na should NA values in both be removed? (default is FALSE)
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
#' x[10] = NA
#' x[15] = NA
#' y[20] = NA
kendallt = function(x, y, remove_both_na = FALSE){
  if (length(x) != length(y)) {
    stop("x and y vector lengths are not the same!")
  }
  pairpoints = combn(length(x), 2)
  
  x_pairs = cbind(x[pairpoints[1, ]], x[pairpoints[2, ]])
  y_pairs = cbind(y[pairpoints[1, ]], y[pairpoints[2, ]])
  
  na_pairs = rowSums(is.na(x_pairs) > 0) | rowSums(is.na(y_pairs) > 0)
  
  basic_x_sign = sign(x_pairs[!na_pairs, 1] - x_pairs[!na_pairs, 2])
  basic_y_sign = sign(y_pairs[!na_pairs, 1] - y_pairs[!na_pairs, 2])
  
  basic_concordant_pair = sum((basic_x_sign * basic_y_sign) > 0)
  basic_discordant_pair = sum((basic_x_sign * basic_y_sign) < 0)
  
  sum_both_na = 0
  
  na_concordant_pair = 0
  na_discordant_pair = 0
  
  
  if (sum(na_pairs) > 0) {
    x_na = x_pairs[na_pairs, ]
    y_na = y_pairs[na_pairs, ]
    
    both_na = (rowSums(is.na(x_na)) == 2) | (rowSums(is.na(y_na))== 2)
    
    x_na = x_na[!both_na, ]
    y_na = y_na[!both_na, ]
    
    if (nrow(x_na) > 0) {
      na_sign = rep(NA, nrow(x_na))
      
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
    
    if (remove_both_na) {
      sum_both_na = sum(both_na) 
    }
  } 
  
  
  sum_concordant = basic_concordant_pair + na_concordant_pair
  sum_discordant = basic_discordant_pair + na_discordant_pair
  n_pairs = ncol(pairpoints) - sum_both_na
  
  k_tau = (sum_concordant - sum_discordant) / n_pairs
  return(k_tau)
}
