#' compute kendall tau
#' 
#' Given two vectors of data, computes the Kendall Tau correlation between them.
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
kendallt = function(x, y, remove_both_na = FALSE){
  if (length(x) != length(y)) {
    stop("x and y vector lengths are not the same!")
  }
  pairpoints = combn(length(x), 2)
  
  x_pairs = cbind(x[pairpoints[1, ]], x[pairpoints[2, ]])
  y_pairs = cbind(y[pairpoints[1, ]], y[pairpoints[2, ]])
  
  x_sign = sign(x_pairs[, 1] - x_pairs[, 2])
  y_sign = sign(y_pairs[, 1] - y_pairs[, 2])
  
  concordant_pair = (x_sign * y_sign) > 0
  discordant_pair = (x_sign * y_sign) < 0
  
  k_tau = (sum(concordant_pair) - sum(discordant_pair)) / ncol(pairpoints)
}
