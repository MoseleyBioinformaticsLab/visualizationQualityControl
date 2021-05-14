#' compute kendall tau
#' 
#' Reference version for IT-kendall-tau. Given two vectors of data, computes the Kendall Tau correlation between them.
#' This version has logic for handling missing data in X and Y.
#' 
#' @param x vector of x data
#' @param y vector of y data
#' @param perspective how to treat missing data, see details
#' 
#' @keywords internal
#' 
#' @return numeric
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
ref_diff_kendallt = function(x, y, perspective = "local", alternative = "two.sided", output = "simple"){
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
  
  min_value = min(c(x, y), na.rm = TRUE)
  na_value = min_value - 0.1
  
  x[is.na(x)] = na_value
  y[is.na(y)] = na_value
  
  if (length(x) < 2) {
    return(NA)
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
  
  x_pair_sign = sign(x_pairs[, 2] - x_pairs[, 1])
  y_pair_sign = sign(y_pairs[, 2] - y_pairs[, 1])
  
  if (perspective == "global") {
    
    sum_concordant = sum((x_pair_sign * y_pair_sign) > 0)
    
    sum_discordant = sum((x_pair_sign * y_pair_sign) < 0)
    
  } else {
    sum_concordant = sum((x_pair_sign * y_pair_sign) > 0)
    
    sum_discordant = sum((x_pair_sign * y_pair_sign) < 0)
  }
  
  x_paired_tie = x_pairs[(x_pair_sign == 0) | (y_pair_sign == 0), ]
  y_paired_tie = y_pairs[(x_pair_sign == 0) | (y_pair_sign == 0), ]
  
  x_ties = (x_paired_tie[, 1] == x_paired_tie[, 2]) & 
    ((x_paired_tie[, 1] != na_value) & (x_paired_tie[, 2] != na_value)) &  
    ((y_paired_tie[, 1] != na_value) & (y_paired_tie[, 2] != na_value)) &
    (y_paired_tie[, 1] != y_paired_tie[, 2])
  
  y_ties = (y_paired_tie[, 1] == y_paired_tie[, 2]) &
    ((y_paired_tie[, 1] != na_value) & (y_paired_tie[, 2] != na_value)) & 
    ((x_paired_tie[, 1] != na_value) & (x_paired_tie[, 2] != na_value)) &
    (x_paired_tie[, 1] != x_paired_tie[, 2])
  
  
  x_na_ties = (x_pairs[, 1] == na_value) & (x_pairs[, 2] == na_value) & ((y_pairs[, 1] != na_value) | (y_pairs[, 2] != na_value))
  sum_x_na_ties = sum(x_na_ties)
  y_na_ties = ((x_paired_tie[, 1] != na_value) | (x_paired_tie[, 2] != na_value)) & (y_paired_tie[, 1] == na_value) & (y_paired_tie[, 2] == na_value) 
  sum_y_na_ties = sum(y_na_ties)
  
  all_na = (x_paired_tie[, 1] == na_value) & (x_paired_tie[, 2] == na_value) & (y_paired_tie[, 1] == na_value) & (y_paired_tie[, 2] == na_value)
  half_sum_na_ties = sum(all_na) / 2
  
  
  if (perspective == "global") {
    sum_x_ties = sum(x_ties) + sum_x_na_ties + half_sum_na_ties
    sum_y_ties = sum(y_ties) + sum_y_na_ties + half_sum_na_ties
  } else {
    sum_x_ties = sum(x_ties)
    sum_y_ties = sum(y_ties)
  }
  
  k_numerator = sum_concordant - sum_discordant
  k_denominator = sum_discordant + sum_concordant + sum_x_ties + sum_y_ties
  k_tau = k_numerator / k_denominator
  
  x_tied_values_t1 <- as.vector(table(x[duplicated(x)]) + 1)
  y_tied_values_t2 <- as.vector(table(y[duplicated(y)]) + 1)
  t_0 <- n * (n - 1)/2
  x_tied_sum_t1 <- sum(x_tied_values_t1 * (x_tied_values_t1 - 1))/2
  y_tied_sum_t2 <- sum(y_tied_values_t2 * (y_tied_values_t2 - 1))/2
  s_adjusted <- k_tau * sqrt((t_0 - x_tied_sum_t1) * (t_0 - y_tied_sum_t2))
  v_0_sum <- n * (n - 1) * sum(2 * n + 5)
  v_t_sum <- sum(x_tied_values_t1 * (x_tied_values_t1 - 1) * (2 * x_tied_values_t1 + 5))
  v_u_sum <- sum(y_tied_values_t2 * (y_tied_values_t2 - 1) * (2 * y_tied_values_t2 + 5))
  v_t1_sum <- sum(x_tied_values_t1 * (x_tied_values_t1 - 1)) * sum(y_tied_values_t2 * (y_tied_values_t2 - 1))
  v_t2_sum <- sum(x_tied_values_t1 * (x_tied_values_t1 - 1) * (x_tied_values_t1 - 2)) *
    sum(y_tied_values_t2 * (y_tied_values_t2 - 1) * (y_tied_values_t2 - 2))
  
  s_adjusted_variance <- (v_0_sum - v_t_sum - v_u_sum) / 18 +
    v_t1_sum / (2 * n * (n - 1)) +
    v_t2_sum / (9 * n * (n - 1) * (n - 2))
  s_adjusted <- sign(s_adjusted) * (abs(s_adjusted) - 1)
  z_b <- s_adjusted / sqrt(s_adjusted_variance)
  p_value <- switch(alternative,
                 "less" = pnorm(z_b),
                 "greater" = pnorm(z_b, lower.tail = FALSE),
                 "two.sided" = 2 * min(pnorm(z_b),
                                       pnorm(z_b, lower.tail = FALSE)))


  
  if (output == "simple") {
    return(c(tau = k_tau,
                pvalue = p_value))
  } else {
    out_data = data.frame(variable = c("n_entry",
                                       "x_ties",
                                       "x_na_ties",
                                       "x_na",
                                       "y_ties",
                                       "y_na_ties",
                                       "y_na",
                                       "half_sum_na_ties",
                                       "sum_concordant",
                                       "sum_discordant",
                                       "sum_numerator",
                                       "sum_denominator",
                                       "S",
                                       "k_tau",
                                       "pvalue"), 
                          value = c(length(x),
                                    sum(x_ties),
                                    sum_x_na_ties,
                                    sum((x == na_value)),
                                    sum(y_ties),
                                    sum_y_na_ties,
                                    sum((y == na_value)),
                                    half_sum_na_ties,
                                    sum_concordant,
                                    sum_discordant,
                                    k_numerator,
                                    k_denominator,
                                    s_adjusted,
                                    k_tau,
                                    p_value))
    return(out_data)
  }
  
}


#' compute kendall tau
#' 
#' Reference version for IT-kendall-tau. Given two vectors of data, computes the Kendall Tau correlation between them.
#' This version has logic for handling missing data in X and Y.
#' 
#' @param x vector of x data
#' @param y vector of y data
#' @param perspective how to treat missing data, see details
#' 
#' @keywords internal
#' 
#' @return numeric
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
ref_kendallt = function(x, y, perspective = "local", output = "simple"){
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
  
  if (length(x) < 2) {
    return(NA)
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
  
  if (perspective == "global") {
    concordant_pairs = 
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2])) |     #3 1
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) |   ## 2
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & !is.na(y_pairs[, 1]) & is.na(y_pairs[, 2])) |                              ## 3
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2])) |                              ## 4
      (!is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2])) |                              ## 5
      (is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) |                              ## 6
      (!is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & !is.na(y_pairs[, 1]) & is.na(y_pairs[, 2])) |                                                         ## 7
      (is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]))
    
  } else {
    concordant_pairs = 
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2])) |     #3 1
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) |   ## 2
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] > x_pairs[, 2]) & !is.na(y_pairs[, 1]) & is.na(y_pairs[, 2])) |                              ## 3
      (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] < x_pairs[, 2]) & is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2])) |                              ## 4
      (!is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] > y_pairs[, 2])) |                              ## 5
      (is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & !is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] < y_pairs[, 2])) 
    
  }
  
  
  
  # xi > xj and yi < yj                 ## 1
  # xi < xj and yi > yj                 ## 2
  # xi > xj and not yi and yj           ## 3
  # xi < xj and yi and not yj           ## 4
  # xi and not xj and yi < yj           ## 5
  # not xi and xj and yi > yj           ## 6
  # xi and not xj and not yi and yj     ## 7
  # not xi and xj and yi and not yj     ## 8
  
  if (perspective == "global") {
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
  
  
  x_ties = (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] == x_pairs[, 2]) & (!is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] != y_pairs[, 2]))) 
  
  y_ties = (!is.na(x_pairs[, 1]) & !is.na(x_pairs[, 2]) & (x_pairs[, 1] != x_pairs[, 2]) & (!is.na(y_pairs[, 1]) & !is.na(y_pairs[, 2]) & (y_pairs[, 1] == y_pairs[, 2]))) 
  
  
  x_na_ties = is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & (!is.na(y_pairs[, 1]) | !is.na(y_pairs[, 2]))
  sum_x_na_ties = sum(x_na_ties)
  y_na_ties = (!is.na(x_pairs[, 1]) | !is.na(x_pairs[, 2])) & is.na(y_pairs[, 1]) & is.na(y_pairs[, 2]) 
  sum_y_na_ties = sum(y_na_ties)
  
  all_na = is.na(x_pairs[, 1]) & is.na(x_pairs[, 2]) & is.na(y_pairs[, 1]) & is.na(y_pairs[, 2])
  half_sum_na_ties = sum(all_na) / 2

  
  sum_concordant = sum(concordant_pairs)
  sum_discordant = sum(discordant_pairs)
  
  if (perspective == "global") {
    sum_x_ties = sum(x_ties) + sum_x_na_ties + half_sum_na_ties
    sum_y_ties = sum(y_ties) + sum_y_na_ties + half_sum_na_ties
  } else {
    sum_x_ties = sum(x_ties)
    sum_y_ties = sum(y_ties)
  }
  
  k_numerator = sum_concordant - sum_discordant
  k_denominator = sum_discordant + sum_concordant + sum_x_ties + sum_y_ties
  k_tau = k_numerator / k_denominator
  

  if (output == "simple") {
    return(k_tau)
  } else {
    out_data = data.frame(variable = c("n_entry",
                                     "x_ties",
                                     "x_na_ties",
                                     "x_na",
                                     "y_ties",
                                     "y_na_ties",
                                     "y_na",
                                     "half_sum_na_ties",
                                     "sum_concordant",
                                     "sum_discordant",
                                     "sum_numerator",
                                     "sum_denominator",
                                     "k_tau"), 
                        value = c(length(x),
                                  sum(x_ties),
                                  sum_x_na_ties,
                                  sum(is.na(x)),
                                  sum(y_ties),
                                  sum_y_na_ties,
                                  sum(is.na(y)),
                                  half_sum_na_ties,
                                  sum_concordant,
                                  sum_discordant,
                                  k_numerator,
                                  k_denominator,
                                  k_tau))
    return(out_data)
  }
  
}

#' information-content-informed kendall tau
#' 
#' Given a data-matrix, computes the information-content-informed (ICI) Kendall-tau-b between
#' all samples.
#' 
#' @param data_matrix samples are rows, features are columns
#' @param exclude_na should NA values be treated as NA?
#' @param exclude_inf should Inf values be treated as NA?
#' @param exclude_0 should zero values be treated as NA?
#' @param zero_value what is the actual zero value?
#' @param perspective how to treat missing data in denominator and ties, see details
#' @param scale_max should everything be scaled compared to the maximum correlation?
#' @param diag_not_na should the diagonal entries reflect how many entries in the sample were "good"?
#' 
#' @details For more details, see the ICI-Kendall-tau vignette:
#'   \href{../doc/ici-kendalltau.html}{\code{vignette("ici-kendalltau", package = "visualizationQualityControl")}}
#' 
#' @return numeric
#' @internal
#' 
visqc_ici_kendallt_ref = function(data_matrix, 
                             exclude_na = TRUE, 
                             exclude_inf = TRUE, 
                             exclude_0 = TRUE, 
                             zero_value = 0, 
                             perspective = "global",
                             scale_max = TRUE,
                             diag_good = TRUE,
                             progress = FALSE){
  
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
  
  exclude_data = data_matrix
  exclude_data[exclude_loc] = NA
  # set everything to NA and let R take care of it
  
  if (ncol(data_matrix) > 2 && progress) {
    prog_bar = knitrProgressBar::progress_estimated(ncol(exclude_data) * (ncol(exclude_data))/ 2)
  } else {
    prog_bar = NULL
  }
  
  cor_matrix = matrix(NA, nrow = ncol(exclude_data), ncol = ncol(exclude_data))
  rownames(cor_matrix) = colnames(cor_matrix) = colnames(exclude_data)
  ntotal = 0
  for (icol in seq(1, ncol(exclude_data))) {
    for (jcol in seq(icol, ncol(exclude_data))) {
      cor_matrix[icol, jcol] = cor_matrix[jcol, icol] = ici_kendallt(exclude_data[, icol], exclude_data[, jcol], perspective = perspective)
      knitrProgressBar::update_progress(prog_bar)
      # ntotal = ntotal + 1
      # message(ntotal)
    }
  }
  
  # calculate the max-cor value for use in scaling across multiple comparisons
  n_observations = nrow(exclude_data)
  n_na = sort(colSums(exclude_loc))
  m_value = floor(sum(n_na[1:2]) / 2)
  n_m = n_observations - m_value
  max_cor_denominator = choose(n_m, 2) + n_observations * m_value
  max_cor_numerator = choose(n_m, 2) + n_observations * m_value + choose(m_value, 2)
  max_cor = max_cor_denominator / max_cor_numerator
  
  if (scale_max) {
    out_matrix = cor_matrix / max_cor
  } else {
    out_matrix = cor_matrix
  }
  
  if (diag_good) {
    n_good = colSums(!exclude_loc)
    diag(out_matrix) = n_good / max(n_good)
  }
  
  return(list(cor = out_matrix, raw = cor_matrix, keep = t(!exclude_loc)))
}


#' information-content-informed kendall tau
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
#' @param scale_max should everything be scaled compared to the maximum correlation?
#' @param diag_not_na should the diagonal entries reflect how many entries in the sample were "good"?
#' @param check_timings should we try to estimate run time for full dataset? (default is FALSE)
#' 
#' @details For more details, see the ICI-Kendall-tau vignette:
#'   \href{../doc/ici-kendalltau.html}{\code{vignette("ici-kendalltau", package = "visualizationQualityControl")}}
#'   
#'   When \code{check_timings = TRUE}, 5 random pairwise comparisons will be run to generate timings on a single core, and then estimates of how long the full set will take are calculated. The data is returned as a data.frame, and will be on the low side, but it should provide you with a good idea of how long your data will take.
#' 
#' @return numeric
#' @import furrr
#' @import future
#' @export
#' 
visqc_ici_kendallt = function(data_matrix, 
                             exclude_na = TRUE, 
                             exclude_inf = TRUE, 
                             exclude_0 = TRUE, 
                             zero_value = 0, 
                             perspective = "global",
                             scale_max = TRUE,
                             diag_good = TRUE,
                             check_timing = FALSE){
  
  # assume row-wise (because that is what the description states), so need to transpose
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
  
  exclude_data = data_matrix
  exclude_data[exclude_loc] = NA
  n_sample = ncol(exclude_data)
  # set everything to NA and let R take care of it
  
  ncore = future::nbrOfWorkers()
  names(ncore) = NULL
  
  pairwise_comparisons = combn(n_sample, 2)
  
  if (!diag_good) {
    extra_comparisons = matrix(rep(seq(1, n_sample), each = 2), nrow = 2, ncol = n_sample, byrow = FALSE)
    pairwise_comparisons = cbind(pairwise_comparisons, extra_comparisons)
  }
  
  n_todo = ncol(pairwise_comparisons)
  n_each = ceiling(n_todo / ncore)
  
  if (check_timing) {
    sample_compare = sample(ncol(pairwise_comparisons), 5)
    tmp_pairwise = pairwise_comparisons[, sample_compare]
    
    run_tmp = check_icikt_timing(exclude_data, tmp_pairwise, perspective, n_todo, ncore)
    return(run_tmp)
  }
  
  split_comparisons = vector("list", ncore)
  start_loc = 1
  
  for (isplit in seq_along(split_comparisons)) {
    stop_loc = min(start_loc + n_each, n_todo)
    
    split_comparisons[[isplit]] = pairwise_comparisons[, start_loc:stop_loc, drop = FALSE]
    start_loc = stop_loc + 1
    
    if (start_loc > n_todo) {
      break()
    }
  }
  
  null_comparisons = purrr::map_lgl(split_comparisons, is.null)
  split_comparisons = split_comparisons[!null_comparisons]
  
  do_split = function(do_comparisons, exclude_data, perspective) {
    #seq_range = seq(in_range[1], in_range[2])
    #print(seq_range)
    tmp_cor = matrix(0, nrow = ncol(exclude_data), ncol = ncol(exclude_data))
    rownames(tmp_cor) = colnames(tmp_cor) = colnames(exclude_data)
    tmp_pval = tmp_cor
    
    for (icol in seq(1, ncol(do_comparisons))) {
      iloc = do_comparisons[1, icol]
      jloc = do_comparisons[2, icol]
      ici_res = ici_kendallt(exclude_data[, iloc], exclude_data[, jloc], perspective = perspective)
      tmp_cor[iloc, jloc] = tmp_cor[jloc, iloc] = ici_res["tau"]
      tmp_pval[iloc, jloc] = tmp_pval[jloc, iloc] = ici_res["pvalue"]
    }
    list(cor = tmp_cor, pval = tmp_pval)
  }
  # we record how much time is actually spent doing ICI-Kt
  # itself, as some of the other operations will add a bit of time
  # 
  t1 = Sys.time()
  split_cor = furrr::future_map(split_comparisons, do_split, exclude_data, perspective)
  t2 = Sys.time()
  t_diff = as.numeric(difftime(t2, t1, units = "secs"))

  cor_matrix = matrix(0, nrow = ncol(exclude_data), ncol = ncol(exclude_data))
  rownames(cor_matrix) = colnames(cor_matrix) = colnames(exclude_data)
  pvalue_matrix = cor_matrix
  for (isplit in split_cor) {
    cor_matrix = cor_matrix + isplit$cor
    pvalue_matrix = pvalue_matrix + isplit$pval
  }
  
  
  # calculate the max-cor value for use in scaling across multiple comparisons
  n_observations = nrow(exclude_data)
  n_na = sort(colSums(exclude_loc))
  m_value = floor(sum(n_na[1:2]) / 2)
  n_m = n_observations - m_value
  max_cor_denominator = choose(n_m, 2) + n_observations * m_value
  max_cor_numerator = choose(n_m, 2) + n_observations * m_value + choose(m_value, 2)
  max_cor = max_cor_denominator / max_cor_numerator
  
  if (scale_max) {
    out_matrix = cor_matrix / max_cor
  } else {
    out_matrix = cor_matrix
  }
  
  if (diag_good) {
    n_good = colSums(!exclude_loc)
    diag(out_matrix) = n_good / max(n_good)
  }
  
  return(list(cor = out_matrix, raw = cor_matrix, pval = pvalue_matrix, keep = t(!exclude_loc),
              run_time = t_diff))
}


check_icikt_timing = function(exclude_data, tmp_pairwise, perspective, n_todo, ncore){
  t_start = Sys.time()
  for (icol in seq_len(ncol(tmp_pairwise))) {
    iloc = tmp_pairwise[1, icol]
    jloc = tmp_pairwise[2, icol]
    tmp_val = ici_kendallt(exclude_data[, iloc], exclude_data[, jloc], perspective = perspective)
  }
  t_stop = Sys.time()
  t_total = as.numeric(difftime(t_stop, t_start, units = "secs"))
  n_comp = ncol(tmp_pairwise)
  
  t_each = t_total / n_comp
  
  t_theoretical = t_each * n_todo
  t_cores = t_theoretical / ncore
  
  data.frame(which =
               c("n_tested",
                 "n_todo",
                 "time_tested",
                 "time_single",
                 "time_all",
                 "time_across_cores",
                 "time_minutes",
                 "time_hours",
                 "time_days"),
             value = c(n_comp,
                n_todo,
                t_total,
                t_each,
                t_theoretical,
                t_cores,
                t_cores / 60,
                t_cores / (60 * 60),
                t_cores / (60 * 60 * 60)))
  
}
