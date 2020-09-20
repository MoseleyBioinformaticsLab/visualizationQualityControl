#' compute kendall tau
#' 
#' Reference version for IT-kendall-tau. Given two vectors of data, computes the Kendall Tau correlation between them.
#' This version has logic for handling missing data in X and Y.
#' 
#' @param x vector of x data
#' @param y vector of y data
#' @param perspective how to treat missing data, see details
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
#' @param scale_max should everything be scaled compared to the maximum correlation?
#' @param diag_not_na should the diagonal entries reflect how many entries in the sample were "good"?
#' 
#' @return numeric
#' @export
#' 
visqc_it_kendallt = function(data_matrix, 
                             exclude_na = TRUE, 
                             exclude_inf = TRUE, 
                             exclude_0 = TRUE, 
                             zero_value = 0, 
                             perspective = "local",
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
      cor_matrix[icol, jcol] = cor_matrix[jcol, icol] = kendallt(exclude_data[, icol], exclude_data[, jcol], perspective = perspective)
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
#' @param scale_max should everything be scaled compared to the maximum correlation?
#' @param diag_not_na should the diagonal entries reflect how many entries in the sample were "good"?
#' 
#' @return numeric
#' @export
#' 
visqc_it_kendallt_splitup = function(data_matrix, 
                             exclude_na = TRUE, 
                             exclude_inf = TRUE, 
                             exclude_0 = TRUE, 
                             zero_value = 0, 
                             perspective = "local",
                             scale_max = TRUE,
                             diag_good = TRUE){
  
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
  n_sample = ncol(exclude_data)
  # set everything to NA and let R take care of it
  
  ncore = future::nbrOfWorkers()
  names(ncore) = NULL
  
  n_each_col = seq(n_sample, 1, -1)
  n_todo = sum(n_each_col)
  n_each = ceiling(n_todo / ncore)
  
  split_indices = vector("list", ncore)
  start_loc = 1
  
  for (isplit in seq_along(split_indices)) {
    curr_cumsum = cumsum(n_each_col[seq(start_loc, n_sample)])
    is_over = sum(curr_cumsum > n_each)
    if (is_over >= 1) {
      stop_loc = min(which(curr_cumsum > n_each)) + start_loc
    } else {
      stop_loc = n_sample
    }
    split_indices[[isplit]] = seq(start_loc, stop_loc)
    if (stop_loc == n_sample) {
      return()
    }
    start_loc = stop_loc + 1
  }
  null_indices = purrr::map_lgl(split_indices, is.null)
  split_indices = split_indices[!null_indices]
  
  do_split = function(seq_range, exclude_data, perspective) {
    #seq_range = seq(in_range[1], in_range[2])
    #print(seq_range)
    tmp_cor = matrix(0, nrow = ncol(exclude_data), ncol = ncol(exclude_data))
    rownames(tmp_cor) = colnames(tmp_cor) = colnames(exclude_data)
    
    for (icol in seq_range) {
      for (jcol in seq(icol, ncol(exclude_data))) {
        #print(c(icol, jcol))
        tmp_cor[icol, jcol] = tmp_cor[jcol, icol] = kendallt(exclude_data[, icol], exclude_data[, jcol], perspective = perspective)
      }
    }
    tmp_cor
  }
  #tictoc::tic()
  split_cor = furrr::future_map(split_indices, do_split, exclude_data, perspective)
  #tictoc::toc()
  
  
  cor_matrix = matrix(0, nrow = ncol(exclude_data), ncol = ncol(exclude_data))
  rownames(cor_matrix) = colnames(cor_matrix) = colnames(exclude_data)
  for (isplit in split_cor) {
    cor_matrix = cor_matrix + isplit
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
