#' Test for left censorship
#' 
#' Does a binomial test to check if the most likely cause of missing values
#' is due to values being below the limit of detection, or coming from a 
#' left-censored distribution.
#' 
#' @param in_data matrix or data.frame of numeric data
#' @param sample_classes which samples are in which class
#' @param global_na what represents zero or missing?
#' 
#' @export
#' @return 
test_left_censorship = function(in_data, sample_classes = NULL, global_na = c(0, NA))
{
  if (is.null(sample_classes)) {
    sample_classes = rep("A", nrow(in_data))
  }
  
  split_indices = split(seq(1, ncol(in_data)), sample_classes)
  in_data_missing = setup_missing_matrix(in_data, global_na)
  count_missing = rowSums(is.na(in_data_missing))
  keep_missing = count_missing > 0
  
  if (sum(keep_missing) == 0) {
    message("No features have missing entries, nothing to do!")
    return(NULL)
  }
  in_data_missing = in_data_missing[keep_missing, , drop = FALSE]
  
  out_up_down = purrr::map(seq_len(nrow(in_data_missing)), \(in_row){
    class_up_down = purrr::map(split_indices, \(in_split){
      tmp_split = in_data_missing[in_row, in_split]
      
      if (sum(is.na(tmp_split)) == 0) {
        return(NULL)
      }
      
      tmp_split = tmp_split[!is.na(tmp_split)]
      median_split = median(tmp_split)
      
      return(as.numeric(tmp_split < median_split))
    })
    return(unlist(class_up_down))
  })
  all_up_down = unlist(out_up_down, use.names = FALSE)
  binom.test(sum(all_up_down), length(all_up_down), p = 0.5, alternative = "two.sided")
}


setup_missing_matrix = function(data_matrix, global_na)
{
  exclude_loc = matrix(FALSE, nrow = nrow(data_matrix), ncol = ncol(data_matrix))
  if (length(global_na) > 0) {
    if (any(is.na(global_na))) {
      exclude_loc[is.na(data_matrix)] = TRUE
      global_na = global_na[!is.na(global_na)]
    }
    if (any(is.infinite(global_na))) {
      exclude_loc[is.infinite(data_matrix)] = TRUE
      global_na = global_na[!is.infinite(global_na)]
    }
  }
  if (length(global_na) > 0) {
    for (ival in global_na) {
      exclude_loc[data_matrix == ival] = TRUE
    }
  }
  out_data = data_matrix
  out_data[exclude_loc] = NA
  out_data
}

add_uniform_noise = function(n_rep, value, sd, use_zero = FALSE){
  n_value = length(value)
  
  n_sd = n_rep * n_value
  
  out_sd = rnorm(n_sd, 0, sd)
  out_sd = matrix(out_sd, nrow = n_value, ncol = n_rep)
  
  if (!use_zero){
    tmp_value = matrix(value, nrow = n_value, ncol = n_rep, byrow = FALSE)
    out_value = tmp_value + out_sd
  } else {
    out_value = out_sd
  }
  
  return(out_value)
}
