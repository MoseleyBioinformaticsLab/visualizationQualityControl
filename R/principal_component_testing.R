#' test scores and outcome
#' 
#' Given a matrix of PCA scores, set of sample attributes, and particular variables to test,
#' goes through and performs an ANOVA of the scores versus the attribute.
#' 
#' @param pca_scores matrix of scores from a PCA decomposition
#' @param sample_info data.frame of sample attributes
#' @param variables which columns to test against the scores
#'
#' @export
#' @return data.frame
visqc_test_pca_scores = function(pca_scores, sample_info, variables){
  pc_test = colnames(pca_scores)
  pc_2_variable = purrr::map_df(variables, function(in_var){
    pc_test = purrr::map_df(pc_test, function(in_pc){
      #message(paste0(in_var, " : ", in_pc))
      tmp_frame = data.frame(y = pca_scores[, in_pc],
                             x = sample_info[[in_var]])
      aov_res = stats::aov(y ~ x, data = tmp_frame)
      tidy_res = broom::tidy(aov_res)[1, ]
      tidy_res$PC = in_pc
      tidy_res$variable = in_var
      tidy_res
    })
    
  })
  pc_2_variable
}


#' test loadings
#' 
#' Given a matrix of loadings for principal components, and a set of components to test,
#' for each loading in each component, generates a null distribution from the other
#' loadings in all the other components, and reports a p-value for that loading.
#' 
#' @param loadings matrix of loadings from pca decomposition
#' @param test_columns names of the columns of the loadings to test
#' @param progress should progress be reported
#' @param direction should direction of loading be tested?
#' 
#' @export
#' @return named list
visqc_test_pca_loadings = function(loadings, test_columns, progress = FALSE, direction = FALSE){
  if (progress) {
    n_loading = nrow(loadings) * length(test_columns)
    progress_bar = knitrProgressBar::progress_estimated(n_loading)
  } else {
    progress_bar = NULL
  }
  loading_pvalues = purrr::map(test_columns, function(in_loading){
    all_cols = colnames(loadings)
    use_loadings = loadings[, in_loading]
    other_loadings = loadings[, !(all_cols %in% in_loading)]
    loading_results = test_loadings(use_loadings, other_loadings, progress_bar = progress_bar, direction = direction)
    loading_df = data.frame(p.value = loading_results, 
                            loading = use_loadings,
                            loading_index = in_loading,
                            stringsAsFactors = FALSE)
    loading_df
  })
  names(loading_pvalues) = test_columns
  loading_pvalues
}

test_loadings = function(use_loadings, other_loadings, progress_bar = NULL, direction = FALSE){
  use_index = seq(1, length(use_loadings))
  purrr::map_dbl(use_index, function(in_row){
    test_value = use_loadings[in_row]
    all_other = (as.vector(other_loadings[-in_row, ]))
    
    if (direction) {
      if (test_value < 0) {
        is_more = sum(all_other < test_value) / sum(all_other < 0)
      } else {
        is_more = sum(all_other > test_value) / sum(all_other >= 0)
      }
    } else {
      all_other = abs(all_other)
      test_value = abs(test_value)
      is_more = sum(all_other > test_value) / length(all_other)
    }
    
    knitrProgressBar::update_progress(progress_bar)
    is_more
  })
}
