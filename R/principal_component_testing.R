#' test scores and outcome
#' 
#' Given a matrix of PCA scores, set of sample attributes to test,
#' goes through and performs an ANOVA of the scores versus the attribute.
#' 
#' @param pca_scores matrix of scores from a PCA decomposition
#' @param sample_info data.frame of sample attributes to test
#'
#' @details
#' For each variable being tested, the function first checks the number of unique values
#' in that variable. If the number of unique values is 1, or equal to the number of rows,
#' then the variable is not tested and is excluded from the results.
#' 
#' 
#' @export
#' @return data.frame
visqc_test_pca_scores = function(pca_scores, sample_info){
  pc_test = colnames(pca_scores)
  pc_2_variable = purrr::imap_dfr(sample_info, function(var_col, var_name){
    pc_test = purrr::map_df(pc_test, function(in_pc){
      #message(paste0(var_name, " : ", in_pc))
      tmp_col = var_col
      tmp_col = tmp_col[!is.na(var_col) | !is.infinite(var_col) | !is.nan(var_col)]
      n_var = length(unique(tmp_col))
      is_1_all = FALSE
      if (n_var == 1) {
        is_1_all = TRUE
      }
      if ((n_var == nrow(sample_info)) && (is.character(var_col))) {
        is_1_all = TRUE
      }

      if (!is_1_all) {
        tmp_frame = data.frame(y = pca_scores[, in_pc],
          x = var_col)
        na_x = is.na(tmp_frame$x) | is.infinite(tmp_frame$x) | is.nan(tmp_frame$x)
        tmp_frame = tmp_frame[!na_x, ]
        aov_res = stats::aov(y ~ x, data = tmp_frame)
        tidy_res = broom::tidy(aov_res)[1, ]
        tidy_res$PC = in_pc
        tidy_res$variable = var_name
        return(tidy_res)
      } else {
        return(NULL)
      }
    })
    
  })
  pc_2_variable
}

#' correlate scores and outcome
#' 
#' Given a matrix of PCA scores, set of sample attributes to test,
#' goes through and performs an ICI-Kt of the scores versus the attribute.
#' 
#' 
#' @param pca_scores the scores matrix to test
#' @param sample_info data.frame of sample attributes to test
#' 
#' Important: All of the attributes must be numeric, or character.
#'   If character, they will be transformed to a factor, and the numeric
#'   factor levels will be used instead.
#'   If missing values are present, that is OK, as long as they are
#'   missing-not-at-random (i.e. missing at the low end of the values).
#'   
#' @export
#' @return data.frame
#' 
visqc_cor_pca_scores = function(pca_scores, sample_info){
  requireNamespace("ICIKendallTau", quietly = TRUE)
  
  pc_test = colnames(pca_scores)
  pc_2_variable = purrr::imap_dfr(sample_info, function(var_col, var_name){
    pc_test = purrr::map_df(pc_test, function(in_pc){
      #message(paste0(var_name, " : ", in_pc))
      if (is.character(var_col)) {
        var_col2 = factor(var_col)
        num_col = as.numeric(unclass(var_col2))
      } else {
        num_col = var_col
      }
      ici_res = ICIKendallTau::ici_kt(pca_scores[, in_pc], num_col, "global")
      tidy_res = data.frame(cor = ici_res["tau"], 
                            pvalue = ici_res["pvalue"],
                            PC = in_pc,
                            variable = var_name)
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
