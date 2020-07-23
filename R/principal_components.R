#' calculate pca contributions
#' 
#' Given a set of PCA scores, calculates their variance contributions, cumulative
#' contributions, and generates a percent label that can be used for labeling plots.
#' 
#' @param pca_scores matrix of scores, columns are each PC
#' 
#' @export
#' @return data.frame
visqc_score_contributions = function(pca_scores){
  pca_vars = apply(pca_scores, 2, var)
  pca_vars = data.frame(pc = names(pca_vars), variance = pca_vars,
                        stringsAsFactors = FALSE)
  pca_vars = dplyr::mutate(pca_vars, percent = variance / sum(variance),
                           cumulative = cumsum(percent))
  pca_labels = purrr::map2_chr(pca_vars$pc, pca_vars$percent, function(.x, .y){
    paste0(.x, " (", format(.y * 100, digits = 2, scientific = FALSE, justify = "none", trim = TRUE), "%)")
  })
  pca_vars$labels = pca_labels
  pca_vars
}
