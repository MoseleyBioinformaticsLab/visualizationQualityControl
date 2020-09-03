#' matching features
#' 
#' For a given set of feature-sample matrices, calculates how many features are
#' in that sample, as well as in common to all other samples, and if provided,
#' in common within and outside the same sample group.
#' 
#' @param feature_matrix the feature to sample matrix.
#' @param zero_value what is the zero value? Default is NA
#' @param groups what are the groups
#' 
#' @export
#' @return data.frame
#' 
count_matching_features = function(feature_matrix, zero_value = NA, groups = NULL){
  ncol_feature = ncol(feature_matrix)
  
  if (!is.null(groups)) {
    nrow_group = nrow(groups)
    if (ncol_feature != nrow_group) {
      stop("Number of rows of groups does not match the number of samples!")
    }
  } else {
    groups = rep("1", ncol_feature)
  }
  
  if (is.na(zero_value)) {
    nonzero_features = !is.na(feature_matrix)
  } else {
    nonzero_features = feature_matrix != zero_value
  }
  
  
  matched = matrix(NA, nrow = ncol_feature, ncol = ncol_feature)
  
  for (icol in seq(1, ncol_feature)) {
    features1 = nonzero_features[, icol]
    for (jcol in (seq(icol, ncol_feature))) {
      features2 = nonzero_features[, jcol]
      matched[icol, jcol] = matched[jcol, icol] = sum(features1 & features2)
    }
  }
  colnames(matched) = colnames(feature_matrix)
  
  n_match = data.frame(sample_id = colnames(matched), self = diag(matched),
                       max_in_group = NA,
                       max_out_group = NA,
                       stringsAsFactors = FALSE)
  
  group_data = split_groups(feature_matrix, groups)
  
  for (igroup in names(group_data$groups)) {
    in_group = group_data$groups[[igroup]]
    other_groups = names(group_data$groups)[!names(group_data$groups) %in% igroup]
    other_samples = unlist(group_data$groups[other_groups], use.names = FALSE)
    
    for (isample in group_data$groups[[igroup]]) {
      other_in = in_group[!(in_group %in% isample)]
      n_match$max_in_group[isample] = max(matched[isample, other_in])
      n_match$max_out_group[isample] = max(matched[isample, other_samples])
    }
  }
  
  return(n_match)
  
}








