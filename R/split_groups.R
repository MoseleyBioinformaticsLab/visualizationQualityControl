#' split groups
#' 
#' Given a matrix and a data.frame, character vector or data.frame of groups,
#' splits the indices / names of the matrix into groups appropriately. This function
#' assumes that the matrix and the groups are in the correct order!!
#' 
#' @param in_matrix the matrix we want to split up
#' @param groups a data.frame, character vector or factor
#' 
#' @export
#' @return list of groups
split_groups = function(in_matrix, groups = NULL){
  ncol_matrix = ncol(in_matrix)
  nrow_matrix = nrow(in_matrix)
  
  if (is.null(groups)) {
    min_matrix = min(c(ncol_matrix, nrow_matrix))
    groups = rep("G", min_matrix)
  }
  
  if (is.data.frame(groups)) {
    n_group = nrow(groups)
  } else {
    n_group = length(groups)
  }
  
  # check if columns or rows matches
  if ((nrow_matrix == n_group) && (ncol_matrix != n_group)) {
    in_matrix = t(in_matrix)
    ncol_matrix = ncol(in_matrix)
  } else if (ncol_matrix == n_group) {
    in_matrix = in_matrix
  } else {
    stop("Rows or columns of the matrix match number of entries in groups! Please check dimensions of both!")
  }
  
  matrix_names <- colnames(in_matrix)
  
  # we use ALL of the columns to generate factors that are used for splitting!
  if (is.data.frame(groups)) {
    if (ncol(groups) > 1) {
      use_groups <- do.call(paste, c(groups, sep="."))
    } else {
      use_groups <- groups[, 1]
    }
    
  } else {
    use_groups = groups
  }
  
  
  split_indices <- split(seq(1, ncol_matrix), use_groups)
  n_indices <- lapply(split_indices, function(x){
    length(x)
  })
  keep_indices <- n_indices > 1
  
  if (sum(!keep_indices) > 0) {
    warning(paste0("Removing groups: ", paste(names(split_indices)[!keep_indices], collapse = ", ")))
  }
  return(list(matrix = in_matrix, groups = split_indices))
}
