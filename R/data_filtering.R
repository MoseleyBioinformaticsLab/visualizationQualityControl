#' keep features with percentage of non-zeros
#' 
#' Given a value matrix (features are columns, samples are rows), and sample classes, 
#' find those things that are not zero in at least a certain 
#' number of one of the classes, and keep them
#' 
#' @param data_matrix the matrix of values to work with 
#' @param sample_classes the classes of each sample
#' @param keep_num what number of samples in each class need a non-zero value (see Details)
#' 
#' @details This function is being deprecated and all code should use the
#'   \code{keep_non_zero_percentage} function instead.
#' 
#' @seealso keep_non_zero_percentage
#' @export
#' @return matrix
filter_non_zero_percentage <- function(data_matrix, sample_classes = NULL, keep_num = 0.75){
  .Deprecated("keep_non_zero_percentage", "visualizationQualityControl")
  keep_non_zero_percentage(data_matrix, sample_classes, keep_num)
}

#' keep features with percentage of non-zeros
#' 
#' Given a value matrix (features are columns, samples are rows), and sample classes, 
#' find those things that are \emph{not zero} in at least a certain 
#' number of samples in one of the classes, and keep those features for further
#' processing.
#' 
#' @param data_matrix the matrix of values to work with 
#' @param sample_classes the classes of each sample
#' @param keep_num what number of samples in each class need a non-zero value (see Details)
#' @param zero_value what number represents zero values
#' @param all is this an either / or OR does it need to be present in all?
#' 
#' @details The number of samples that must be non-zero can be expressed either as a whole
#'   number (that is greater than one), or as a fraction that will be be multiplied 
#'   by the number of samples in each class to get the lower limits for each of the classes.
#' 
#' @export
#' @return matrix
keep_non_zero_percentage <- function(data_matrix, sample_classes = NULL, keep_num = 0.75, zero_value = 0, all = FALSE){
  stopifnot(ncol(data_matrix) != 0)
  stopifnot(nrow(data_matrix) != 0)
  stopifnot(keep_num >= 0)
  stopifnot(!(keep_num >= nrow(data_matrix)))
  
  if (is.null(sample_classes)) {
    sample_classes <- rep("A", nrow(data_matrix))
  }
  uniq_classes <- unique(sample_classes)
  class_index <- lapply(uniq_classes, function(x){sample_classes %in% x})
  names(class_index) <- uniq_classes
  
  if (keep_num < 1) {
    min_notzero <- sapply(class_index, function(x){round(sum(x) * keep_num)})
  } else {
    min_notzero <- sapply(class_index, function(x){round(keep_num)})
  }
  
  if (all) {
    min_pass <- length(min_notzero)
  } else {
    min_pass <- 1
  }
  
  # what is this doing?
  # For each of the features (columns), check how many non-zero entries
  # there are for each class. If the number is greater than the limit
  # in at least one of the classes, keep that feature. This allows easy filtering
  # of those features that have more than the specified number of zeros in all
  # classes.
  # 
  # iterate over columns (features)
  has_min <- apply(data_matrix, 2, function(in_col){
    # how many are non-zero in each class
    n_pass <- sapply(class_index, function(index){sum(in_col[index] > zero_value)})
    # is the minimum reached in at least one class
    keep <- sum(n_pass >= min_notzero) >= min_pass
    keep
  })
  data_matrix[, has_min]
}
