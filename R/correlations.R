#' information volume
#' 
#' calculates the information volume weight
#' 
#' @param in_x which things are in X
#' @param in_y which things are in Y
#' 
#' @details The information volume is based on how many things are actually
#'   not missing in both samples. For weighting of correlations, this value
#'   will be scaled by the maximum observed value across calculated correlations.
#' 
#' @export
#' @return numeric
information_volume <- function(in_x, in_y){
  in_both <- sum(in_x & in_y)
  in_both
}

both_samples <- function(in_x, in_y){
  in_both <- sum(in_x & in_y)
  in_both
}

either_samples <- function(in_x, in_y){
  in_either <- sum(in_x | in_y)
  in_either
}

neither_sample <- function(in_x, in_y){
  not_in_both <- sum(!in_x & !in_y)
  not_in_both
}

calculate_jaccard <- function(in_both, in_either){
  in_both / in_either
}

calculate_information_consistency <- function(in_both, not_in_both, total_features){
  (in_both + not_in_both) / total_features
}


#' globally information theoretic weighted pairwise correlation
#'
#' given a data matrix, calculate \emph{inter-row} correlation values that are weighted
#' by information content and information consistency.
#'
#' @param data_matrix the data
#' @param use what data to use, see \code{cor}
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#' @param exclude_0 should 0 values be excluded (default FALSE)
#' @param zero_value what value represents zero (default is 0)
#' @param method which method of correlation to use
#' @param scale_information_content should the information content be scaled by the maximum value?
#'
#' @details The function returns a named list with:
#'   \describe{
#'     \item{cor}{the weighted correlation matrix}
#'     \item{keep}{a logical matrix indicating which points passed filtering}
#'     \item{raw}{the raw correlations}
#'     \item{content}{the weighting calculated from information content (weighted jaccard)}
#'     \item{consistency}{weighted information consistency}
#'     }
#'
#'
#' @return list
#' @export
globally_it_weighted_pairwise_correlation <- function(data_matrix, use = "pairwise.complete.obs", 
                                        exclude_na = TRUE, exclude_inf = TRUE, 
                                        exclude_0 = FALSE, zero_value = 0, method = "pearson",
                                        scale_information_content = TRUE){
  correlation_values <- pairwise_correlation(data_matrix, use = use, exclude_na = exclude_na,
                                             exclude_inf = exclude_inf, exclude_0 = exclude_0,
                                             zero_value = zero_value, method = method)
  cor_values <- correlation_values$cor
  keep_matrix <- t(correlation_values$keep)
  
  in_both <- in_either <- not_in_both <- matrix(0, nrow = ncol(keep_matrix), ncol = ncol(keep_matrix))
  
  for (icol in seq_len(ncol(keep_matrix))) {
    for (jcol in seq_len(ncol(keep_matrix))) {
      #print(c(icol, jcol))
      in_both[icol, jcol] <- in_both[jcol, icol] <- both_samples(keep_matrix[, icol], keep_matrix[, jcol])
      in_either[icol, jcol] <- in_either[jcol, icol] <- either_samples(keep_matrix[, icol], keep_matrix[, jcol])
      not_in_both[icol, jcol] <- not_in_both[jcol, icol] <- neither_sample(keep_matrix[, icol], keep_matrix[, jcol])
    }
  }
  
  total_features <- ncol(data_matrix)
  
  jaccard_info <- calculate_jaccard(in_both, in_either)
  if (scale_information_content) {
    off_diag_values <- jaccard_info[lower.tri(jaccard_info)]
    max_jaccard <- max(off_diag_values)
    information_content <- jaccard_info / max_jaccard
    diag(information_content) <- 1
  } else {
    information_content <- jaccard_info
  }
  
  information_consistency <- calculate_information_consistency(in_both, not_in_both, total_features)
  diag(information_consistency) <- diag(in_both) / total_features
  
  global_correlation <- cor_values * information_content * information_consistency
  
  return(list(cor = global_correlation, raw = cor_values,
       keep = t(keep_matrix), content = information_content,
       consistency = information_consistency))
}


#' locally information theoretic weighted pairwise correlation
#'
#' given a data matrix, calculate \emph{inter-row} correlation values that are weighted
#' by information content where only the pairwise information (jaccard index) is considered.
#'
#' @param data_matrix the data
#' @param use what data to use, see \code{cor}
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#' @param exclude_0 should 0 values be excluded (default FALSE)
#' @param zero_value what value represents zero (default is 0)
#' @param method which method of correlation to use
#'
#' @details The function returns a named list with:
#'   \describe{
#'     \item{cor}{the weighted correlation matrix}
#'     \item{keep}{a logical matrix indicating which points passed filtering}
#'     \item{raw}{the raw correlations}
#'     \item{jaccard}{the jaccard information index used to weight the values}
#'     }
#'
#'
#' @return list
#' @export
locally_it_weighted_pairwise_correlation <- function(data_matrix, use = "pairwise.complete.obs", 
                                        exclude_na = TRUE, exclude_inf = TRUE, 
                                        exclude_0 = FALSE, zero_value = 0, method = "pearson"){
  correlation_values <- pairwise_correlation(data_matrix, use = use, exclude_na = exclude_na,
                                             exclude_inf = exclude_inf, exclude_0 = exclude_0,
                                             zero_value = zero_value, method = method)
  cor_values <- correlation_values$cor
  keep_matrix <- t(correlation_values$keep)
  
  in_both <- in_either <- matrix(0, nrow = ncol(keep_matrix), ncol = ncol(keep_matrix))
  
  for (icol in seq_len(ncol(keep_matrix))) {
    for (jcol in seq_len(ncol(keep_matrix))) {
      #print(c(icol, jcol))
      in_both[icol, jcol] <- in_both[jcol, icol] <- both_samples(keep_matrix[, icol], keep_matrix[, jcol])
      in_either[icol, jcol] <- in_either[jcol, icol] <- either_samples(keep_matrix[, icol], keep_matrix[, jcol])
    }
  }
  
  jaccard_info <- calculate_jaccard(in_both, in_either)
  
  local_correlation <- cor_values * jaccard_info
  
  return(list(cor = local_correlation, raw = cor_values,
              keep = t(keep_matrix), jaccard = jaccard_info))
}

#' correspondence
#' 
#' @param in_x which things are in X
#' @param in_y which things are in Y
#' @param not_both also include those things NOT in both
#' 
#' @details Correspondence is based on how many things agree compared to how
#'   many things do not agree. So by default it is the ratio of
#'   
#'     things observed in both / things observed in either
#'     
#'   This should handle well most types of data, including sparse data sets.
#'   However, there are cases where the fact that values are missing in both
#'   things being compared counts as information. This alternative method is
#'   triggered by setting \code{not_both = TRUE}, resulting in the calculation:
#'   
#'     (things observed in both + things not observed in both) / total # of things
#'     
#'   This value will be scaled to the largest observed correspondence for weighting
#'   as well.
#' 
#' @export
#' @return numeric
correspondence <- function(in_x, in_y, not_both = FALSE){
  in_both <- sum(in_x & in_y)
  
  if (not_both) {
    not_in_both <- sum(!in_x & !in_y)
    corresp_ratio <- (in_both + not_in_both) / length(in_x)
  } else {
    in_either <- sum(in_x | in_y)
    corresp_ratio <- in_both / in_either
  }
  
  corresp_ratio
}


#' calculate weights
#' 
#' calculates the weights for correlations
#' 
#' @param nonzero_loc the non-zero location logical matrix
#' @param not_both consider things FALSE in both as information
#' 
#' @export
#' @return list
calculate_weights <- function(nonzero_loc, not_both = FALSE){
  info_matrix <- cor_matrix <- matrix(0, ncol(nonzero_loc), ncol(nonzero_loc))
  for (icol in seq_len(ncol(nonzero_loc) - 1)) {
    for (jcol in seq(icol + 1, ncol(nonzero_loc))) {
      #print(c(icol, jcol))
      info_matrix[icol, jcol] <- info_matrix[jcol, icol] <- information_volume(nonzero_loc[, icol], nonzero_loc[, jcol])
      cor_matrix[icol, jcol] <- cor_matrix[jcol, icol] <- correspondence(nonzero_loc[, icol], nonzero_loc[, jcol], not_both)
      
    }
  }
  info_weight <- info_matrix / max(info_matrix)
  
  if (not_both) {
    cor_weight <- cor_matrix / max(cor_matrix)
  } else {
    cor_weight <- cor_matrix
  }

  diag(info_weight) <- diag(cor_weight) <- 1
  
  n_keep <- colSums(nonzero_loc)
  keep_frac <- n_keep / max(n_keep)
  
  diag(info_weight) <- keep_frac
  
  return(list(info = info_weight, correspondence = cor_weight))
}

#' pairwise correlation
#'
#' given a data matrix, calculate \emph{inter-row} correlation values.
#'
#' @param data_matrix the data
#' @param use what data to use, see \code{cor}
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#' @param exclude_0 should 0 values be excluded (default FALSE)
#' @param zero_value what value represents zero (default is 0)
#' @param method which method of correlation to use
#' @param weight should the correlations be weighted by information content and correspondence?
#' @param not_both should correspondence include things FALSE in both?
#'
#' @details The function returns a named list with:
#'   \describe{
#'     \item{cor}{the correlation matrix}
#'     \item{keep}{a logical matrix indicating which points passed filtering}
#'     }
#'
#'
#' @return list
#' @export
pairwise_correlation <- function(data_matrix, use = "pairwise.complete.obs", 
                                 exclude_na = TRUE, exclude_inf = TRUE, 
                                 exclude_0 = FALSE, zero_value = 0, method = "pearson"){

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

  # set everything to NA and let R take care of it
  data_matrix[exclude_loc] <- NA

  calc_cor <- cor(data_matrix, use = use, method = method)
  
  return(list(cor = calc_cor, keep = t(!exclude_loc)))
}

#' pairwise correlation counts
#'
#' given the \code{keep} entry from \code{pairwise_correlation}, get the number
#' of things that would have been used for each of the pairwise comparisons. Note
#' that pairwise comparisons of the rows are used.
#'
#' @param keep_matrix \code{logical} matrix
#'
#' @export
#' @return matrix
pairwise_correlation_counts <- function(keep_matrix){
  stopifnot(is.logical(keep_matrix))

  n_entry <- nrow(keep_matrix)

  out_count <- matrix(0, n_entry, n_entry)

  for (i_entry in seq(1, (n_entry - 1))) {
    for (j_entry in seq(i_entry, n_entry)) {
      out_count[i_entry, j_entry] <- out_count[j_entry, i_entry] <-
        sum(keep_matrix[i_entry, ] & keep_matrix[j_entry, ])
    }
  }
  out_count
}

#' pairwise correlation keep both zero
#'
#' given a data matrix, calculate inter-row correlation values, but if both
#' entries have zero, keep it
#'
#' @param data_matrix the data
#' @param use what data to use, "pairwise" or "complete"
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#' @param exclude_0 should 0 values be excluded (default FALSE)
#' @param method which method of correlation to use
#'
#' @return matrix
#' @export
pairwise_correlation_both0 <- function(data_matrix, use = "pairwise", exclude_na = TRUE, exclude_inf = TRUE, exclude_0 = FALSE, method = "pearson"){
  n_entry <- nrow(data_matrix)
  out_cor <- matrix(0, nrow = nrow(data_matrix), ncol = nrow(data_matrix))
  rownames(out_cor) <- colnames(out_cor) <- rownames(data_matrix)
  diag(out_cor) <- 1

  na_loc <- matrix(FALSE, nrow = n_entry, ncol = ncol(data_matrix))
  inf_loc <- na_loc
  zero_loc <- na_loc

  if (exclude_na){
    na_loc <- is.na(data_matrix)
  }

  if (exclude_inf){
    inf_loc <- is.infinite(data_matrix)
  }

  if (exclude_0){
    zero_loc <- data_matrix == 0
  }

  exclude_loc <- na_loc | zero_loc | inf_loc

  if (use == "complete"){
    keep_vals <- apply(exclude_loc, 2, function(x){sum(!x) == n_entry})
    keep_vals <- matrix(keep_vals, nrow = n_entry, ncol = ncol(data_matrix), byrow = FALSE)
  } else {
    keep_vals <- !exclude_loc
  }

  for (i in seq(1, n_entry)){
    for (j in seq(i, n_entry)){
      both_zero <- zero_loc[i, ] & zero_loc[j, ]
      use_vals <- (keep_vals[i, ] & keep_vals[j, ]) | both_zero

      if (sum(use_vals) > 1){
        #print(c(i, j))
        out_cor[i, j] <- out_cor[j, i] <- cor(data_matrix[i, use_vals], data_matrix[j, use_vals], method = method)
      }
    }
  }
  out_cor
}


#' pairwise correlation multicore
#'
#' given a data matrix, calculate inter-row correlation values
#'
#' @param data_matrix the data
#' @param use what data to use, "pairwise" or "complete"
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#' @param exclude_0 should 0 values be excluded (default FALSE)
#' @param method which method of correlation to use
#'
#' @return matrix
#' @export
pairwise_correlation_multicore <- function(data_matrix, use = "pairwise", exclude_na = TRUE, exclude_inf = TRUE, exclude_0 = FALSE, method = "pearson"){
  n_entry <- nrow(data_matrix)
  out_cor <- matrix(0, nrow = nrow(data_matrix), ncol = nrow(data_matrix))
  rownames(out_cor) <- colnames(out_cor) <- rownames(data_matrix)
  diag(out_cor) <- 1

  na_loc <- matrix(FALSE, nrow = n_entry, ncol = ncol(data_matrix))
  inf_loc <- na_loc
  zero_loc <- na_loc

  if (exclude_na){
    na_loc <- is.na(data_matrix)
  }

  if (exclude_inf){
    inf_loc <- is.infinite(data_matrix)
  }

  if (exclude_0){
    zero_loc <- data_matrix == 0
  }

  exclude_loc <- na_loc | zero_loc | inf_loc

  if (use == "complete"){
    keep_vals <- apply(exclude_loc, 2, function(x){sum(!x) == n_entry})
    keep_vals <- matrix(keep_vals, nrow = n_entry, ncol = ncol(data_matrix), byrow = FALSE)
  } else {
    keep_vals <- !exclude_loc
  }

  all_comparisons <- expand.grid(seq(1, n_entry), seq(1, n_entry))
  all_comparisons <- all_comparisons[(all_comparisons[,2] > all_comparisons[,1]), ]
  all_comparisons <- as.matrix(all_comparisons)

  correlation_function <- function(x){
    n_vals <- length(x)
    out_value <- rep(0, n_vals)
    do_comparison <- all_comparisons[x, ]

    for (i_val in seq_len(n_vals)){
      i <- do_comparison[1]
      j <- do_comparison[2]

      use_vals <- keep_vals[i, ] & keep_vals[j, ]

      if (sum(use_vals) > 1){
        #print(c(i, j))
        out_value[i_val] <- cor(data_matrix[i, use_vals], data_matrix[j, use_vals], method = method)
      }
    }

    out_value
  }

  all_cors <- parallel::pvec(seq(1, nrow(all_comparisons)), correlation_function)

  out_cor[all_comparisons[,1], all_comparisons[,2]] <- all_cors
  out_cor[all_comparisons[,2], all_comparisons[,1]] <- all_cors
  out_cor
}

#' pairwise distance
#'
#' given a data matrix, calculate inter-row distances
#'
#' @param data_matrix the data
#' @param use what data to use, "pairwise" or "complete"
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#' @param exclude_0 should 0 values be excluded (default FALSE)
#'
#' @return matrix
#' @export
pairwise_distance <- function(data_matrix, use = "pairwise", exclude_na = TRUE, exclude_inf = TRUE, exclude_0 = FALSE){
  n_entry <- nrow(data_matrix)
  out_cor <- matrix(0, nrow = nrow(data_matrix), ncol = nrow(data_matrix))
  rownames(out_cor) <- colnames(out_cor) <- rownames(data_matrix)
  diag(out_cor) <- 1

  na_loc <- matrix(FALSE, nrow = n_entry, ncol = ncol(data_matrix))
  inf_loc <- na_loc
  zero_loc <- na_loc

  if (exclude_na){
    na_loc <- is.na(data_matrix)
  }

  if (exclude_inf){
    inf_loc <- is.infinite(data_matrix)
  }

  if (exclude_0){
    zero_loc <- data_matrix == 0
  }

  exclude_loc <- na_loc | zero_loc | inf_loc

  if (use == "complete"){
    keep_vals <- apply(exclude_loc, 2, function(x){sum(!x) == n_entry})
    keep_vals <- matrix(keep_vals, nrow = n_entry, ncol = ncol(data_matrix), byrow = FALSE)
  } else {
    keep_vals <- !exclude_loc
  }

  for (i in seq(1, n_entry)){
    for (j in seq(i, n_entry)){
      use_vals <- keep_vals[i, ] & keep_vals[j, ]

      if (sum(use_vals) > 1){
        #print(c(i, j))
        out_cor[j, i] <- out_cor[i, j] <- sqrt(sum((data_matrix[i, use_vals] - data_matrix[j, use_vals]) ^ 2))
      }
    }
  }
  out_cor
}

#' pairwise non-zero
#'
#' given a data matrix, how many entries are non-zero on a pairwise basis
#'
#' @param data_matrix the data
#' @param use how to work with the data
#' @param exclude_na should NA values be excluded (default TRUE)
#' @param exclude_inf should Inf values be excluded (default TRUE)
#'
#' @return matrix
#' @export
pairwise_nonzero <- function(data_matrix, use = "pairwise", exclude_na = TRUE, exclude_inf = TRUE){
  n_entry <- nrow(data_matrix)
  out_cor <- matrix(0, nrow = nrow(data_matrix), ncol = nrow(data_matrix))
  rownames(out_cor) <- colnames(out_cor) <- rownames(data_matrix)
  diag(out_cor) <- 1

  na_loc <- matrix(FALSE, nrow = n_entry, ncol = ncol(data_matrix))
  inf_loc <- na_loc

  if (exclude_na){
    na_loc <- is.na(data_matrix)
  }

  if (exclude_inf){
    inf_loc <- is.infinite(data_matrix)
  }

  zero_loc <- data_matrix == 0

  exclude_loc <- na_loc | inf_loc

  keep_vals <- !exclude_loc

  for (i in seq(1, n_entry)){
    for (j in seq(i, n_entry)){
      use_vals <- keep_vals[i, ] & keep_vals[j, ]

      zero_both <- zero_loc[i, ] & zero_loc[j, ]

      sum_non_zero <- sum(use_vals) - sum(zero_both)
      out_cor[j, i] <- out_cor[i, j] <- sum_non_zero
    }
  }
  out_cor
}

#' calculate median correlations
#'
#' Given a correlation matrix and optionally the sample class information,
#' calculates the median correlations of each sample to all other samples in
#' the same class. May be useful for determining outliers.
#'
#' @param cor_matrix the sample - sample correlations
#' @param sample_classes the sample classes as a character or factor
#' @param between_classes should the between class correlations be evaluated?
#'
#' @return data.frame
#' @export
#'
#' @details The data.frame may have 5 columns, first three are always present,
#'   the second two come up if \code{between_classes = TRUE}:
#' \describe{
#'   \item{med_cor}{the median correlation with other samples}
#'   \item{sample_id}{the sample id, either the rowname or an index}
#'   \item{sample_class}{the class of the sample. If not provided, set to "C1"}
#'   \item{compare_class}{the class of the other sample}
#'   \item{plot_class}{\code{sample_class::compare_class} for easy grouping}
#' }
#'
median_correlations <- function(cor_matrix, sample_classes = NULL, between_classes = FALSE){
  stopifnot(nrow(cor_matrix) == ncol(cor_matrix))
  n_sample <- nrow(cor_matrix)

  if (is.null(sample_classes)) {
    use_classes <- factor(rep("C1", n_sample))
  } else if (is.factor(sample_classes)) {
    use_classes <- sample_classes
  } else {
    sample_classes <- factor(sample_classes)
    use_classes <- sample_classes
  }

  if (is.null(rownames(cor_matrix))) {
    stopifnot(rownames(cor_matrix) == colnames(cor_matrix))
    sample_id <- paste0("S", seq(1, n_sample))
    rownames(cor_matrix) <- colnames(cor_matrix) <- sample_id
  } else {
    sample_id <- rownames(cor_matrix)
  }

  split_classes <- split(sample_id, use_classes)
  names(use_classes) <- sample_id

  sample_median_cor <- lapply(sample_id, function(in_sample){
    #message(in_sample)
    sample_loc <- which(colnames(cor_matrix) %in% in_sample)
    
    sample_cor <- cor_matrix[sample_loc, -sample_loc]
    
    cor_by_class <- lapply(split_classes, function(class_ids){
      #message(class_ids[1])
      class_samples <- setdiff(class_ids, in_sample)
      if (length(class_samples) > 0) {
        med_cor = median(sample_cor[class_samples])
      } else {
        med_cor = NA
      }
        data.frame(sample_id = in_sample,
                   med_cor = med_cor,
                   sample_class = as.character(use_classes[in_sample]),
                   compare_class = as.character(unique(use_classes[class_ids])),
                   stringsAsFactors = FALSE)
      })
      
    cor_by_class <- do.call(rbind, cor_by_class)
  })
  sample_median_cor <- do.call(rbind, sample_median_cor)
  
  if (between_classes) {
    sample_median_cor$plot_class <- paste0(sample_median_cor$sample_class, "::", sample_median_cor$compare_class)
  } else {
    sample_median_cor <- sample_median_cor[(sample_median_cor$sample_class == sample_median_cor$compare_class), ]
    sample_median_cor$compare_class <- NULL
  }

  rownames(sample_median_cor) <- NULL

  sample_median_cor
}


# calculates if something is an outlier
.calc_outlier <- function(data, n_trim, n_sd, remove_0){
  outlier_data <- apply(data, 2, function(x){
    is_bad <- is.infinite(x) | is.na(x) | is.nan(x)
    if (remove_0) {
      all_bad <- is_bad | (x == 0)
    } else {
      all_bad <- is_bad
    }
    y <- x[!all_bad]
    n_y <- length(y)
    y <- sort(y)
    y_start <- n_trim + 1
    y_end <- n_y - (n_trim + 1)
    if ((y_end <= y_start) || (y_end <= 0)) {
      y_start <- 1
      y_end <- n_y
    }
    y <- y[y_start:y_end]
    n_y <- length(y)
    if (n_y >= 3) {
      y_mean <- mean(y)
      y_sd <- sd(y)
      y_lo <- y_mean - (n_sd * y_sd)
      y_hi <- y_mean + (n_sd * y_sd)

      x_out <- !((x >= y_lo) & (x <= y_hi))

      # set those things that were outliers & bad in beginning to FALSE so
      # we don't overestimate the outlier proportion
      x_out[!(x_out & !all_bad)] <- FALSE
    } else {
      x_out <- rep(FALSE, length(x))
    }

    x_out
  })
  return(outlier_data)
}


#' fraction of outliers
#'
#' Calculates the fraction of entries in each sample that are more than \code{X}
#' standard deviations from the trimmed mean. See Details.
#'
#' @param data the data matrix (samples are rows, columns are features)
#' @param sample_classes the sample classes
#' @param n_trim how many features to trim at each end (default is 3)
#' @param n_sd how many SD before treated as outlier (default is 5)
#' @param remove_0 should zeros be removed before calculating? (default is FALSE)
#'
#' @details Based on the Gerlinski paper \href{https://dx.doi.org/10.1093/bioinformatics/btv425}{link}
#' for each feature (in a sample class), take the range across all the samples,
#' remove the \code{n_trim} lowest and highest values, and calculate the \code{mean}
#' and \code{sd}, and the actual upper and lower ranges of \code{n_sd} from the
#' \code{mean}. For each sample and feature, determine if \emph{within} or \emph{outside}
#' that limit. Fraction is reported as the number of features outside the range.
#'
#' Returns a \code{data.frame} with:
#' \describe{
#'   \item{sample_id}{the sample id, \code{rownames} are used if available, otherwise
#'   this is an index}
#'   \item{sample_class}{the class of the sample if \code{sample_classes} were provided,
#'   otherwise given a default of "C1"}
#'   \item{frac}{the actual outlier fraction calculated for that sample}
#' }
#'
#' @export
#' @return data.frame
outlier_fraction <- function(data, sample_classes = NULL, n_trim = 3,
                             n_sd = 5, remove_0 = FALSE){
  n_sample <- nrow(data)

  if (is.null(sample_classes)) {
    use_classes <- factor(rep("C1", n_sample))
  } else if (is.factor(sample_classes)) {
    use_classes <- sample_classes
  } else {
    sample_classes <- factor(sample_classes)
    use_classes <- sample_classes
  }

  split_classes <- split(seq(1, n_sample), use_classes)

  if (is.null(rownames(data))) {
    sample_names <- paste0("S", seq(1, n_sample))
  } else {
    sample_names <- rownames(data)
  }

  frac_outlier_class <- lapply(names(split_classes), function(class_name){
    class_index <- split_classes[[class_name]]
    is_outlier <- .calc_outlier(data[class_index, , drop = FALSE], n_trim, n_sd, remove_0)
    frac_outlier <- rowSums(is_outlier) / ncol(data)
    data.frame(sample_id = sample_names[class_index], sample_class = class_name, frac = frac_outlier)
  })

  frac_outlier <- do.call(rbind, frac_outlier_class)
  frac_outlier
}

#' determine outliers
#' 
#' @param median_correlations median correlations
#' @param outlier_fraction outlier fractions
#' 
#' @export
#' @return data.frame
determine_outliers = function(median_correlations, outlier_fraction,
                              cor_weight = 1, frac_weight = 1){
  full_data = dplyr::left_join(median_correlations, outlier_fraction, by = "sample_id", suffix = c(".cor", ".frac"))
  
  if ((nrow(full_data) != nrow(median_correlations)) || (nrow(full_data) != nrow(outlier_fraction))) {
    stop("Samples in median_correlations and outlier_fraction don't match up!")
  }
  
  if (!all(full_data$sample_class.cor == full_data$sample_class.frac)) {
    stop("The sample classes (sample_class) from median_correlations and outlier_fraction are NOT all the same!")
  }
  
  data_score = (cor_weight * log(1 - full_data$med_cor)) + (frac_weight * log1p(full_data$frac))
  names(data_score) = full_data$sample_id
  split_score = split(data_score, full_data$sample_class.frac)
  out_each = purrr::map(split_score, function(in_scores){
    score_out = boxplot.stats(in_scores)$out
    names(in_scores)[in_scores %in% score_out]
  })
  all_out = unlist(out_each)
  
  full_data$score = data_score
  full_data$outlier = FALSE
  full_data$outlier[full_data$sample_id %in% all_out] = TRUE
  
  full_data
  
}

#' grp_cor_data
#'
#' Example data used for demonstrating \code{median correlation}. A list with
#' 2 named entries:
#'
#' \describe{
#'   \item{data}{a data matrix with 100 rows and 20 columns}
#'   \item{class}{a character vector of 20 entries denoting classes}
#' }
#'
#' The data comes from two groups of samples, where there is ~0.85 correlation
#' within each group, and ~0.53 correlation between groups.
#'
#' @format List with 2 entries, data and class
#' @source Robert M Flight
"grp_cor_data"

#' grp_exp_data
#'
#' Example data that requires log-transformation before doing PCA or other QC.
#' A list with 2 named entries:
#'
#' \describe{
#'   \item{data}{a data matrix with 1000 rows and 20 columns}
#'   \item{class}{a character vector of 20 entries denoting classes}
#' }
#'
#' The data comes from two groups of samples, where there is ~0.80 correlation
#' within each group, and ~0.38 correlation between groups.
#'
#' @format List with 2 entries, data and class
#' @source Robert M Flight
"grp_exp_data"
