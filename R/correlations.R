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
#' @details The data.frame has 6 columns:
#' \describe{
#'   \item{med_cor}{the median correlation with other samples}
#'   \item{sample_id}{the sample id, either the rowname or an index}
#'   \item{sample_class}{the class of the sample. If not provided, set to "C1"}
#'   \item{compare_class}{the class of the other sample}
#'   \item{plot_class}{\code{sample_class::compare_class} for easy grouping}
#'   \item{group}{whether the median is \emph{within} or \emph{between} groups of samples}
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
  
  sample_median_cor$plot_class <- paste0(sample_median_cor$sample_class, "::", sample_median_cor$compare_class)
  sample_median_cor <- sample_median_cor[(sample_median_cor$sample_class == sample_median_cor$compare_class), ]
  group = rep("within", nrow(sample_median_cor))
  group[sample_median_cor$sample_class != sample_median_cor$compare_class] = "between"
  
  sample_median_cor$group = group

  if (!between_classes) {
    sample_median_cor = sample_median_cor[sample_median_cor$group %in% "within", ]
  }
  rownames(sample_median_cor) <- NULL

  sample_median_cor
}

#' calculate median class correlations
#'
#' Given a correlation matrix the sample class information,
#' calculates the median correlations of the samples within
#' the class and between classes.
#'
#' @param cor_matrix the sample - sample correlations
#' @param sample_classes the sample classes as a character or factor
#'
#' @return matrix
#' @export
#'
#'
median_class_correlations <- function(cor_matrix, sample_classes = NULL){
  stopifnot(nrow(cor_matrix) == ncol(cor_matrix))
  n_sample <- nrow(cor_matrix)
  
  if (is.null(sample_classes)) {
    stop("You didn't provide sample classes to work with.")
  }
  
  if (is.null(rownames(cor_matrix))) {
    stopifnot(rownames(cor_matrix) == colnames(cor_matrix))
    sample_id <- paste0("S", seq(1, n_sample))
    rownames(cor_matrix) <- colnames(cor_matrix) <- sample_id
  } else {
    sample_id <- rownames(cor_matrix)
  }
  
  names(sample_classes) = rownames(cor_matrix)
  sample_classes = sort(sample_classes)
  
  all_comp = combn(names(sample_classes), 2)
  comp_df = data.frame(s1 = all_comp[1, ],
                       s2 = all_comp[2, ],
                       class1 = sample_classes[all_comp[1, ]],
                       class2 = sample_classes[all_comp[2, ]],
                       cor = NA)
  
  for (irow in seq_len(nrow(comp_df))) {
    comp_df[irow, "cor"] = cor_matrix[comp_df[irow, "s1"],
                                      comp_df[irow, "s2"]]
  }
  
  med_df = comp_df %>%
    dplyr::group_by(class1, class2) %>%
    dplyr::summarise(median = median(cor), .groups = "keep") %>%
    dplyr::ungroup() %>%
    data.frame()
  
  out_classes = unique(sample_classes)
  out_median = matrix(NA, nrow = length(out_classes),
                      ncol = length(out_classes))
  rownames(out_median) = colnames(out_median) = out_classes
  
  for (irow in seq_len(nrow(med_df))) {
    class1 = med_df[irow, "class1"]
    class2 = med_df[irow, "class2"]
    out_median[class1, class2] = med_df[irow, "median"]
  }
  
  out_median
}

# calculates if something is an outlier
.calc_outlier <- function(data, n_trim, n_sd, remove_missing){
  outlier_data <- apply(data, 1, function(x){
    is_bad <- is.infinite(x) | is.na(x) | is.nan(x)
    if (length(remove_missing) > 0) {
      all_bad <- is_bad | (x %in% remove_missing)
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
#' @param data the data matrix (samples are columns, rows are features)
#' @param sample_classes the sample classes
#' @param n_trim how many features to trim at each end (default is 3)
#' @param n_sd how many SD before treated as outlier (default is 5)
#' @param remove_missing what missing values be removed before calculating? (default is NA)
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
                             n_sd = 5, remove_missing = NA){
  n_sample <- ncol(data)

  if (is.null(sample_classes)) {
    use_classes <- factor(rep("C1", n_sample))
  } else if (is.factor(sample_classes)) {
    use_classes <- sample_classes
  } else {
    sample_classes <- factor(sample_classes)
    use_classes <- sample_classes
  }

  split_classes <- split(seq(1, n_sample), use_classes)

  if (is.null(colnames(data))) {
    sample_names <- paste0("S", seq(1, n_sample))
  } else {
    sample_names <- colnames(data)
  }

  frac_outlier_class <- lapply(names(split_classes), function(class_name){
    class_index <- split_classes[[class_name]]
    is_outlier <- .calc_outlier(data[, class_index, drop = FALSE], n_trim, n_sd, remove_missing)
    if (nrow(is_outlier) != nrow(data)) {
      is_outlier = t(is_outlier)
    }
    frac_outlier <- colSums(is_outlier) / nrow(data)
    data.frame(sample_id = sample_names[class_index], sample_class = class_name, frac = frac_outlier)
  })

  frac_outlier <- do.call(rbind, frac_outlier_class)
  frac_outlier
}

#' determine outliers
#' 
#' @param median_correlations median correlations
#' @param outlier_fraction outlier fractions
#' @param cor_weight how much weight for the correlation score?
#' @param frac_weight how much weight for the outlier fraction?
#' @param only_high should only things at the low end of score be removed?
#' 
#' @details For outlier sample detection, one should 
#'   first generate median correlations using
#'   `median_correlations`, and outlier fractions using
#'   `outlier_fraction`. If you only have one or the other,
#'   than you should use named arguments to only pass the one
#'   or the other. 
#'   
#'   Alternatively, you can change the weighting used for median correlations
#'   or outlier fraction, including setting them to 0.
#' 
#' @export
#' @return data.frame
determine_outliers = function(median_correlations = NULL, outlier_fraction = NULL,
                              cor_weight = 1, frac_weight = 1, only_high = TRUE){
  
  if (!is.null(median_correlations) && !is.null(outlier_fraction)) {
    full_data = dplyr::left_join(median_correlations, outlier_fraction, by = "sample_id", suffix = c(".cor", ".frac"))
    
    if ((nrow(full_data) != nrow(median_correlations)) || (nrow(full_data) != nrow(outlier_fraction))) {
      stop("Samples in median_correlations and outlier_fraction don't match up!")
    }
    
    if (!all(full_data$sample_class.cor == full_data$sample_class.frac)) {
      stop("The sample classes (sample_class) from median_correlations and outlier_fraction are NOT all the same!")
    }
    full_data$sample_class = full_data$sample_class.cor
    full_data$sample_class.cor = NULL
    full_data$sample_class.frac = NULL
  }
  
  if (is.null(outlier_fraction) && !is.null(median_correlations)) {
    full_data = median_correlations
    full_data$frac = 0
  }
  
  if (!is.null(outlier_fraction) && is.null(median_correlations)) {
    full_data = outlier_fraction
    full_data$med_cor = 0
  }
  
  data_score = (cor_weight * log(1 - full_data$med_cor)) + (frac_weight * log1p(full_data$frac))
  names(data_score) = full_data$sample_id
  split_score = split(data_score, full_data$sample_class)
  out_each = purrr::map(split_score, function(in_scores){
    score_out = boxplot.stats(in_scores)$out
    names(in_scores)[in_scores %in% score_out]
  })
  all_out = unlist(out_each)
  
  full_data$score = data_score
  full_data$outlier = FALSE
  full_data$outlier[full_data$sample_id %in% all_out] = TRUE
  
  if (only_high) {
    split_data = split(full_data, full_data$sample_class)
    full_data = purrr::map(split_data, \(in_data){
      mean_score = mean(in_data$score)
      wrong_side = in_data |>
        dplyr::filter(score < mean_score, outlier) |>
        dplyr::pull(sample_id)
      in_data$outlier[in_data$sample_id %in% wrong_side] = FALSE
      in_data
    }) |>
      dplyr::bind_rows()
  }
  
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
