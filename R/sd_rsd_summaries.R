#' summarize data
#'
#' summarizes a matrix or data.frame, where columns are samples
#' and rows are features
#'
#' @param in_data matrix or data.frame
#' @param sample_classes which samples are in which class
#' @param avg_function which function to use for summary
#' @param log_transform apply a log-transform to the mean
#' @param remove_missing remove missing values before summarizing
#'
#' @return data.frame
#' @export
summarize_data <- function(in_data, sample_classes=NULL, avg_function = mean,
                           log_transform = FALSE, remove_missing = NA){
  if (is.null(sample_classes)){
    sample_classes <- rep("A", nrow(in_data))
  }

  n_feature <- nrow(in_data)
  split_indices <- split(seq(1, ncol(in_data)), sample_classes)

  split_values <- lapply(split_indices, function(use_index){
    tmp_mean_sd <- apply(in_data[use_index, , drop = FALSE], 1, function(x){
      if (length(remove_missing) > 0) {
        x <- x[!(x %in% remove_missing)]
      }
      c(avg_function(x), sd(x), max(x) - min(x))
    })
    tmp_mean_sd <- t(tmp_mean_sd)
    colnames(tmp_mean_sd) <- c("mean", "sd", "diff")
    data.frame(mean = tmp_mean_sd[, "mean"],
               var = c(tmp_mean_sd[, "sd"],
                       tmp_mean_sd[, "sd"] / tmp_mean_sd[, "mean"],
                       tmp_mean_sd[, "diff"]),
               type = rep(c("sd", "rsd", "diff"), each = n_feature))
  })

  out_data <- do.call(rbind, split_values)
  class_rep <- vapply(split_values, nrow, numeric(1))
  out_data$class <- rep(names(split_values), times = class_rep)
  out_data$type <- factor(out_data$type, ordered = TRUE, levels = c("diff", "sd", "rsd"))

  if (is.function(log_transform)) {
    out_data$log_mean <- log_transform(out_data$mean)
  }
  out_data
}




#' calculate F-ratio
#'
#' given a data matrix of samples (columns) and features (rows), and a vector of classes (character or factor),
#' calculate an F-ratio for each feature.
#'
#' @param data the data matrix, with samples (columns) and features (rows)
#' @param data_classes what are the classes of the samples (columns)
#'
#' @return vector
#' @export
calculate_fratio <- function(data, data_classes){
  if(is.character(data_classes)){
    data_classes <- factor(data_classes)
  }

  all_means <- rowMeans(data)

  split_indices <- split(seq(1, ncol(data)), data_classes)
  n_sample <- ncol(data)
  n_group <- length(split_indices)
  split_data <- lapply(split_indices, function(in_index){data[, in_index, drop = FALSE]})

  group_means <- lapply(split_data, rowMeans)
  group_var <- lapply(split_data, function(in_data){apply(in_data, 1, var)})
  group_count <- lapply(split_data, ncol)

  weight_var <- function(count, var, sub1 = TRUE){
    if (sub1){
      count <- count - 1
    }
    count * var
  }
  within_var <- Map(weight_var, group_count, group_var)
  within_var <- do.call(rbind, within_var)
  within_var <- colSums(within_var) / (n_sample - n_group)

  between_var <- lapply(group_means, function(in_mean){(in_mean - all_means)^2})
  between_var <- Map(weight_var, group_count, between_var, FALSE)
  between_var <- do.call(rbind, between_var)
  between_var <- colSums(between_var) / (n_group - 1)

  f_ratio <- between_var / within_var
  f_ratio
}

#' calculate values from summaries
#'
#' given a data.frame of means and variances, calculate mean sd at low end and
#' mean rsd at high end.
#'
#' @param data data.frame of means and variances
#' @param low_cut means <= this value used for average sd
#' @param hi_cut means >= this value used for average rsd
#'
#' @import dplyr
#' @export
#' @return vector
calc_sd_rsd <- function(data, low_cut, hi_cut = NULL){
  if (is.null(hi_cut)){
    hi_cut <- min(data[, "mean"]) > low_cut
  }

  sd_mn <- filter(data, mean <= low_cut, type == "sd") %>% summarise(., mean(var))
  rsd_mn <- filter(data, mean >= hi_cut, type == "rsd") %>% summarise(., mean(var))

  return(c(sd = sd_mn, rsd = rsd_mn))
}

#' calculate values from summaries v2
#'
#' given a data.frame of means and variances, use a two step non-linear least squares.
#' The first step is done on the mean vs sd, then the estimates are used in a second
#' that estimates them using the mean vs rsd.
#'
#' @param data data.frame of means and variances
#' @param ... other nls parameters
#'
#' @import dplyr
#' @export
#' @return vector
calc_sd_rsd_nls <- function(data, ...){
  nl_sd <- filter(data, type == "sd") %>% nls(var ~ B + A*mean, data = ., start = list(A = 0, B = 0), ...)
  nl_rsd <- filter(data, type == "rsd") %>% nls(var ~ ((A * mean) + B) / mean, data = ., start = list(A = coef(nl_sd)["A"], B = coef(nl_sd)["B"]), ...)

  return(c(additive = coef(nl_rsd)["B"], proportional = coef(nl_rsd)["A"]))
}
