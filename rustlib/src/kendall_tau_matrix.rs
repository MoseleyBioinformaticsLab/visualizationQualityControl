use crate::kendall_tau::Perspective;

/// information-content-informed kendall tau
/// Given a data-matrix, computes the information-content-informed (ICI) Kendall-tau-b between
/// all samples.
///
/// - `data_matrix`: samples are rows, features are columns
/// - `exclude_na`: should NA values be treated as NA?
/// - `exclude_inf`: should Inf values be treated as NA?
/// - `exclude_0`: should zero values be treated as NA?
/// - `zero_value`: what is the actual zero value?
/// - `perspective`: how to treat missing data in denominator and ties, see details
/// - `scale_max`: should everything be scaled compared to the maximum correlation?
/// - `diag_not_na`: should the diagonal entries reflect how many entries in the sample were "good"?
///
/// Default call: `function(data_matrix, exclude_na = TRUE, exclude_inf = TRUE,
///     exclude_0 = TRUE, zero_value = 0, perspective = "global", scale_max = TRUE,
///     diag_good = TRUE, progress = FALSE)`
pub fn visqc_ici_kendall_tau(data_matrix,
exclude_na:bool,
exclude_inf: bool,
exclude_0: bool,
zero_value: f64,
perspective: Perspective,
scale_max:bool,
diag_good: bool,
progress: bool){
    unimplemented!()
}

/*
visqc_ici_kendallt = function(data_matrix,
                             exclude_na = TRUE,
                             exclude_inf = TRUE,
                             exclude_0 = TRUE,
                             zero_value = 0,
                             perspective = "global",
                             scale_max = TRUE,
                             diag_good = TRUE,
                             progress = FALSE){

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

  exclude_data = data_matrix
  exclude_data[exclude_loc] = NA
  # set everything to NA and let R take care of it

  if (ncol(data_matrix) > 2 && progress) {
    prog_bar = knitrProgressBar::progress_estimated(ncol(exclude_data) * (ncol(exclude_data))/ 2)
  } else {
    prog_bar = NULL
  }

  cor_matrix = matrix(NA, nrow = ncol(exclude_data), ncol = ncol(exclude_data))
  rownames(cor_matrix) = colnames(cor_matrix) = colnames(exclude_data)
  ntotal = 0
  for (icol in seq(1, ncol(exclude_data))) {
    for (jcol in seq(icol, ncol(exclude_data))) {
      cor_matrix[icol, jcol] = cor_matrix[jcol, icol] = ici_kendallt(exclude_data[, icol], exclude_data[, jcol], perspective = perspective)
      knitrProgressBar::update_progress(prog_bar)
      # ntotal = ntotal + 1
      # message(ntotal)
    }
  }

  # calculate the max-cor value for use in scaling across multiple comparisons
  n_observations = nrow(exclude_data)
  n_na = sort(colSums(exclude_loc))
  m_value = floor(sum(n_na[1:2]) / 2)
  n_m = n_observations - m_value
  max_cor_denominator = choose(n_m, 2) + n_observations * m_value
  max_cor_numerator = choose(n_m, 2) + n_observations * m_value + choose(m_value, 2)
  max_cor = max_cor_denominator / max_cor_numerator

  if (scale_max) {
    out_matrix = cor_matrix / max_cor
  } else {
    out_matrix = cor_matrix
  }

  if (diag_good) {
    n_good = colSums(!exclude_loc)
    diag(out_matrix) = n_good / max(n_good)
  }

  return(list(cor = out_matrix, raw = cor_matrix, keep = t(!exclude_loc)))
}

*/