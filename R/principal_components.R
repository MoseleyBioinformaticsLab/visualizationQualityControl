#' calculate pca contributions
#'
#' Given a set of PCA scores, calculates their variance contributions, cumulative
#' contributions, and generates a percent label that can be used for labeling plots.
#'
#' @param pca_scores matrix of scores, columns are each PC
#'
#' @export
#' @return data.frame
visqc_score_contributions = function(pca_scores) {
  pca_vars = apply(pca_scores, 2, var)
  pca_vars = data.frame(
    pc = names(pca_vars),
    variance = pca_vars,
    stringsAsFactors = FALSE
  )
  pca_vars = dplyr::mutate(
    pca_vars,
    percent = variance / sum(variance),
    cumulative = cumsum(percent)
  )
  pca_labels = purrr::map2_chr(
    pca_vars$pc,
    pca_vars$percent,
    function(.x, .y) {
      paste0(
        .x,
        " (",
        format(
          .y * 100,
          digits = 2,
          scientific = FALSE,
          justify = "none",
          trim = TRUE
        ),
        "%)"
      )
    }
  )
  pca_vars$labels = pca_labels
  pca_vars
}

#' run pca
#'
#' Given a matrix of values, runs principal component analysis (PCA)
#' using `prcomp`.
#'
#' @param data_matrix matrix of values
#' @param center should the data be centered (TRUE)
#' @param scale should the data be scaled to unit variance (FALSE)
#' @param ... other arguments to `prcomp`
#'
#' @details runs PCA via prcomp, turns the scores into a data.frame,
#' labelled with their percent variance contribution if `labelled` is
#' available, and returns the scores, the loadings, and the
#' variance contributions.
#'
#' @export
#' @seealso [stats::prcomp]
#' @return `visqc_run_pca` returns a list containing
#'   \item{`scores`}{a data.frame of scores, each column is a principal component.}
#'   \item{`loadings`}{a data.frame of the loadings, each column is a principal component.}
#'   \item{`variances`}{a data.frame of variance information, each row is a principal component.}
#'
visqc_run_pca = function(data_matrix, center = TRUE, scale = FALSE, ...) {
  p_results = prcomp(data_matrix, center = center, scale. = scale, ...)

  p_scores = as.data.frame(p_results$x)
  p_loadings = as.data.frame(p_results$rotation)

  p_var = visqc_score_contributions(p_scores)

  if (nzchar(find.package("labelled", quiet = TRUE))) {
    labelled::var_label(p_scores) = p_var$labels
  }

  return(list(scores = p_scores, loadings = p_loadings, variances = p_var))
}
