# run pca

Given a matrix of values, runs principal component analysis (PCA) using
`prcomp`.

## Usage

``` r
visqc_run_pca(data_matrix, center = TRUE, scale = FALSE, ...)
```

## Arguments

- data_matrix:

  matrix of values

- center:

  should the data be centered (TRUE)

- scale:

  should the data be scaled to unit variance (FALSE)

- ...:

  other arguments to `prcomp`

## Value

`visqc_run_pca` returns a list containing

- `scores`:

  a data.frame of scores, each column is a principal component.

- `loadings`:

  a data.frame of the loadings, each column is a principal component.

- `variance`:

  a data.frame of variance information, each row is a principal
  component.

## Details

runs PCA via prcomp, turns the scores into a data.frame, labelled with
their percent variance contribution if `labelled` is available, and
returns the scores, the loadings, and the variance contributions.

## See also

[stats::prcomp](https://rdrr.io/r/stats/prcomp.html)
