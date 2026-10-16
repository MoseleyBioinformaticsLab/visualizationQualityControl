# test loadings

Given a matrix of loadings for principal components, and a set of
components to test, for each loading in each component, generates a null
distribution from the other loadings in all the other components, and
reports a p-value for that loading.

## Usage

``` r
visqc_test_pca_loadings(
  loadings,
  test_columns,
  progress = FALSE,
  direction = FALSE
)
```

## Arguments

- loadings:

  matrix of loadings from pca decomposition

- test_columns:

  names of the columns of the loadings to test

- progress:

  should progress be reported

- direction:

  should direction of loading be tested?

## Value

named list
