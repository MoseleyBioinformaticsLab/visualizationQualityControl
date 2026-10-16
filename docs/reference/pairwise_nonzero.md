# pairwise non-zero

given a data matrix, how many entries are non-zero on a pairwise basis

## Usage

``` r
pairwise_nonzero(
  data_matrix,
  use = "pairwise",
  exclude_na = TRUE,
  exclude_inf = TRUE
)
```

## Arguments

- data_matrix:

  the data

- use:

  how to work with the data

- exclude_na:

  should NA values be excluded (default TRUE)

- exclude_inf:

  should Inf values be excluded (default TRUE)

## Value

matrix
