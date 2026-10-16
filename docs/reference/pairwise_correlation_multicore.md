# pairwise correlation multicore

given a data matrix, calculate inter-row correlation values

## Usage

``` r
pairwise_correlation_multicore(
  data_matrix,
  use = "pairwise",
  exclude_na = TRUE,
  exclude_inf = TRUE,
  exclude_0 = FALSE,
  method = "pearson"
)
```

## Arguments

- data_matrix:

  the data

- use:

  what data to use, "pairwise" or "complete"

- exclude_na:

  should NA values be excluded (default TRUE)

- exclude_inf:

  should Inf values be excluded (default TRUE)

- exclude_0:

  should 0 values be excluded (default FALSE)

- method:

  which method of correlation to use

## Value

matrix
