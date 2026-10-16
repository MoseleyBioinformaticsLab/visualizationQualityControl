# pairwise correlation keep both zero

given a data matrix, calculate inter-row correlation values, but if both
entries have zero, keep it

## Usage

``` r
pairwise_correlation_both0(
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
