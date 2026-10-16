# pairwise distance

given a data matrix, calculate inter-row distances

## Usage

``` r
pairwise_distance(
  data_matrix,
  use = "pairwise",
  exclude_na = TRUE,
  exclude_inf = TRUE,
  exclude_0 = FALSE
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

## Value

matrix
