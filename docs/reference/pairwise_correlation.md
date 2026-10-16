# pairwise correlation

given a data matrix, calculate *inter-row* correlation values.

## Usage

``` r
pairwise_correlation(
  data_matrix,
  use = "pairwise.complete.obs",
  exclude_na = TRUE,
  exclude_inf = TRUE,
  exclude_0 = FALSE,
  zero_value = 0,
  method = "pearson"
)
```

## Arguments

- data_matrix:

  the data

- use:

  what data to use, see `cor`

- exclude_na:

  should NA values be excluded (default TRUE)

- exclude_inf:

  should Inf values be excluded (default TRUE)

- exclude_0:

  should 0 values be excluded (default FALSE)

- zero_value:

  what value represents zero (default is 0)

- method:

  which method of correlation to use

- weight:

  should the correlations be weighted by information content and
  correspondence?

- not_both:

  should correspondence include things FALSE in both?

## Value

list

## Details

The function returns a named list with:

- cor:

  the correlation matrix

- keep:

  a logical matrix indicating which points passed filtering
