# locally information theoretic weighted pairwise correlation

given a data matrix, calculate *inter-row* correlation values that are
weighted by information content where only the pairwise information
(jaccard index) is considered.

## Usage

``` r
locally_it_weighted_pairwise_correlation(
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

## Value

list

## Details

The function returns a named list with:

- cor:

  the weighted correlation matrix

- keep:

  a logical matrix indicating which points passed filtering

- raw:

  the raw correlations

- jaccard:

  the jaccard information index used to weight the values
