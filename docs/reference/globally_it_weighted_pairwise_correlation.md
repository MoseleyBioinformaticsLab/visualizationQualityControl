# globally information theoretic weighted pairwise correlation

given a data matrix, calculate *inter-row* correlation values that are
weighted by information content and information consistency.

## Usage

``` r
globally_it_weighted_pairwise_correlation(
  data_matrix,
  use = "pairwise.complete.obs",
  exclude_na = TRUE,
  exclude_inf = TRUE,
  exclude_0 = FALSE,
  zero_value = 0,
  method = "pearson",
  scale_information_content = TRUE
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

- scale_information_content:

  should the information content be scaled by the maximum value?

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

- content:

  the weighting calculated from information content (weighted jaccard)

- consistency:

  weighted information consistency
