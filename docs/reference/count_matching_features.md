# matching features

For a given set of feature-sample matrices, calculates how many features
are in that sample, as well as in common to all other samples, and if
provided, in common within and outside the same sample group.

## Usage

``` r
count_matching_features(feature_matrix, zero_value = NA, groups = NULL)
```

## Arguments

- feature_matrix:

  the feature to sample matrix.

- zero_value:

  what is the zero value? Default is NA

- groups:

  what are the groups

## Value

data.frame
