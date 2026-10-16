# test scores and outcome

Given a matrix of PCA scores, set of sample attributes to test, goes
through and performs an ANOVA of the scores versus the attribute.

## Usage

``` r
visqc_test_pca_scores(pca_scores, sample_info)
```

## Arguments

- pca_scores:

  matrix or data.frame of scores from a PCA decomposition

- sample_info:

  data.frame of sample attributes to test

## Value

data.frame

## Details

For each variable being tested, the function first checks the number of
unique and missing values in that variable. If they are all missing, or
the unique values is 1, or equal to the number of rows, then the
variable is not tested and NA values will be returned for that variable.
