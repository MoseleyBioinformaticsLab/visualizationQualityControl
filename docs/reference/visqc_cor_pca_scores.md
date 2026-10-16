# correlate scores and outcome

Given a matrix of PCA scores, set of sample attributes to test, goes
through and performs an ICI-Kt of the scores versus the attribute.

## Usage

``` r
visqc_cor_pca_scores(pca_scores, sample_info)
```

## Arguments

- pca_scores:

  the scores data.frame or matrix to test

- sample_info:

  data.frame of sample attributes to test

  Important: All of the attributes must be numeric, or character. If
  character, they will be transformed to a factor, and the numeric
  factor levels will be used instead. If missing values are present,
  that is OK, as long as they are missing-not-at-random (i.e. missing at
  the low end of the values).

## Value

data.frame
