# determine outliers

determine outliers

## Usage

``` r
determine_outliers(
  median_correlations = NULL,
  outlier_fraction = NULL,
  cor_weight = 1,
  frac_weight = 1,
  only_high = TRUE
)
```

## Arguments

- median_correlations:

  median correlations

- outlier_fraction:

  outlier fractions

- cor_weight:

  how much weight for the correlation score?

- frac_weight:

  how much weight for the outlier fraction?

- only_high:

  should only things at the low end of score be removed?

## Value

data.frame

## Details

For outlier sample detection, one should first generate median
correlations using `median_correlations`, and outlier fractions using
`outlier_fraction`. If you only have one or the other, than you should
use named arguments to only pass the one or the other.

Alternatively, you can change the weighting used for median correlations
or outlier fraction, including setting them to 0.
