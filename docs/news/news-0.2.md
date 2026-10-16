# Version 0.2

## visualizationQualityControl 0.2.1

- `pairwise_correlation` now uses `cor` internally directly, whereas
  previously it did a `for` loop to allow pairwise comparisons. This
  makes the correlations 3x faster.

- `count` has been removed from the list returned by
  `pairwise_correlation`

- new function `pairwise_correlation_count` to get the counts in each
  pairwise comparison

## visualizationQualityControl 0.2.3

- A bug was discovered in `median_correlations` that meant the wrong
  sample ids might be added to the output data, making detection of real
  problems difficult

## visualizationQualityControl 0.2.5

- Added two functions, `information_volume` and `correspondence` to
  calculate weights based on the amount of things that are non-zero in
  both things when doing pairwise correlation.

- Added logical argument `weight` to `pairwise_correlation` to weight
  the correlations. If `weight = TRUE`, the diagonal will not be **1**
  anymore, but instead will reflect how many features out of the total
  are in that sample.

## visualizationQualityControl 0.2.18

- Augmented correlations (`weight = TRUE`) should be much more useful
  and interpretable.

- `information_volume` and `correspondence` calculations improved.
  Namely that `information_volume` is being scaled by the maximum.

- `correspondence` by default **does not** consider presence of zeros in
  both samples to be informative, this can be changed by setting
  `not_both = TRUE`. The default is more useful in cases where there are
  lots of features and the data is sparse, and zeros are likely to
  happen by chance.

- In addition to returning the `cor` matrix and `keep` matrix,
  `pairwise_correlations` now returns the `raw` correlations, and the
  weighting matrices `info` and `correspondence` so that each one can be
  examined.

- The diagonal of `info` weighting corresponds to how many features a
  sample has compared to the sample with the most features.
