# Version 0.4

## visualizationQualityControl 0.4.0

- Moving all of the ICI-Kendall-tau code into it’s own package,
  [ICIKendallTau](ICIKendallTau::ici_kendalltau). This reduces the
  dependencies necessary if all you want is to run a fast Kendall-tau.

## visualizationQualityControl 0.4.7

- Updated determine_outliers to be able to use either the output from
  median_correlations or outlier_fraction singly or together. If using
  one or the other alone, I suggest explicitly naming the arguments so
  that the correct entry is set to NULL and the other one used.
- Updated the README to show using ici_kendalltau instead of the
  it_weighted_correlation.
- Updated tests, and moved to testthat v 3.
- Updated pkgdown for rendering the help site.
