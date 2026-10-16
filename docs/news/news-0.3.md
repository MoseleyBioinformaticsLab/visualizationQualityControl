# Version 0.3

## visualizationQualityControl 0.3.0

- `median_correlations` gains a new argument, `between_classes` to
  generate the median values to samples in other classes. This causes
  the appearance of two more columns when set to TRUE. The default is
  FALSE, so hopefully this does not cause current code to misbehave, but
  I’ve bumped the version number as a warning.

## visualizationQualityControl 0.3.2

- `keep_non_zero_percentage` gains an argument, `all`, that defaults to
  `FALSE` to keep previous behavior. Setting `all = TRUE` means that the
  value must be non-zero in at least X% of **all** of the sample
  classes.

## visualizationQualityControl 0.3.16

- Removed requirement for `ggbiplot`, instead we added a function for
  calculating the variances of each of the PCs in the scores.

- updated the vignette accordingly.

- Now using `globally_it_weighted_correlation` and
  `locally_it_weighted_correlation` instead of `pairwise_correlation`.

## visualizationQualityControl 0.3.85

- Added a function for calculating the information-content-informed
  Kendall-tau correlation, `ici_kendallt`, and variants around
  calculating all pairwise correlations between samples;
  `visqc_ici_kendallt` and `visqc_ici_kendallt_splitup` for parallel
  processing.

## visualizationQualityControl 0.3.96

- Now throw an error if X and Y are not the same length in
  `ici_kendallt`.

## visualizationQualityControl 0.3.100

- Making the splitup version of ICI-Kendall-tau the “implementation”
  (`visqc_ici_kendallt`), and using a single core if the user doesn’t
  setup a “plan” first. A reference version still exists so we can run
  tests against it, but it is no longer exported for general users.

- Also inlined the C++ sign function, which gave us another 3X speedup
  on my 8 core machine on a larger test data set.
