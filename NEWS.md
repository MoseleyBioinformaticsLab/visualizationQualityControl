# visualizationQualityControl 0.5.3

* `median_correlations` now has more columns, regardless of the setting for `between_classes`, to provide consistency.
* `median_correlations` output gains an extra column, `group` to allow easier filtering of **within** and **between** class comparisons.

# visualizationQualityControl 0.5.2

* Fixed behavior of `visqc_test_pca_scores` to skip testing attributes that either have a single value, or the number of numeric values matches the number of samples.

# visualizationQualityControl 0.5.0

* Changed some functions to treat columns as samples and rows as features:
  * `keep_non_zero_percentage`
  * `summarize_data`
  * `calculate_fratio`
  
* Added `keep_non_missing_percentage`, which allows using multiple values to represent missingnes.
* Made `summarize_data` handle possible missing values.
* Removed correlation calculation functions, those have been superseded by [ICIKendallTau](https://moseleybioinformaticslab.github.io/ICIKendallTau/).

# visualizationQualityControl 0.4.11

* Added a new argument `only_high` to `determine_outliers` to only look at the high end of the score distribution for outliers, as sometimes `boxplot.stats` will pick up outliers at the low end as well.

# visualizationQualityControl 0.4.10

* Updated the quality_control vignette to use ICIKendallTau instead of other correlation measures.

# visualizationQualityControl 0.4.9

* Windows and Mac binaries are now available via r-universe, and installation instructions are updated to reflect that.

# visualizationQualityControl 0.4.7

* Updated determine_outliers to be able to use either the output from median_correlations or outlier_fraction singly or together. 
If using one or the other alone, I suggest explicitly naming the arguments so that the correct entry is set to NULL and the other one used.
* Updated the README to show using ici_kendalltau instead of the it_weighted_correlation.
* Updated tests, and moved to testthat v 3.
* Updated pkgdown for rendering the help site.

# visualizationQualityControl 0.4.0

* Moving all of the ICI-Kendall-tau code into it's own package, [ICIKendallTau](https://moseleybioinformaticslab.github.io/ICIKendallTau). This reduces the dependencies necessary if all you want is to run a fast Kendall-tau.

# visualizationQualityControl 0.3.100

* Making the splitup version of ICI-Kendall-tau the "implementation" (`visqc_ici_kendallt`), and using a single core if the user doesn't setup a "plan" first.
A reference version still exists so we can run tests against it, but it is no longer exported for general users.

* Also inlined the C++ sign function, which gave us another 3X speedup on my 8 core machine on a larger test data set.


# visualizationQualityControl 0.3.96

* Now throw an error if X and Y are not the same length in `ici_kendallt`.

# visualizationQualityControl 0.3.85

* Added a function for calculating the information-content-informed Kendall-tau
correlation, `ici_kendallt`, and variants around calculating all pairwise correlations
between samples; `visqc_ici_kendallt` and `visqc_ici_kendallt_splitup` for parallel processing.

# visualizationQualityControl 0.3.16

* Removed requirement for `ggbiplot`, instead we added a function for calculating
the variances of each of the PCs in the scores.

* updated the vignette accordingly.

* Now using `globally_it_weighted_correlation` and `locally_it_weighted_correlation`
instead of `pairwise_correlation`.

# visualizationQualityControl 0.3.2

* `keep_non_zero_percentage` gains an argument, `all`, that defaults to `FALSE`
to keep previous behavior. Setting `all = TRUE` means that the value must be
non-zero in at least X% of **all** of the sample classes.

# visualizationQualityControl 0.3.0

* `median_correlations` gains a new argument, `between_classes` to generate the
median values to samples in other classes. This causes the appearance of two
more columns when set to TRUE. The default is FALSE, so hopefully this does not
cause current code to misbehave, but I've bumped the version number as a warning.

# visualizationQualityControl 0.2.18

* Augmented correlations (`weight = TRUE`) should be much more useful and interpretable.

* `information_volume` and `correspondence` calculations improved. Namely that
`information_volume` is being scaled by the maximum. 

* `correspondence` by default **does not** consider presence of zeros in both
samples to be informative, this can be changed by setting `not_both = TRUE`. The
default is more useful in cases where there are lots of features and the data is
sparse, and zeros are likely to happen by chance.

* In addition to returning the `cor` matrix and `keep` matrix, `pairwise_correlations`
now returns the `raw` correlations, and the weighting matrices `info` and `correspondence`
so that each one can be examined.

* The diagonal of `info` weighting corresponds to how many features a sample has
compared to the sample with the most features.

# visualizationQualityControl 0.2.5

* Added two functions, `information_volume` and `correspondence` to calculate
weights based on the amount of things that are non-zero in both things when
doing pairwise correlation.

* Added logical argument `weight` to `pairwise_correlation` to weight the correlations. If `weight = TRUE`, the diagonal will not be **1** anymore, but instead will reflect how many features out of the total are in that sample.

# visualizationQualityControl 0.2.3

* A bug was discovered in `median_correlations` that meant the wrong sample ids
might be added to the output data, making detection of real problems difficult

# visualizationQualityControl 0.2.1

* `pairwise_correlation` now uses `cor` internally directly, whereas previously
it did a `for` loop to allow pairwise comparisons. This makes the correlations
3x faster.

* `count` has been removed from the list returned by `pairwise_correlation`

* new function `pairwise_correlation_count` to get the counts in each pairwise
comparison
 
# visualizationQualityControl 0.1.1

* Changed correlation function to return a list instead of a matrix. This
list contains the correlations (`cor`), counts in each correlation (`count`),
and which points passed the criteria (`keep`).
