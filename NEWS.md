# v 0.3.85

* Added a function for calculating the information-content-informed Kendall-tau
correlation, `ici_kendallt`, and variants around calculating all pairwise correlations
between samples; `visqc_ici_kendallt` and `visqc_ici_kendallt_splitup` for parallel processing.

# v 0.3.16

* Removed requirement for `ggbiplot`, instead we added a function for calculating
the variances of each of the PCs in the scores.

* updated the vignette accordingly.

* Now using `globally_it_weighted_correlation` and `locally_it_weighted_correlation`
instead of `pairwise_correlation`.

# v 0.3.2

* `keep_non_zero_percentage` gains an argument, `all`, that defaults to `FALSE`
to keep previous behavior. Setting `all = TRUE` means that the value must be
non-zero in at least X% of **all** of the sample classes.

# v 0.3.0

* `median_correlations` gains a new argument, `between_classes` to generate the
median values to samples in other classes. This causes the appearance of two
more columns when set to TRUE. The default is FALSE, so hopefully this does not
cause current code to misbehave, but I've bumped the version number as a warning.

# v 0.2.18

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

# v 0.2.5

* Added two functions, `information_volume` and `correspondence` to calculate
weights based on the amount of things that are non-zero in both things when
doing pairwise correlation.

* Added logical argument `weight` to `pairwise_correlation` to weight the correlations. If `weight = TRUE`, the diagonal will not be **1** anymore, but instead will reflect how many features out of the total are in that sample.

# v 0.2.3

* A bug was discovered in `median_correlations` that meant the wrong sample ids
might be added to the output data, making detection of real problems difficult

# v 0.2.1

* `pairwise_correlation` now uses `cor` internally directly, whereas previously
it did a `for` loop to allow pairwise comparisons. This makes the correlations
3x faster.

* `count` has been removed from the list returned by `pairwise_correlation`

* new function `pairwise_correlation_count` to get the counts in each pairwise
comparison
 
# v 0.1.1

* Changed correlation function to return a list instead of a matrix. This
list contains the correlations (`cor`), counts in each correlation (`count`),
and which points passed the criteria (`keep`).
