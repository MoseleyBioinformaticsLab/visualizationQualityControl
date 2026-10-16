# fraction of outliers

Calculates the fraction of entries in each sample that are more than `X`
standard deviations from the trimmed mean. See Details.

## Usage

``` r
outlier_fraction(
  data,
  sample_classes = NULL,
  n_trim = 3,
  n_sd = 5,
  remove_missing = NA
)
```

## Arguments

- data:

  the data matrix (samples are columns, rows are features)

- sample_classes:

  the sample classes

- n_trim:

  how many features to trim at each end (default is 3)

- n_sd:

  how many SD before treated as outlier (default is 5)

- remove_missing:

  what missing values be removed before calculating? (default is NA)

## Value

data.frame

## Details

Based on the Gerlinski paper
[link](https://dx.doi.org/10.1093/bioinformatics/btv425) for each
feature (in a sample class), take the range across all the samples,
remove the `n_trim` lowest and highest values, and calculate the `mean`
and `sd`, and the actual upper and lower ranges of `n_sd` from the
`mean`. For each sample and feature, determine if *within* or *outside*
that limit. Fraction is reported as the number of features outside the
range.

Returns a `data.frame` with:

- sample_id:

  the sample id, `rownames` are used if available, otherwise this is an
  index

- sample_class:

  the class of the sample if `sample_classes` were provided, otherwise
  given a default of "C1"

- frac:

  the actual outlier fraction calculated for that sample
