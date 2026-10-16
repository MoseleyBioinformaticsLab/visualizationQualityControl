# calculate median correlations

Given a correlation matrix and optionally the sample class information,
calculates the median correlations of each sample to all other samples
in the same class. May be useful for determining outliers.

## Usage

``` r
median_correlations(cor_matrix, sample_classes = NULL, between_classes = FALSE)
```

## Arguments

- cor_matrix:

  the sample - sample correlations

- sample_classes:

  the sample classes as a character or factor

- between_classes:

  should the between class correlations be evaluated?

## Value

data.frame

## Details

The data.frame has 6 columns:

- med_cor:

  the median correlation with other samples

- sample_id:

  the sample id, either the rowname or an index

- sample_class:

  the class of the sample. If not provided, set to "C1"

- compare_class:

  the class of the other sample

- plot_class:

  `sample_class::compare_class` for easy grouping

- group:

  whether the median is *within* or *between* groups of samples
