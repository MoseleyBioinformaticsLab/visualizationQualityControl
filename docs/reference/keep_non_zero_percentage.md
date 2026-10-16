# keep features with percentage of non-zeros

Given a value matrix (features are rows, samples are columns), and
sample classes, find those things that are *not zero* in at least a
certain number of samples in one of the classes, and keep those features
for further processing.

## Usage

``` r
keep_non_zero_percentage(
  data_matrix,
  sample_classes = NULL,
  keep_num = 0.75,
  zero_value = 0,
  all = FALSE
)
```

## Arguments

- data_matrix:

  the matrix of values to work with

- sample_classes:

  the classes of each sample

- keep_num:

  what number of samples in each class need a non-zero value (see
  Details)

- zero_value:

  what number represents zero values

- all:

  is this an either / or OR does it need to be present in all?

## Value

logical

## Details

The number of samples that must be non-zero can be expressed either as a
whole number (that is greater than one), or as a fraction that will be
be multiplied by the number of samples in each class to get the lower
limits for each of the classes.
