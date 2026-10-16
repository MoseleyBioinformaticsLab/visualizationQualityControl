# keep features with percentage of non-zeros

Given a value matrix (features are columns, samples are rows), and
sample classes, find those things that are not zero in at least a
certain number of one of the classes, and keep them

## Usage

``` r
filter_non_zero_percentage(data_matrix, sample_classes = NULL, keep_num = 0.75)
```

## Arguments

- data_matrix:

  the matrix of values to work with

- sample_classes:

  the classes of each sample

- keep_num:

  what number of samples in each class need a non-zero value (see
  Details)

## Value

matrix

## Details

This function is being deprecated and all code should use the
`keep_non_zero_percentage` function instead.

## See also

keep_non_zero_percentage
