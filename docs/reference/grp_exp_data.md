# grp_exp_data

Example data that requires log-transformation before doing PCA or other
QC. A list with 2 named entries:

## Usage

``` r
data(grp_exp_data)
```

## Format

List with 2 entries, data and class

## Source

Robert M Flight

## Details

- data:

  a data matrix with 1000 rows and 20 columns

- class:

  a character vector of 20 entries denoting classes

The data comes from two groups of samples, where there is ~0.80
correlation within each group, and ~0.38 correlation between groups.
