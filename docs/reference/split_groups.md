# split groups

Given a matrix and a data.frame, character vector or data.frame of
groups, splits the indices / names of the matrix into groups
appropriately. This function assumes that the matrix and the groups are
in the correct order!!

## Usage

``` r
split_groups(in_matrix, groups = NULL)
```

## Arguments

- in_matrix:

  the matrix we want to split up

- groups:

  a data.frame, character vector or factor

## Value

list of groups
