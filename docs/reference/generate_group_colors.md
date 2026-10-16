# create set of disjoint colors

When multiple sample classes need to be visualized on a heatmap, it is
useful to be able to distinguish them by color. This function generates
a set of colors for sample classes

## Usage

``` r
generate_group_colors(n_group, randomize = NULL)
```

## Arguments

- n_group:

  how many groups should there be colors for

- randomize:

  should colors be randomized? (default is `NULL`). See *details*.

## Details

the default for `randomize` is `NULL`, so that reordering the colors
randomly is decided purely based on the number of colors requested.
Currently, that cutoff is *5* colors, less than that the colors will
always be in the same order, for *5* colors or more, they will be in a
scrambled order, different each time unless `set.seed` is used. If
`randomize` is `TRUE` or `FALSE`, then it overrides the defaults.
