# information volume

calculates the information volume weight

## Usage

``` r
information_volume(in_x, in_y)
```

## Arguments

- in_x:

  which things are in X

- in_y:

  which things are in Y

## Value

numeric

## Details

The information volume is based on how many things are actually not
missing in both samples. For weighting of correlations, this value will
be scaled by the maximum observed value across calculated correlations.
