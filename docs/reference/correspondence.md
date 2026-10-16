# correspondence

correspondence

## Usage

``` r
correspondence(in_x, in_y, not_both = FALSE)
```

## Arguments

- in_x:

  which things are in X

- in_y:

  which things are in Y

- not_both:

  also include those things NOT in both

## Value

numeric

## Details

Correspondence is based on how many things agree compared to how many
things do not agree. So by default it is the ratio of

things observed in both / things observed in either

This should handle well most types of data, including sparse data sets.
However, there are cases where the fact that values are missing in both
things being compared counts as information. This alternative method is
triggered by setting `not_both = TRUE`, resulting in the calculation:

(things observed in both + things not observed in both) / total \# of
things

This value will be scaled to the largest observed correspondence for
weighting as well.
