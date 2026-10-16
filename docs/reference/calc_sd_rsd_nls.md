# calculate values from summaries v2

given a data.frame of means and variances, use a two step non-linear
least squares. The first step is done on the mean vs sd, then the
estimates are used in a second that estimates them using the mean vs
rsd.

## Usage

``` r
calc_sd_rsd_nls(data, ...)
```

## Arguments

- data:

  data.frame of means and variances

- ...:

  other nls parameters

## Value

vector
