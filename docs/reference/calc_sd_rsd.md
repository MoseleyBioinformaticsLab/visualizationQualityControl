# calculate values from summaries

given a data.frame of means and variances, calculate mean sd at low end
and mean rsd at high end.

## Usage

``` r
calc_sd_rsd(data, low_cut, hi_cut = NULL)
```

## Arguments

- data:

  data.frame of means and variances

- low_cut:

  means \<= this value used for average sd

- hi_cut:

  means \>= this value used for average rsd

## Value

vector
