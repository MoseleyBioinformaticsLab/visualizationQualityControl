# summarize data

summarizes a matrix or data.frame, where columns are samples and rows
are features

## Usage

``` r
summarize_data(
  in_data,
  sample_classes = NULL,
  avg_function = mean,
  log_transform = FALSE,
  remove_missing = NA
)
```

## Arguments

- in_data:

  matrix or data.frame

- sample_classes:

  which samples are in which class

- avg_function:

  which function to use for summary

- log_transform:

  apply a log-transform to the mean

- remove_missing:

  remove missing values before summarizing

## Value

data.frame
