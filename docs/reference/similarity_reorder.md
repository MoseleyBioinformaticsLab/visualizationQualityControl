# cluster and reorder

given a matrix (maybe a distance matrix), cluster and then re-order
using dendsort.

## Usage

``` r
similarity_reorder(
  similarity_matrix,
  matrix_indices = NULL,
  transform = "none",
  hclust_method = "complete",
  dendsort_type = "min"
)
```

## Arguments

- similarity_matrix:

  matrix of similarities

- matrix_indices:

  indices to reorder

- transform:

  should a transformation be applied to the data first

- hclust_method:

  which method for clustering should be used?

- dendsort_type:

  how should the reordering be done? (default is "min")

## Value

a [dendrogram](https://rdrr.io/r/stats/dendrogram.html) object. To get
the order use `order.dendogram`.
