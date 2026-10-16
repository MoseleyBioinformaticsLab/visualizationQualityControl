# reorder by sample class

to avoid spurious visualization problems, it is useful in a heatmap
visualization to reorder the samples within each sample class. This
function uses hierarchical clustering and dendsort to sort entries in a
distance matrix.

## Usage

``` r
similarity_reorderbyclass(
  similarity_matrix,
  sample_classes = NULL,
  transform = "none",
  hclust_method = "complete",
  dendsort_type = "min"
)
```

## Arguments

- similarity_matrix:

  matrix of similarities between objects

- sample_classes:

  data.frame or factor denoting classes

- transform:

  a transformation to apply to the data

- hclust_method:

  which method for clustering should be used

- dendsort_type:

  how should dendsort do reordering?

## Value

a list containing the reordering of the matrix in a:

1.  dendrogram

2.  numeric vector

3.  character vector (will be NULL if rownames are not set on the
    matrix)

## Details

The `similarity_matrix` should be either a square matrix of similarity
values or a distance matrix of class `dist`. If your matrix does not
encode a "true" distance, you can use a `transform` to turn it into a
true `distance` (for example, if you have correlation, then a distance
would be 1 - correlation, use "sub_1" as the transform argument).

The `sample_classes` should be either a data.frame or factor argument.
If a data.frame is passed, all columns of the data.frame will be pasted
together to create a factor for splitting the data into groups. If the
rownames of the data.frame do not correspond to the rownames or colnames
of the matrix, then it is assumed that the ordering in the matrix and
the data.frame are identical.

## Examples

``` r
library(visualizationQualityControl)
set.seed(1234)
mat <- matrix(rnorm(100, 2, sd = 0.5), 10, 10)
rownames(mat) <- colnames(mat) <- letters[1:10]
neworder <- similarity_reorderbyclass(mat)
mat[neworder$indices, neworder$indices]
#>           a         d        g        b        c        f        e        h
#> a 1.3964671 2.5511488 2.328294 1.761404 2.067044 1.096984 2.724748 2.003446
#> d 0.8271511 1.7493710 1.665183 2.032229 2.229795 1.492519 1.859688 2.324143
#> g 1.7126300 0.9099802 1.430696 1.744495 2.287378 2.823909 1.446341 1.304650
#> b 2.1387146 1.7622035 3.274496 1.500807 1.754657 1.708962 1.465679 1.772266
#> c 2.5422206 1.6452800 1.982620 1.611873 1.779726 1.445555 1.572318 1.816738
#> f 2.2530279 1.4161904 2.888542 1.944857 1.275898 2.281528 1.515743 1.923301
#> e 2.2145623 1.1854533 1.996198 2.479747 1.653140 1.918845 1.502830 3.035135
#> h 1.7266841 1.3295034 2.683914 1.544402 1.488172 1.613323 1.374007 1.638209
#> i 1.7177740 1.8528531 2.664782 1.581414 1.992431 2.802955 1.738086 2.129131
#> j 1.5549811 1.7670512 2.168236 3.207918 1.532026 1.421096 1.751575 1.841470
#>          i        j
#> a 1.911105 1.973421
#> d 1.913106 2.500757
#> g 2.274999 1.432696
#> b 1.915003 2.127598
#> c 1.313849 2.852982
#> f 2.348804 2.177775
#> e 2.425116 1.752208
#> h 1.798634 2.439102
#> i 1.904203 2.486458
#> j 1.402736 3.060559

sample_class <- data.frame(grp = rep(c("grp1", "grp2"), each = 5), stringsAsFactors = FALSE)
rownames(sample_class) <- rownames(mat)
neworder2 <- similarity_reorderbyclass(mat, sample_class[, "grp", drop = FALSE])

# if there is a class with only one member, it is dropped, with a warning
sample_class[10, "grp"] = "grp3"
neworder3 <- similarity_reorderbyclass(mat, sample_class[, "grp", drop = FALSE])
#> Warning: Removing groups: grp3
neworder3$indices # 10 should be missing
#> [1] 1 4 5 2 3 6 8 7 9

mat[neworder2$indices, neworder2$indices]
#>           a         d        e        b        c        i        j        g
#> a 1.3964671 2.5511488 2.724748 1.761404 2.067044 1.911105 1.973421 2.328294
#> d 0.8271511 1.7493710 1.859688 2.032229 2.229795 1.913106 2.500757 1.665183
#> e 2.2145623 1.1854533 1.502830 2.479747 1.653140 2.425116 1.752208 1.996198
#> b 2.1387146 1.7622035 1.465679 1.500807 1.754657 1.915003 2.127598 3.274496
#> c 2.5422206 1.6452800 1.572318 1.611873 1.779726 1.313849 2.852982 1.982620
#> i 1.7177740 1.8528531 1.738086 1.581414 1.992431 1.904203 2.486458 2.664782
#> j 1.5549811 1.7670512 1.751575 3.207918 1.532026 1.402736 3.060559 2.168236
#> g 1.7126300 0.9099802 1.446341 1.744495 2.287378 2.274999 1.432696 1.430696
#> f 2.2530279 1.4161904 1.515743 1.944857 1.275898 2.348804 2.177775 2.888542
#> h 1.7266841 1.3295034 1.374007 1.544402 1.488172 1.798634 2.439102 2.683914
#>          f        h
#> a 1.096984 2.003446
#> d 1.492519 2.324143
#> e 1.918845 3.035135
#> b 1.708962 1.772266
#> c 1.445555 1.816738
#> i 2.802955 2.129131
#> j 1.421096 1.841470
#> g 2.823909 1.304650
#> f 2.281528 1.923301
#> h 1.613323 1.638209
cbind(neworder$names, neworder2$names)
#>       [,1] [,2]
#>  [1,] "a"  "a" 
#>  [2,] "d"  "d" 
#>  [3,] "g"  "e" 
#>  [4,] "b"  "b" 
#>  [5,] "c"  "c" 
#>  [6,] "f"  "i" 
#>  [7,] "e"  "j" 
#>  [8,] "h"  "g" 
#>  [9,] "i"  "f" 
#> [10,] "j"  "h" 
```
