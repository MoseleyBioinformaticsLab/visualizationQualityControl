# easier heatmaps

rolls some of the common `Heatmap` options into a single function call
to make life easier when creating lots of heatmaps. **Note:** clustering
of rows and columns is disabled, it is expected that you are reordering
the matrix beforehand, or passing in `column_order` and `row_order` as
arguments to be passed to `Heatmap` (see example). Matrices can be
reordered using
[`similarity_reorderbyclass`](similarity_reorderbyclass.md), and nice
class colors generated using
[`generate_group_colors`](generate_group_colors.md)

## Usage

``` r
visqc_heatmap(
  matrix_data,
  color_values,
  title = "",
  row_color_data = NULL,
  row_color_list = NULL,
  col_color_data = NULL,
  col_color_list = NULL,
  ...
)
```

## Arguments

- matrix_data:

  the matrix you want to plot as a heatmap

- color_values:

  the color mapping of values to colors (see Details)

- title:

  what do the values represent

- row_color_data:

  data for row annotations

- row_color_list:

  list for row annotations

- col_color_data:

  data for column annotations

- col_color_list:

  list for column annotations

- ...:

  other `Heatmap` parameters

## Details

This function uses the `ComplexHeatmap` package to produce heatmaps with
complex row- and column-color annotations. Both `row_color_data` and
`col_color_data` should be `data.frame`'s where each column describes
meta-data about the rows or columns of the matrix. The `row_color_list`
and `col_color_list` provide the mapping of color to annotation, where
each `list` entry should be a named vector of colors, with the list
entry corresponding to a column entry in the data.frame, and the names
of the colors corresponding to annotations in that column.

## Examples

``` r
 
if (FALSE) { # \dontrun{
library(circlize)
data(grp_cor)
data(grp_info)
colormap <- colorRamp2(c(0, 1), c("black", "white"))

annotation_color <- c(grp1 = "green", grp2 = "red", set1 = "blue",
                      set2 = "yellow")

row_data <- grp_info[, "grp", drop = FALSE]
col_data <- grp_info[, "set", drop = FALSE]
row_annotation = list(grp = annotation_color[1:2])
col_annotation = list(set = annotation_color[3:4])

visqc_heatmap(grp_cor, colormap, row_color_data = row_data, row_color_list = row_annotation,
                 col_color_data = col_data, col_color_list = col_annotation)
                 
reorder_sim <- similarity_reorderbyclass(grp_cor, transform = "sub_1")
visqc_heatmap(grp_cor, colormap, "reorder1", row_data, row_annotation, col_data, col_annotation,
                 column_order = reorder_sim$indices, row_order = reorder_sim$indices)

sample_classes <- grp_info[, "grp", drop = FALSE]
reorder_sim2 <- similarity_reorderbyclass(grp_cor, sample_classes, "sub_1")
visqc_heatmap(grp_cor, colormap, "reorder2", row_data, row_annotation, col_data, col_annotation,
                 column_order = reorder_sim2$indices, row_order = reorder_sim2$indices)
} # }
```
