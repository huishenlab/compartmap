# Show differential compartments as overlaid `ggplot2::geom_rect()` over significant distances

Show differential compartments as overlaid
[`ggplot2::geom_rect()`](https://ggplot2.tidyverse.org/reference/geom_tile.html)
over significant distances

## Usage

``` r
dc_plotter(
  ccall_pd,
  md,
  select,
  type = c("line", "bar"),
  show_md = TRUE,
  alpha_level = 0.05,
  shade = "maroon",
  alpha = 0.5,
  ylim = c(-0.5, 0.5),
  xlim = NULL,
  label_ids = TRUE,
  color = NULL,
  fill = NULL,
  linewidth = 0.5
)
```

## Arguments

- ccall_pd:

  The `@df` slot of a `MultiCompartmapCall` object

- md:

  data.table output from
  [`compute_mahalanobis()`](https://huishenlab.github.io/compartmap/dev/reference/compute_mahalanobis.md)

- select:

  Column name or index to compute differential compartments on Set this
  to plot the compartment scores of all input columns but show
  differential compartments based only on the provided columns.

- show_md:

  Whether to plot the Mahalanobis distance

- alpha_level:

  Significance level to use (default: 0.05)

- shade:

  Color of the `geom_rect` used to shade differential bins

- alpha:

  Transparency of the `geom_rect`

- ylim:

  The y-axis limits

- xlim:

  The x-axis limits

- label_ids:

  Whether to differential bin ID labels

- color:

  The line colors used for the each sample

- fill:

  The colors used for positive and negative values in bar plots

- linewidth:

  The width of lines in the line plot

## Examples

``` r
v <- c(1:4, 5, 7:9)
get_sequential_idx(v)
#>    start_idx end_idx dc_id
#>        <num>   <num> <int>
#> 1:         1       5     1
#> 2:         1       5     1
#> 3:         1       5     1
#> 4:         1       5     1
#> 5:         1       5     1
#> 6:         7       9     2
#> 7:         7       9     2
#> 8:         7       9     2
```
