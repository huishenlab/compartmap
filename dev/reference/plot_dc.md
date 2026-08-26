# Plot differential compartments based on Mahalanobis distance

Plot differential compartments based on Mahalanobis distance

## Usage

``` r
plot_dc(
  x,
  type = c("line", "bar"),
  cov_method = c("base", "robust", "mcd"),
  alpha_level = 0.05,
  show_md = TRUE,
  shade = "maroon",
  alpha = 0.5,
  ylim = c(-0.1, 0.1),
  xlim = NULL,
  label_ids = TRUE,
  select = NULL,
  color = NULL,
  fill = NULL,
  linewidth = 0.5
)
```

## Arguments

- x:

  A `MultiCompartmentCall` object

- cov_method:

  Method to compute covariance. "base":
  [`stats::cov`](https://rdrr.io/r/stats/cor.html),

- alpha_level:

  Significance level to use (default: 0.05)

- show_md:

  Whether to plot the Mahalanobis distance

- shade:

  Color of the `geom_rect`

- alpha:

  Transparency of the `geom_rect`

- ylim:

  The y-axis limits

- xlim:

  The x-axis limits

- label_ids:

  Whether to differential bin ID labels

- select:

  Column name or index to compute differential compartments on. Set this
  to plot the compartment scores of all object columns but show
  differential compartments based only on the provided columns.

- color:

  The line colors used for the each sample

- fill:

  The colors used for positive and negative values in bar plots

- linewidth:

  The width of lines in the line plot
