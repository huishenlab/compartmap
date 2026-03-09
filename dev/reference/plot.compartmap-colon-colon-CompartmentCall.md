# Plot singular values of a `CompartmentCall` object

Plot singular values of a `CompartmentCall` object

## Usage

``` r
# S3 method for class '`compartmap::CompartmentCall`'
plot(
  x,
  ...,
  type = "line",
  label_coords = FALSE,
  res = "kb",
  width = 0.5,
  ylim = c(-1, 1)
)
```

## Arguments

- x:

  `CompartmentCall` or `CompartmapCall` object

- ...:

  Placeholder for the `plot` generic - arguments have not effect

- type:

  Whether to plot the singular values as `"line"` or `"bar"` plots.

- label_coords:

  Label the x-axis with genomic coordinates. Uses a numeric index when
  set to `FALSE`. Using coordinate labels can severely crowd the x-axis,
  especially with Kb-resolution calls.

- res:

  The resolution to round the genomic coordinates to (kilobase: "kb" or
  megabase: "mb")

- width:

  The width of the `geom_line` if `type = "line"` or the width of the
  bar if `type = "bar"` in the plot

- ylim:

  Upper and lower bound for the y-axis
