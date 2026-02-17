# Plot a denoised correlation matrix

Plot a denoised correlation matrix

## Usage

``` r
plotCorMatrix(
  denoised.cor.mat,
  midpoint = 0.3,
  return.plot.obj = FALSE,
  uppertri = FALSE,
  lowertri = FALSE
)
```

## Arguments

- denoised.cor.mat:

  The denoised correlation matrix object from getDenoisedMatrix

- midpoint:

  The midpoint for the coloring (default is 0.3)

- return.plot.obj:

  Whether to return the ggplot object

- uppertri:

  Whether to keep the upper triangle of the matrix

- lowertri:

  Whether to keep the lower triangle of the matrix

## Value

Either a ggplot object or plot

## Examples

``` r
dummy <- matrix(rnorm(10000), ncol=25)
set.seed(1000)
my_plot <- plotCorMatrix(dummy, return.plot.obj = TRUE)
```
