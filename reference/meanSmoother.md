# Windowed mean smoother

Windowed mean smoother

## Usage

``` r
meanSmoother(mat, k = 1, iters = 2, delta = 0, weights = NULL)
```

## Arguments

- mat:

  Input data matrix: samples are columns, regions/loci are rows

- k:

  Number of windows to use (default = 1 to smooth with 1 window on
  either side of a position)

- iters:

  Number of iterations to smooth (default is 2)

- delta:

  Convergence threshhold (overrides iter if \> 0; default is 0)

- weights:

  Weights, if using any (NULL)

## Value

     Smoothed data matrix

## Examples

``` r
dummy <- matrix(rnorm(10000), ncol = 25)
smooth.dummy <- meanSmoother(dummy)
smooth.dummy <- meanSmoother(dummy, iters = 3)
smooth.dummy <- meanSmoother(dummy, delta = 1e-3)
```
