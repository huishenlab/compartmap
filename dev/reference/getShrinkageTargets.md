# Get the specified samples to shrink towards instead of the global mean

Get the specified samples to shrink towards instead of the global mean

## Usage

``` r
getShrinkageTargets(obj, group)
```

## Arguments

- obj:

  Input matrix

- group:

  Sample names, colnames or column indices to use for targeted shrinkage

## Value

A matrix composed of samples to shrink towards

## Examples

``` r
dummy <- matrix(rnorm(1000), ncol=25)
dummy.sub <- getShrinkageTargets(dummy, group = c(1,5,8,10))
```
