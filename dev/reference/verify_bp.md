# Verify that the input BiocParallelParam is valid

Verify that the input BiocParallelParam is valid

## Usage

``` r
verify_bp(bp)
```

## Arguments

- bp:

  A BiocParallelParam or list of 2 BiocParallelParam objects

## Value

TRUE if the total `bpnworkers` in the input does not exceed available
resources as defined by
[`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)
