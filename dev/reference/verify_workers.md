# Verify that requested thread count is not higher than available

Verify that requested thread count is not higher than available

## Usage

``` r
verify_workers(n_workers)
```

## Arguments

- n_workers:

  The number of workers to check availability

## Value

TRUE if the requested `thread_count` does not exceed available resources
as defined by
[`parallelly::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)
