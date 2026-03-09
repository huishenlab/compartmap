# Filter to bins with call values greater than or equal to a threshold value

Filter to bins with call values greater than or equal to a threshold
value

## Usage

``` r
filter(x, threshold = 0.02)
```

## Arguments

- x:

  A `CompartmentCall` object

- threshold:

  The absolute value to use for filtering. Rows where any value is less
  than this threshold are dropped
