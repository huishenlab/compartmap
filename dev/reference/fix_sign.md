# Correct the sign of the compartment call vector

Correct the sign of the compartment call vector

## Usage

``` r
fix_sign(x, na.rm = FALSE)
```

## Arguments

- x:

  A `CompartmentCall` object

- na.rm:

  Whether to remove NAs. The presence of NAs can reduce the accuracy of
  the sign flip as the open vs closed gene densities cannot be computed
  correctly
