# Normalize quantiles From limma version 3.64.3

Normalize quantiles From limma version 3.64.3

## Usage

``` r
normalize_quantiles(A, ties = TRUE)
```

## Arguments

- A:

  Numeric matrix of compartment call singular values

- ties:

  logical. If TRUE, ties in each column of A are treated in careful way.
  tied values will be normalized to the mean of the corresponding pooled
  quantiles.
