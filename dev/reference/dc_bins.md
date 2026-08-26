# Get bins with significant Mahalanobis distances

Get bins with significant Mahalanobis distances

## Usage

``` r
dc_bins(gr, md, alpha_level = 0.05)
```

## Arguments

- gr:

  Matrix of singular values

- md:

  data.table output from
  [`compute_mahalanobis()`](https://huishenlab.github.io/compartmap/dev/reference/compute_mahalanobis.md)

- alpha_level:

  Significance level to use (default: 0.05)
