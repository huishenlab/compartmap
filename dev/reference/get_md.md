# Compute Mahalanobis distances and pvalues - unfiltered

Compute Mahalanobis distances and pvalues - unfiltered

## Usage

``` r
get_md(x, cov_method = c("base", "robust", "mcd"))
```

## Arguments

- x:

  A `MultiCompartmentCall` object

- cov_method:

  Method to compute covariance. "base":
  [`stats::cov`](https://rdrr.io/r/stats/cor.html), "robust":
  [`robust::covRob()`](https://rdrr.io/pkg/robust/man/covRob.html),
  "mcd": `robust::covRob(estim = "mcd")`

- alpha_level:

  Significance level to use (default: 0.05)
