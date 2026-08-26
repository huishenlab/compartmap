# Compute differential compartments based on Mahalanobis distance

Compute differential compartments based on Mahalanobis distance

## Usage

``` r
get_dc(x, cov_method = c("base", "robust", "mcd"), alpha_level = 0.05)
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
