# Get the Mahalanobis distances for the compartment scores

Get the Mahalanobis distances for the compartment scores

## Usage

``` r
compute_mahalanobis(
  mat,
  cov_method = c("base", "robust", "mcd"),
  alpha_level = 0.05,
  fdr = TRUE
)
```

## Arguments

- mat:

  Matrix of singular values

- cov_method:

  Method to compute covariance. "base":
  [`stats::cov`](https://rdrr.io/r/stats/cor.html), "robust":
  [`robust::covRob()`](https://rdrr.io/pkg/robust/man/covRob.html),
  "mcd": `robust::covRob(estim = "mcd")`

- alpha_level:

  Significance level to use (default: 0.05)

- fdr:

  Whether to perform Benjamini-Hochberg false discovery correction

## Details

The resulting data frame contains the mahalanobis distances and p-values
computed based on the squared mahalanobis distances on the Chi-squared
distribution.
