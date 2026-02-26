# James-Stein estimator

James-Stein estimator

## Usage

``` r
.jse(x, grand.mean = NULL, targets = NULL)
```

## Arguments

- x:

  input vector of binned means across samples

- grand.mean:

  The global mean across samples

- targets:

  Samples to shrink towards

  \\\hat{\theta}\_{JS+} = \left(1 - \frac{(m -
  3)\sigma^2}{\|\|\textbf{y} - \nu\|\|^2}\right)\\
