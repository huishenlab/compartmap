# Agresti-Coull confidence interval for a binomial proportion

Agresti-Coull confidence interval for a binomial proportion

## Usage

``` r
agrestiCoullCI(n1, n0, q)
```

## Arguments

- n1:

  number of successes/ones

- n0:

  number of failures/zeroes

- q:

  quantile for eventual CI (e.g. 0.95 for a 95 percent binomial CI)

## Value

the approximate (q x 100) percent confidence interval for (p\|n1,n0,q)

## Details

\\z\_\alpha = \Phi^{-1}(1 - \frac{\alpha}{2})\\

\\\tilde{n} = n\_{\text{successes}} + n\_{\text{failures}} +
z^2\_\alpha\\

\\\tilde{p} = \frac{1}{\tilde{n}}(n\_{\text{success}} +
\frac{z^2\_\alpha}{2})\\

\\p \approx \tilde{p} \pm z\_\alpha \times
\sqrt{\frac{\tilde{p}}{\tilde{n}} \times (1 - \tilde{p})}\\

## Examples

``` r
agrestiCoullCI(10, 3, 0.95)
#>    conf.est conf.est.lowerCI conf.est.upperCI
#> 1 0.7078205         0.490628        0.9250129
```
