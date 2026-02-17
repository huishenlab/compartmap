# p_tilde in AC

p_tilde in AC

## Usage

``` r
.p_approx(n1, n0, q)
```

## Arguments

- n1:

  number of successes/ones

- n0:

  number of failures/zeroes

- q:

  quantile for eventual CI (e.g. 0.95 for a 95 percent binomial CI)

## Value

the approximate success probability for a smoothed CIs

## Details

\\\tilde{p} = \frac{1}{\tilde{n}}(n\_{\text{success}} +
\frac{z^2\_\alpha}{2})\\
