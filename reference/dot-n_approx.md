# n_tilde in AC

n_tilde in AC

## Usage

``` r
.n_approx(n1, n0, q)
```

## Arguments

- n1:

  number of successes/ones

- n0:

  number of failures/zeroes

- q:

  quantile for eventual CI (e.g. 0.95 for a 95 percent binomial CI)

## Value

the effective sample size for smoothed CIs

## Details

\\\tilde{n} = n\_{\text{successes}} + n\_{\text{failures}} +
z^2\_\alpha\\
