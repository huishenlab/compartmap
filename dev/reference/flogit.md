# Helper function: squeezed logit

Helper function: squeezed logit

## Usage

``` r
flogit(p, sqz = 0.000001)
```

## Arguments

- p:

  a vector of values between 0 and 1 inclusive

- sqz:

  the amount by which to 'squeeze', default is .000001

## Value

       a vector of values between -Inf and +Inf

## Examples

``` r
p <- runif(n = 1000)
summary(p)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.0002544 0.2430427 0.5285060 0.5059655 0.7759185 0.9997144 

sqz <- 1 / (10**6)
x <- flogit(p, sqz = sqz)
summary(x)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -8.27550 -1.13607  0.11415  0.02439  1.24204  8.15963 

all(abs(p - fexpit(x, sqz = sqz)) < sqz)
#> [1] TRUE
all(abs(p - fexpit(flogit(p, sqz = sqz), sqz = sqz)) < sqz)
#> [1] TRUE
```
