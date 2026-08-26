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
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> 0.001476 0.259948 0.513944 0.504526 0.745252 0.998696 

sqz <- 1 / (10**6)
x <- flogit(p, sqz = sqz)
summary(x)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -6.516553 -1.046238  0.055792  0.000998  1.073445  6.641115 

all(abs(p - fexpit(x, sqz = sqz)) < sqz)
#> [1] TRUE
all(abs(p - fexpit(flogit(p, sqz = sqz), sqz = sqz)) < sqz)
#> [1] TRUE
```
