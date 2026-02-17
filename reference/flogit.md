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
#> 0.0002739 0.2348171 0.4932884 0.4984524 0.7651745 0.9994931 

sqz <- 1 / (10**6)
x <- flogit(p, sqz = sqz)
summary(x)
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -8.20172 -1.18131 -0.02685 -0.03337  1.18126  7.58617 

all(abs(p - fexpit(x, sqz = sqz)) < sqz)
#> [1] TRUE
all(abs(p - fexpit(flogit(p, sqz = sqz), sqz = sqz)) < sqz)
#> [1] TRUE
```
