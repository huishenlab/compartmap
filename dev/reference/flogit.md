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
#> 0.0002544 0.2452321 0.5274872 0.5070668 0.7782408 0.9997144 

sqz <- 1 / (10**6)
x <- flogit(p, sqz = sqz)
summary(x)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> -8.2755 -1.1242  0.1101  0.0348  1.2554  8.1596 

all(abs(p - fexpit(x, sqz = sqz)) < sqz)
#> [1] TRUE
all(abs(p - fexpit(flogit(p, sqz = sqz), sqz = sqz)) < sqz)
#> [1] TRUE
```
