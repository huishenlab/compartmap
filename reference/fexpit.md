# Helper function: expanded expit

Helper function: expanded expit

## Usage

``` r
fexpit(x, sqz = 0.000001)
```

## Arguments

- x:

  a vector of values between -Inf and +Inf

- sqz:

  the amount by which we 'squoze', default is .000001

## Value

       a vector of values between 0 and 1 inclusive

## Examples

``` r
x <- rnorm(n = 1000)
summary(x)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> -3.476112 -0.701113 -0.041293 -0.001381  0.697089  3.090699 

sqz <- 1 / (10**6)
p <- fexpit(x, sqz = sqz)
summary(p)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>  0.0300  0.3316  0.4897  0.4989  0.6675  0.9565 

all((abs(x - flogit(p)) / x) < sqz)
#> [1] TRUE
all(abs(x - flogit(fexpit(x))) < sqz)
#> [1] TRUE
```
