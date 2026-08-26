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
#>     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#> -3.38866 -0.72615 -0.03265 -0.04192  0.62820  2.97110 

sqz <- 1 / (10**6)
p <- fexpit(x, sqz = sqz)
summary(p)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.03265 0.32604 0.49184 0.49080 0.65208 0.95125 

all((abs(x - flogit(p)) / x) < sqz)
#> [1] TRUE
all(abs(x - flogit(fexpit(x))) < sqz)
#> [1] TRUE
```
