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
#> -3.583908 -0.705677 -0.006277 -0.029761  0.578246  3.273838 

sqz <- 1 / (10**6)
p <- fexpit(x, sqz = sqz)
summary(p)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.02702 0.33055 0.49843 0.49349 0.64066 0.96352 

all((abs(x - flogit(p)) / x) < sqz)
#> [1] TRUE
all(abs(x - flogit(fexpit(x))) < sqz)
#> [1] TRUE
```
