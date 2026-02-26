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
#>       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#> -3.5839083 -0.7056768  0.0002002 -0.0274849  0.5874546  3.2738377 

sqz <- 1 / (10**6)
p <- fexpit(x, sqz = sqz)
summary(p)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.02702 0.33055 0.50005 0.49406 0.64278 0.96352 

all((abs(x - flogit(p)) / x) < sqz)
#> [1] TRUE
all(abs(x - flogit(fexpit(x))) < sqz)
#> [1] TRUE
```
