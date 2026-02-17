# Fisher's Z transformation

`fisherZ` returns (squeezed) Fisher's Z transformed Pearson's r

## Usage

``` r
fisherZ(cormat)
```

## Arguments

- cormat:

  Pearson correlation matrix

## Value

Fisher Z transformed Pearson correlations

## Examples

``` r
#Generate a random binary (-1, 1) matrix
(mat <- matrix(sample(c(1,-1), 25, replace = TRUE), ncol = 5))
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1   -1    1    1    1
#> [2,]    1    1   -1   -1    1
#> [3,]    1    1    1   -1    1
#> [4,]   -1   -1    1    1    1
#> [5,]    1    1    1   -1   -1

#Correct matrix diag
diag(mat) <- 1

#Transform
fisherZ(mat)
#>           [,1]      [,2]      [,3]      [,4]     [,5]
#> [1,]  7.254329 -7.254329  7.254329  7.254329 7.254329
#> [2,]  7.254329  7.254329 -7.254329 -7.254329 7.254329
#> [3,]  7.254329  7.254329  7.254329 -7.254329 7.254329
#> [4,] -7.254329 -7.254329  7.254329  7.254329 7.254329
#> [5,]  7.254329  7.254329  7.254329 -7.254329 7.254329
```
