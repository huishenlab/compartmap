# Inverse Fisher's Z transformation

`ifisherZ` returns the inverse (squeezed) Fisher's Z transformed
Pearson's r.

## Usage

``` r
ifisherZ(zmat)
```

## Arguments

- zmat:

  matrix of Fisher's Z transformed Pearson correlations or an
  eignevector

## Value

Back transformed Fisher's Z Pearson correlations

## Examples

``` r
# Generate a random binary (-1, 1) matrix
(mat <- matrix(sample(c(1,-1), 25, replace = TRUE), ncol = 5))
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   -1   -1   -1    1   -1
#> [2,]    1    1    1    1    1
#> [3,]    1    1    1   -1   -1
#> [4,]   -1    1   -1   -1    1
#> [5,]    1   -1   -1   -1   -1

# Correct matrix diag
diag(mat) <- 1

# Transform
(mat.transform <- fisherZ(mat))
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,]  7.254329 -7.254329 -7.254329  7.254329 -7.254329
#> [2,]  7.254329  7.254329  7.254329  7.254329  7.254329
#> [3,]  7.254329  7.254329  7.254329 -7.254329 -7.254329
#> [4,] -7.254329  7.254329 -7.254329  7.254329  7.254329
#> [5,]  7.254329 -7.254329 -7.254329 -7.254329  7.254329

#Back transform
ifisherZ(mat.transform)
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,]  0.999999 -0.999999 -0.999999  0.999999 -0.999999
#> [2,]  0.999999  0.999999  0.999999  0.999999  0.999999
#> [3,]  0.999999  0.999999  0.999999 -0.999999 -0.999999
#> [4,] -0.999999  0.999999 -0.999999  0.999999  0.999999
#> [5,]  0.999999 -0.999999 -0.999999 -0.999999  0.999999
```
