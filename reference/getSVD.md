# Compute the SVD of a matrix using irlba

Compute the SVD of a matrix using irlba

## Usage

``` r
getSVD(matrix, k = 1, sing.vec = c("left", "right"))
```

## Arguments

- matrix:

  A p x n input matrix

- k:

  Number of singular vectors to return

- sing.vec:

  Whether to return the right or left singular vector

## Value

A singular vector or matrix with sign corresponding to positive values

## Examples

``` r
dummy <- matrix(rnorm(10000), ncol=25)
sing_vec <- getSVD(dummy, k = 1, sing.vec = "right")
```
