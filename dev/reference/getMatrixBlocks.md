# Get chunked sets of row-wise or column-wise indices of a matrix

Get chunked sets of row-wise or column-wise indices of a matrix

## Usage

``` r
getMatrixBlocks(mat, chunk.size = 100000, chunk.by = "row")
```

## Arguments

- mat:

  Input matrix

- chunk.size:

  The size of the chunks to use for coercion

- chunk.by:

  Whether to chunk in a `"row"`- or `"column"`-wise fashion

## Value

A set of chunked indices

## Examples

``` r
# make a sparse binary matrix
library(Matrix)
m <- 100
n <- 1000
mat <- round(matrix(runif(m * n), m, n))
mat.sparse <- Matrix(mat, sparse = TRUE)

# get row-wise chunks of 10
chunks <- getMatrixBlocks(mat.sparse, chunk.size = 10)
#> Breaking into row chunks.
```
