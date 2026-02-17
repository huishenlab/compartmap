# Convert a sparse matrix to a dense matrix in a block-wise fashion

Convert a sparse matrix to a dense matrix in a block-wise fashion

## Usage

``` r
sparseToDenseMatrix(
  mat,
  blockwise = TRUE,
  chunk.by = "row",
  chunk.size = 100000,
  parallel = FALSE,
  cores = 2
)
```

## Arguments

- mat:

  Input sparse matrix

- blockwise:

  Whether to do the coercion in a block-wise manner

- chunk.by:

  Whether to chunk by `"row"` or `"column"`

- chunk.size:

  The size of the chunks to use for coercion

- parallel:

  Whether to perform the coercion in parallel

- cores:

  The number of cores to use in the parallel coercion

## Value

A dense matrix of the same dimensions as the input

## Examples

``` r
# make a sparse binary matrix
library(Matrix)
m <- 100
n <- 1000
mat <- round(matrix(runif(m * n), m, n))
mat.sparse <- Matrix(mat, sparse = TRUE)

# coerce back
mat.dense <- sparseToDenseMatrix(mat.sparse, chunk.size = 10)
#> Breaking into row chunks.

# make sure they are the same dimensions
dim(mat) == dim(mat.dense)
#> [1] TRUE TRUE

# make sure they are the same numerically
all(mat == mat.dense)
#> [1] TRUE
```
