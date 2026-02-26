# Transform/normalize compartment calls using TF-IDF on HDF5-backed objects

Transform/normalize compartment calls using TF-IDF on HDF5-backed
objects

## Usage

``` r
hdf5TFIDF(h5, scale.factor = 100000, return.dense = FALSE, return.se = FALSE)
```

## Arguments

- h5:

  SummarizedExperiment object, DelayedMatrix, or a normal matrix

- scale.factor:

  Scaling factor for the term-frequency (TF)

- return.dense:

  Whether to return a dense, in memory matrix

- return.se:

  Whether to return the TF-IDF matrix as a new assay in the
  SummarizedExperiment

## Value

A TF-IDF transformed matrix of the same dimensions as the input

## Examples

``` r
m <- 1000
n <- 100
mat <- round(matrix(runif(m * n), m, n))
# Input needs to be a tall matrix
tfidf <- hdf5TFIDF(mat)
#> Computing term frequency.
#> Computing inverse document frequency.
#> TF-IDF
```
