# Get the open and closed compartment calls based on sign of singular values

Get the open and closed compartment calls based on sign of singular
values

## Usage

``` r
extractOpenClosed(gr, cutoff = 0, assay = c("rna", "atac", "array"))
```

## Arguments

- gr:

  Input GRanges with associated mcols that represent singular values

- cutoff:

  Threshold to define open and closed states

- assay:

  The type of assay we are working with

## Value

A vector of binary/categorical compartment states

## Examples

``` r
dummy <- matrix(rnorm(10000), ncol = 25)
sing_vec <- getSVD(dummy, k = 1, sing.vec = "right")
```
