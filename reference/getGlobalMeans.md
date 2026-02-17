# Get the global means of a matrix

Get the global means of a matrix

## Usage

``` r
getGlobalMeans(obj, targets = NULL, assay = c("atac", "rna", "array"))
```

## Arguments

- obj:

  Input SummarizedExperiment object

- targets:

  Column names or indices to indicate which samples to shrink towards

- assay:

  What type of assay the data are from

## Value

A vector of global or targeted means

## Examples

``` r
data("k562_scrna_chr14", package = "compartmap")
scrna.global.means <- getGlobalMeans(k562_scrna_chr14, assay = "rna")
```
