# Get the chromosomes from an object

Get the chromosomes from an object

## Usage

``` r
getChrs(obj)
```

## Arguments

- obj:

  Input SummarizedExperiment object

## Value

A character vector of chromosomes present in an object

## Examples

``` r
data("k562_scrna_chr14", package = "compartmap")
getChrs(k562_scrna_chr14)
#> [1] "chr14"
```
