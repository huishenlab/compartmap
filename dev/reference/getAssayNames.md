# Get the assay names from a SummarizedExperiment object

Get the assay names from a SummarizedExperiment object

## Usage

``` r
getAssayNames(se)
```

## Arguments

- se:

  Input SummarizedExperiment object

## Value

The names of the assays

## Examples

``` r
data("k562_scrna_chr14", package = "compartmap")
getAssayNames(k562_scrna_chr14)
#> [1] "counts"
```
