# Gather open sea CpG from a GRanges of CpG islands

This function accepts a GRanges input of CpG islands that can be derived
from UCSC table browser and rtracklayer::import(yourbed.bed). It resizes
the intevals to create 4kb flanking regions around CpG islands. Open sea
regions, as defined by [Fortin by and Hansen (Genome Biology,
2015)](https://doi.org/10.1186/s13059-015-0741-y), are outside this
flanking regions and obtained as their complement.

## Usage

``` r
getOpenSeas(gr)
```

## Arguments

- gr:

  Input GRanges of CpG islands

## Value

GRanges object that can be used with filterOpenSea()

## Examples

``` r
#cpgi <- rtracklayer::import(system.file("inst/extdata/mm10_cpgi.bed", package = "compartmap"))
#opensea_cpg <- getOpenSeas(cpgi)
```
