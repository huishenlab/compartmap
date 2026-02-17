# Get a GRanges object from bundled compartmap genomes

Get a GRanges object from bundled compartmap genomes

## Usage

``` r
getGenome(genome = c("hg19", "hg38", "mm9", "mm10"), type = "genome")
```

## Arguments

- genome:

  The desired genome to use ("hg19", "hg38", "mm9", "mm10")

- type:

  The type of data - full genome or open sea regions

## Value

Granges of the genome

## Examples

``` r
hg19 <- getGenome(genome = "hg19")
```
