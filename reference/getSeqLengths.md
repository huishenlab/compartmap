# Get the seqlengths of a chromosome from a given genome's GRanges

The goal for this function is to eliminate the need to lug around large
packages when we only want seqlengths for things.

## Usage

``` r
getSeqLengths(genome.gr, chr = "chr14")
```

## Arguments

- genome.gr:

  A GRanges object of the genome (from
  [`getGenome()`](https://huishenlab.github.io/compartmap/reference/getGenome.md))

- chr:

  What chromosome to extract the seqlengths of

## Value

The seqlengths of a specific chromosome

## Examples

``` r
hg19.chr14.seqlengths <- getSeqLengths(getGenome('hg19'), chr = "chr14")
```
