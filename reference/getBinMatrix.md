# Generate bins for A/B compartment estimation

Generate bins across a user defined chromosome for A/B compartment
estimation. A/B compartment estimation can be used for non-supported
genomes if chr.end is set.

## Usage

``` r
getBinMatrix(
  mat,
  genloc,
  chr = "chr1",
  chr.start = 0,
  chr.end = NULL,
  res = 100000,
  FUN = sum,
  genome = c("hg19", "hg38", "mm9", "mm10")
)
```

## Arguments

- mat:

  A p x n matrix where p (rows) = loci and n (columns) = samples/cells

- genloc:

  GRanges object that contains corresponding genomic locations of the
  loci

- chr:

  Chromosome to be analyzed

- chr.start:

  Starting position (in bp) to be analyzed

- chr.end:

  End position (in bp) to be analyzed

- res:

  Binning resolution (in bp)

- FUN:

  Function to be used to summarize information within a bin

- genome:

  Genome corresponding to the input data ("hg19", "hg38", "mm9", "mm10")

## Value

A list object to pass to getCorMatrix

## Details

This function is used to generate a list object to be passed to
getCorMatrix

## Examples

``` r
library(GenomicRanges)

# Generate random genomic intervals of 1-1000 bp on chr1-22
# Modified from https://www.biostars.org/p/225520/
genome.gr <- getGenome("hg19")
random_genomic_int <- data.frame(chr = rep("chr14", 100))
random_genomic_int$start <- apply(random_genomic_int, 1, function(x) {
  round(runif(1, 0, getSeqLengths(genome.gr, chr = x)[[1]]), 0)
})
random_genomic_int$end <- random_genomic_int$start + runif(1, 1, 1000)
random_genomic_int$strand <- "*"

# Generate random counts
counts <- rnbinom(1000, 1.2, 0.4)

# Build random counts for 10 samples
count.mat <- matrix(sample(counts, nrow(random_genomic_int) * 10, replace = FALSE), ncol = 10)
colnames(count.mat) <- paste0("sample_", seq(1:10))

# Bin counts
bin.counts <- getBinMatrix(
  count.mat,
  makeGRangesFromDataFrame(random_genomic_int),
  chr = "chr14",
  genome = "hg19"
)
#> 1074 bins created...
```
