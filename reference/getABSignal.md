# Calculate Pearson correlations of smoothed eigenvectors

This function is used to generate a list x to be passed to getABSignal

## Usage

``` r
getABSignal(x, squeeze = FALSE, assay = c("rna", "atac", "array"))
```

## Arguments

- x:

  A list object from getCorMatrix

- squeeze:

  Whether squeezing was used (implies Fisher's Z transformation)

- assay:

  What kind of assay are we working on ("array", "atac", "array")

## Value

A list x to pass to getABSignal

## Examples

``` r
library(SummarizedExperiment)
library(BiocSingular)

#Generate random genomic intervals of 1-1000 bp on chr1-22
#Modified from https://www.biostars.org/p/225520/
random_genomic_int <- data.frame(chr = rep("chr14", 100))
random_genomic_int$start <- apply(random_genomic_int, 1, function(x) {
  round(runif(1, 0, getSeqLengths(getGenome("hg19"), chr = x)[[1]]), 0)
})
random_genomic_int$end <- random_genomic_int$start + runif(1, 1, 1000)
random_genomic_int$strand <- "*"

#Generate random counts
counts <- rnbinom(1000, 1.2, 0.4)

#Build random counts for 10 samples
count.mat <- matrix(sample(counts, nrow(random_genomic_int) * 10, replace = FALSE), ncol = 10)
colnames(count.mat) <- paste0("sample_", seq(1:10))

#Bin counts
bin.counts <- getBinMatrix(
  count.mat,
  makeGRangesFromDataFrame(random_genomic_int),
  chr = "chr14",
  genome = "hg19"
)
#> 1074 bins created...

#Calculate correlations
bin.cor.counts <- getCorMatrix(bin.counts)
#> Calculating correlations...
#> Done...

#Get A/B signal
absignal <- getABSignal(bin.cor.counts)
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
```
