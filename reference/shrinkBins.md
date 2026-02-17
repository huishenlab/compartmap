# Employ an eBayes shrinkage approach for bin-level estimates for A/B inference

`shrinkBins` returns shrunken bin-level estimates

## Usage

``` r
shrinkBins(
  x,
  original.x,
  prior.means = NULL,
  chr = NULL,
  res = 1000000,
  targets = NULL,
  jse = TRUE,
  assay = c("rna", "atac", "array"),
  genome = c("hg19", "hg38", "mm9", "mm10")
)
```

## Arguments

- x:

  Input SummarizedExperiment object

- original.x:

  Full sample set SummarizedExperiment object

- prior.means:

  The means of the bin-level prior distribution

- chr:

  The chromosome to operate on

- res:

  Resolution to perform the binning

- targets:

  The column/sample/cell names to shrink towards

- jse:

  Whether to use a James-Stein estimator (default is TRUE)

- assay:

  What assay type this is ("rna", "atac", "array")

- genome:

  What genome are we working with ("hg19", "hg38", "mm9", "mm10")

## Value

A list object to pass to getCorMatrix

## Details

This function computes shrunken bin-level estimates using a James-Stein
estimator (JSE), reformulated as an eBayes procedure. JSE can be used
only if at least 4 targets are provided - any less and `shrinkBins` will
fall back to using Bayes rule which will probably not be great but it
won't explode and may provide some reasonable results anyway

## Examples

``` r
data("k562_scrna_chr14", package = "compartmap")
shrunken.bin.scrna <- shrinkBins(
  x = k562_scrna_chr14,
  original.x = k562_scrna_chr14,
  chr = "chr14", assay = "rna"
)
#> 108 bins created...
```
