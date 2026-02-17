# Wrapper to denoise a correlation matrix using a Random Matrix Theory approach

Wrapper to denoise a correlation matrix using a Random Matrix Theory
approach

## Usage

``` r
getDenoisedCorMatrix(
  obj,
  res = 1000000,
  chr = "chr14",
  genome = c("hg19", "hg38", "mm9", "mm10"),
  iter = 2,
  targets = NULL,
  prior.means = NULL,
  assay = c("rna", "atac", "array")
)
```

## Arguments

- obj:

  SummarizedExperiment object with rowRanges for each feature and
  colnames

- res:

  The resolution desired (default is a megabase 1e6)

- chr:

  Which chromosome to perform the denoising

- genome:

  Which genome (default is hg19)

- iter:

  How many iterations to perform denoising

- targets:

  Samples/cells to shrink towards

- prior.means:

  The means of the bin-level prior distribution (default will compute
  them for you)

- assay:

  What assay type this is ("rna", "atac")

## Value

A denoised correlation matrix object for plotting with plotCorMatrix

## Examples

``` r
data("k562_scrna_chr14", package = "compartmap")
denoised_cor_mat <- getDenoisedCorMatrix(k562_scrna_chr14, genome = "hg19", assay = "rna")
#> Shrinking bins with the JSE.
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Denoising the correlation matrix using RMT.
#> Iterative denoising. Iteration: 2
```
