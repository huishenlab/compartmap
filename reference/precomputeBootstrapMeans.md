# Pre-compute the global means for bootstrapping compartments

Pre-compute the global means for bootstrapping compartments

## Usage

``` r
precomputeBootstrapMeans(
  obj,
  BPPARAM,
  targets = NULL,
  num.bootstraps = 100,
  assay = c("atac", "rna", "array")
)
```

## Arguments

- obj:

  Input SummarizedExperiment object

- BPPARAM:

  BiocParallelParam for parallelizing computation

- targets:

  Optional targets to shrink towards

- num.bootstraps:

  The number of bootstraps to compute

- assay:

  What type of assay the data are from

## Value

A matrix of bootstrapped global means

## Examples

``` r
data("k562_scrna_chr14", package = "compartmap")
scrna.bootstrap.global.means <- precomputeBootstrapMeans(
  k562_scrna_chr14,
  BPPARAM = BiocParallel::SerialParam(),
  assay = "rna",
  num.bootstraps = 2
)
```
