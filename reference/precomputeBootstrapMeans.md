# Pre-compute the global means for bootstrapping compartments

Pre-compute the global means for bootstrapping compartments

## Usage

``` r
precomputeBootstrapMeans(
  obj,
  targets = NULL,
  num.bootstraps = 100,
  assay = c("atac", "rna", "array"),
  parallel = FALSE,
  num.cores = 1
)
```

## Arguments

- obj:

  Input SummarizedExperiment object

- targets:

  Optional targets to shrink towards

- num.bootstraps:

  The number of bootstraps to compute

- assay:

  What type of assay the data are from

- parallel:

  Whether to run in parallel

- num.cores:

  How many cores to use for parallel processing

## Value

A matrix of bootstrapped global means

## Examples

``` r
data("k562_scrna_chr14", package = "compartmap")
scrna.bootstrap.global.means <- precomputeBootstrapMeans(
  k562_scrna_chr14,
  assay = "rna",
  num.bootstraps = 2
)
#> Working on bootstrap 1
#> Working on bootstrap 2
```
