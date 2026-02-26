# Non-parametric bootstrapping of compartments and summarization of bootstraps/compute confidence intervals

Non-parametric bootstrapping of compartments and summarization of
bootstraps/compute confidence intervals

## Usage

``` r
bootstrapCompartments(
  obj,
  original.obj,
  BPPARAM,
  bootstrap.samples = 1000,
  chr = "chr14",
  group = FALSE,
  assay = c("rna", "atac", "array"),
  targets = NULL,
  res = 1000000,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  q = 0.95,
  svd = NULL,
  bootstrap.means = NULL
)
```

## Arguments

- obj:

  List object of computed compartments for a sample with 'pc' and 'gr'
  as elements

- original.obj:

  The original, full input SummarizedExperiment of all samples/cells

- BPPARAM:

  BiocParallelParam for parallelizing bootstrapping

- bootstrap.samples:

  How many bootstraps to run

- chr:

  Which chromosome to operate on

- group:

  Whether this is for group-level inference

- assay:

  What sort of assay are we working on

- targets:

  Targets to shrink towards

- res:

  The compartment resolution

- genome:

  What genome are we working on

- q:

  What sort of confidence intervals are we computing (e.g. 0.95 for 95
  percentCI)

- svd:

  The original compartment calls as a GRanges object

- bootstrap.means:

  Pre-computed bootstrap means matrix

## Value

Compartment estimates with summarized bootstraps and confidence
intervals

## Examples

``` r
# this needs a good example
```
