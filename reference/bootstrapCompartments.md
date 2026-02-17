# Non-parametric bootstrapping of compartments and summarization of bootstraps/compute confidence intervals

Non-parametric bootstrapping of compartments and summarization of
bootstraps/compute confidence intervals

## Usage

``` r
bootstrapCompartments(
  obj,
  original.obj,
  bootstrap.samples = 1000,
  chr = "chr14",
  assay = c("rna", "atac", "array"),
  parallel = TRUE,
  cores = 2,
  targets = NULL,
  res = 1000000,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  q = 0.95,
  svd = NULL,
  group = FALSE,
  bootstrap.means = NULL
)
```

## Arguments

- obj:

  List object of computed compartments for a sample with 'pc' and 'gr'
  as elements

- original.obj:

  The original, full input SummarizedExperiment of all samples/cells

- bootstrap.samples:

  How many bootstraps to run

- chr:

  Which chromosome to operate on

- assay:

  What sort of assay are we working on

- parallel:

  Whether to run the bootstrapping in parallel

- cores:

  How many cores to use for parallel processing

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

- group:

  Whether this is for group-level inference

- bootstrap.means:

  Pre-computed bootstrap means matrix

## Value

Compartment estimates with summarized bootstraps and confidence
intervals

## Examples

``` r
# this needs a good example
```
