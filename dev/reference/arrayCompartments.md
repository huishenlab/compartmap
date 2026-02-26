# Estimate A/B compartments from methylation array data

`arrayCompartments` returns estimated A/B compartments from methylation
array data.

## Usage

``` r
arrayCompartments(
  obj,
  res = 1000000,
  chr = NULL,
  group = FALSE,
  targets = NULL,
  bootstrap = TRUE,
  num.bootstraps = 1000,
  preprocess = TRUE,
  array.type = c("hm450", "EPIC"),
  genome = c("hg19", "hg38", "mm9", "mm10"),
  other = NULL,
  boot.parallel = TRUE,
  BPPARAM = bpparam()
)
```

## Arguments

- obj:

  Input SummarizedExperiment object

- res:

  Compartment resolution in bp

- chr:

  What chromosome to work on (leave as NULL to run on all chromosomes)

- group:

  Whether to treat this as a group set of samples

- targets:

  Samples/cells to shrink towards

- bootstrap:

  Whether we should perform bootstrapping of inferred compartments

- num.bootstraps:

  How many bootstraps to run

- preprocess:

  Whether to preprocess the arrays prior to compartment inference

- array.type:

  What type of array is this ("hm450", "EPIC")

- genome:

  What genome to work on ("hg19", "hg38", "mm9", "mm10")

- other:

  Another arbitrary genome to compute compartments on

- boot.parallel:

  Whether to run the bootstrapping in parallel. See details.

- BPPARAM:

  BiocParallelParam object to use for parallelization. See details.

## Value

A RaggedExperiment of inferred compartments

## Details

compartmap uses `BiocParallel` to parallelize operations in four
configurations. The default setting is to parallelize across columns but
not bootstraps using the thread count as reported by
[`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html),
which is usually two cores fewer than the number of available cores.
Parallel bootstrapping is disabled by default to avoid nested
parallelism issues but can be done independent of column-wise
parallelization.

### Available configurations

#### Serial bootstrapping

- Serially with just one core: `BPPARAM = BiocParallel::SerialParam()`

- Parallel across columns and serially across bootstraps:
  `BPPARAM = BiocParallel::MulticoreParam(n)` where `n` is the number of
  threads to use

See
[`?BiocParallel::BiocParallelParam`](https://rdrr.io/pkg/BiocParallel/man/BiocParallelParam-class.html)
for other parallel backends. Parallel backends may also be passed to
[`BiocParallel::register()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
to make them available to `bpparam()`.

#### Parallel bootstrapping

Set `boot.parallel = TRUE` for one the these configurations:

- Serially across columns and parallel across bootstraps: Set \`BPPARAM
  = list(SerialParam(), MulticoreParam(n))'

- Parallel across both columns and bootstraps: Set
  `BPPARAM = list(MulticoreParam(outer), MulticoreParam(inner))` where
  `outer` is the thread count for column-wise operations and `inner` the
  thread count for bootstrapping. The required number of threads is
  given by

`( outer * inner ) + outer`

which is more easily calculated as `outer * (inner + 1)`.

We recommend using an explicit list of two BiocParallelParam backends
over relying on `register()` and `bpparam()` for parallelizing across
bootstraps. With nested `bplapply` calls, the registered backend is used
for both the outer and inner parallel loops. On a system with 8
available threads if the registered backend asks for 4 workers, it will
try to use 20 threads in the nested loops. Instead to use all 8 cores,
set `BPPARAM = list(MulticoreParam(2), MulticoreParam(3))`.

#### Load balancing

Unless you have only 1 chromosome or are not bootstrapping/not
bootstrapping in parallel, you can use nested parallelism. If you are
working on just 1 chromosome, put all cores into the inner bootstrapping
backend. Conversely with multiple chromosomes without bootstrapping, put
all available workers in the outer loop.

In general, use more 'outer' workers, which loop over chromosomes when
`group = TRUE` and cells when `group = FALSE`, than 'inner' workers that
loop over bootstraps. Using 8 outer and 7 inner workers is faster than 7
outer and 8 inner.

When `group = FALSE`, use `MulticoreParam()` only on the outer workers.
We find that parallelizing at both column and bootstrap levels with the
single-cell inference is slower than only parallelizing at the
column-level.

With `group = TRUE`, minimize the difference between the two worker
counts: with 64 total cores, doing 8 outer and 7 inner is faster than 16
outer and 3 inner.

## Examples

``` r
if (requireNamespace("minfi", quietly = TRUE)) {
  data("array_data_chr14", package = "compartmap")
  array_compartments <- arrayCompartments(
    array.data.chr14,
    chr="chr14",
    group=TRUE,
    bootstrap=FALSE,
    genome="hg19",
    array.type="hm450",
    BPPARAM = BiocParallel::SerialParam()
  )
}
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
#> Filtering to open sea CpG loci...
#> Dropping samples with >80% NAs.
#> Imputing missing data with kNN.
#> Cluster size 3332 broken into 518 2814 
#> Done cluster 518 
#> Cluster size 2814 broken into 969 1845 
#> Done cluster 969 
#> Cluster size 1845 broken into 600 1245 
#> Done cluster 600 
#> Done cluster 1245 
#> Done cluster 1845 
#> Done cluster 2814 
#> INFO [2026-02-26 16:11:58] Computing group level compartments
#> INFO [2026-02-26 16:11:58] 
```
