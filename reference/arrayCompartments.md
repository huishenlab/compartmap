# Estimate A/B compartments from methylation array data

`arrayCompartments` returns estimated A/B compartments from methylation
array data.

## Usage

``` r
arrayCompartments(
  obj,
  res = 1000000,
  chr = NULL,
  targets = NULL,
  preprocess = TRUE,
  parallel = TRUE,
  cores = 2,
  bootstrap = TRUE,
  num.bootstraps = 1000,
  boot.parallel = TRUE,
  boot.cores = 2,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  group = FALSE,
  other = NULL,
  array.type = c("hm450", "EPIC")
)
```

## Arguments

- obj:

  Input SummarizedExperiment object

- res:

  Compartment resolution in bp

- chr:

  What chromosome to work on (leave as NULL to run on all chromosomes)

- targets:

  Samples/cells to shrink towards

- preprocess:

  Whether to preprocess the arrays prior to compartment inference

- parallel:

  Whether to run samples in parallel

- cores:

  How many cores to use when running samples in parallel

- bootstrap:

  Whether we should perform bootstrapping of inferred compartments

- num.bootstraps:

  How many bootstraps to run

- boot.parallel:

  Whether to run the bootstrapping in parallel

- boot.cores:

  How many cores to use for the bootstrapping

- genome:

  What genome to work on ("hg19", "hg38", "mm9", "mm10")

- group:

  Whether to treat this as a group set of samples

- other:

  Another arbitrary genome to compute compartments on

- array.type:

  What type of array is this ("hm450", "EPIC")

## Value

A RaggedExperiment of inferred compartments

## Examples

``` r
if (requireNamespace("minfi", quietly = TRUE)) {
  data("array_data_chr14", package = "compartmap")
  array_compartments <- arrayCompartments(
    array.data.chr14,
    chr="chr14",
    parallel=FALSE,
    bootstrap=FALSE,
    genome="hg19",
    array.type="hm450"
  )
}
#> Setting options('download.file.method.GEOquery'='auto')
#> Setting options('GEOquery.inmemory.gpl'=FALSE)
#> Filtering to open sea CpG loci...
#> Converting to squeezed M-values.
#> Imputing missing values.
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
#> Working on naive.1
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
#> Working on rTreg.2
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
#> Working on act_naive.3
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
#> Working on naive.4
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
#> Working on act_naive.5
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
#> Working on act_rTreg.6
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
#> Working on naive.7
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
#> Working on rTreg.8
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
#> Working on act_naive.9
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
#> Working on act_rTreg.10
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
#> Working on birth.11
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
```
