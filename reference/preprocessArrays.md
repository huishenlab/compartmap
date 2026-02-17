# Preprocess arrays for compartment inference

Preprocess arrays for compartment inference

## Usage

``` r
preprocessArrays(
  obj,
  genome = c("hg19", "hg38", "mm9", "mm10"),
  other = NULL,
  array.type = c("hm450", "EPIC")
)
```

## Arguments

- obj:

  Input SummarizedExperiment

- genome:

  What genome are we working on ("hg19", "hg38", "mm9", "mm10")

- other:

  Another arbitrary genome to compute compartments on

- array.type:

  What type of array is this ("hm450", "EPIC")

## Value

A preprocessed SummarizedExperiment to compute compartments

## Examples

``` r
if (requireNamespace("minfiData", quietly = TRUE)) {
  grSet <- minfi::preprocessNoob(minfiData::RGsetEx.sub) |>
    minfi::ratioConvert() |>
    minfi::mapToGenome()
  preprocessArrays(grSet)
}
#> Loading required package: IlluminaHumanMethylation450kmanifest
#> Loading required package: minfi
#> Loading required package: Biostrings
#> Loading required package: XVector
#> 
#> Attaching package: ‘Biostrings’
#> The following object is masked from ‘package:base’:
#> 
#>     strsplit
#> Loading required package: bumphunter
#> Loading required package: foreach
#> Loading required package: iterators
#> Loading required package: parallel
#> Loading required package: locfit
#> locfit 1.5-9.12   2025-03-05
#> Loading required package: IlluminaHumanMethylation450kanno.ilmn12.hg19
#> Filtering to open sea CpG loci...
#> Converting to squeezed M-values.
#> class: GenomicRatioSet 
#> dim: 47 6 
#> metadata(0):
#> assays(3): Beta M CN
#> rownames(47): cg01003813 cg01051089 ... cg01757887 cg03930849
#> rowData names(0):
#> colnames(6): 5723646052_R02C02 5723646052_R04C01 ... 5723646053_R05C02
#>   5723646053_R06C02
#> colData names(13): Sample_Name Sample_Well ... Basename filenames
#> Annotation
#>   array: IlluminaHumanMethylation450k
#>   annotation: ilmn12.hg19
#> Preprocessing
#>   Method: NA
#>   minfi version: NA
#>   Manifest version: NA
```
