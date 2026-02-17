# Remove columns/cells/samples with NAs exceeding a threshold. See `cleanAssay()`

Remove columns/cells/samples with NAs exceeding a threshold. See
[`cleanAssay()`](https://huishenlab.github.io/compartmap/reference/cleanAssay.md)

## Usage

``` r
cleanAssayCols(se, na.max = 0.8, assay = c("array", "bisulfite"))
```

## Arguments

- se:

  Input SummarizedExperiment object

- na.max:

  The maximum number of NAs allowed as a fraction

- assay:

  The type of assay we are working with

## Value

A filtered matrix

## Examples

``` r
if (requireNamespace("minfi", quietly = TRUE)) {
  data("array_data_chr14", package = "compartmap")
  compartmap:::cleanAssayCols(array.data.chr14, assay = "array")
}
#> class: GenomicRatioSet 
#> dim: 15078 11 
#> metadata(1): SNPs
#> assays(2): Beta CN
#> rownames(15078): cg07948143 cg20517615 ... cg23902114 cg24746738
#> rowData names(0):
#> colnames(11): naive.1 rTreg.2 ... act_rTreg.10 birth.11
#> colData names(10): Sample_Name Sample_Well ... Basename filenames
#> Annotation
#>   array: IlluminaHumanMethylation450k
#>   annotation: ilmn12.hg19
#> Preprocessing
#>   Method: SeSAMe (type I)
#>   minfi version: 1.27.8
#>   Manifest version: 0.6.0
```
