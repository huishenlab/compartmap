# Filter to open sea CpG loci

Filter to open sea CpG loci

## Usage

``` r
filterOpenSea(obj, genome = c("hg19", "hg38", "mm10", "mm9"), other = NULL)
```

## Arguments

- obj:

  Input SummarizedExperiment or GRanges object

- genome:

  Which genome to filter

- other:

  GRanges of open sea regions (TODO)

## Value

Filtered to open sea CpG loci

## Examples

``` r
if (requireNamespace("minfi", quietly = TRUE)) {
  data("array_data_chr14", package = "compartmap")
  opensea <- filterOpenSea(array.data.chr14, genome = "hg19")
}
#> Filtering to open sea CpG loci...
```
