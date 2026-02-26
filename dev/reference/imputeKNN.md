# Impute missing values/NAs with KNN

Impute missing values/NAs with KNN

## Usage

``` r
imputeKNN(
  obj,
  rowmax = 0.5,
  colmax = 0.8,
  k = 10,
  maxp = 1500,
  in.place = TRUE,
  drop.sparse.samps = TRUE,
  assay = c("array", "atac", "bisulfite")
)
```

## Arguments

- obj:

  Input SummarizedExperiment object

- rowmax:

  Maximum fraction of NAs that can exist in a row

- colmax:

  Maximum fraction of NAs that can exist in a column/sample

- k:

  Number of neighbors to be used in the imputation

- maxp:

  Largest block of regions/loci imputed using KNN

- in.place:

  Whether to modify the Beta/counts in place (default: TRUE)

- drop.sparse.samps:

  Whether to drop samples that are too sparse (default: TRUE)

- assay:

  The type of assay ("array", "bisulfite")

## Value

Imputed data matrix that is added to the assays slot

## Examples

``` r
if (requireNamespace("minfi", quietly = TRUE)) {
  data("array_data_chr14", package = "compartmap")
  #impute
  imputed <- imputeKNN(array.data.chr14, assay = "array")
}
#> Dropping samples with >80% NAs.
#> Imputing missing data with kNN.
#> Cluster size 12972 broken into 7075 5897 
#> Cluster size 7075 broken into 5100 1975 
#> Cluster size 5100 broken into 1885 3215 
#> Cluster size 1885 broken into 687 1198 
#> Done cluster 687 
#> Done cluster 1198 
#> Done cluster 1885 
#> Cluster size 3215 broken into 1144 2071 
#> Done cluster 1144 
#> Cluster size 2071 broken into 1250 821 
#> Done cluster 1250 
#> Done cluster 821 
#> Done cluster 2071 
#> Done cluster 3215 
#> Done cluster 5100 
#> Cluster size 1975 broken into 1292 683 
#> Done cluster 1292 
#> Done cluster 683 
#> Done cluster 1975 
#> Done cluster 7075 
#> Cluster size 5897 broken into 2293 3604 
#> Cluster size 2293 broken into 995 1298 
#> Done cluster 995 
#> Done cluster 1298 
#> Done cluster 2293 
#> Cluster size 3604 broken into 2213 1391 
#> Cluster size 2213 broken into 1116 1097 
#> Done cluster 1116 
#> Done cluster 1097 
#> Done cluster 2213 
#> Done cluster 1391 
#> Done cluster 3604 
#> Done cluster 5897 
```
