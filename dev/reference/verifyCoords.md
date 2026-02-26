# Throw error if assay does not contain coordinates

Throw error if assay does not contain coordinates

## Usage

``` r
verifyCoords(obj)
```

## Arguments

- obj:

  Input object

## Examples

``` r
data("k562_scrna_chr14", package = "compartmap")
compartmap:::verifyCoords(k562_scrna_chr14)
```
