# Summarize the bootstrap compartment estimates and compute Agresti-Coull confidence intervals

Summarize the bootstrap compartment estimates and compute Agresti-Coull
confidence intervals

## Usage

``` r
summarizeBootstraps(boot.list, est.ab, q = 0.95, assay = c("rna", "atac"))
```

## Arguments

- boot.list:

  List of bootstraps as GRanges objects

- est.ab:

  The original compartment calls

- q:

  Confidence interval to compute (0.95 for 95 percent CI)

- assay:

  Type of assay we are working with

## Value

A GRanges object with bootstraps summarized
