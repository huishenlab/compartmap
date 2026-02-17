# A wrapper function to generate a GRanges object of chromatin domain inflection points

A wrapper function to generate a GRanges object of chromatin domain
inflection points

## Usage

``` r
getDomainInflections(
  gr,
  what = "score",
  res = 1000000,
  chrs = c(paste0("chr", 1:22), "chrX"),
  genome = c("hg19", "hg38", "mm9", "mm10")
)
```

## Arguments

- gr:

  Input GRanges object with mcols column corresponding to chromatin
  domains

- what:

  The name of the column containing the chromatin domain information

- res:

  What resolution the domains were called

- chrs:

  Which chromosomes to work on

- genome:

  Which genome does the input data come from

## Value

A GRanges object of compartment inflection points

## Examples

``` r
data("k562_scrna_chr14", package = "compartmap")
chr14_domains <- scCompartments(k562_scrna_chr14,
  res = 1e6, genome = "hg19",
  group = TRUE, bootstrap = FALSE
)
#> Assuming we want to process all chromosomes.
#> Computing compartments for chr14
#> 108 bins created...
#> Calculating correlations...
#> Done...
#> Calculating eigenvectors.
#> Smoothing eigenvector.
#> Done smoothing.
chr14_domain_inflections <- getDomainInflections(chr14_domains, what = "pc")
#> Tiling genome.
#> Warning: GRanges object contains 23 out-of-bound ranges located on sequences chr1,
#>   chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13,
#>   chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, and chrX. Note
#>   that ranges located on a sequence whose length is unknown (NA) or on a
#>   circular sequence are not considered out-of-bound (use seqlengths() and
#>   isCircular() to get the lengths and circularity flags of the underlying
#>   sequences). You can use trim() to trim these ranges. See
#>   ?`trim,GenomicRanges-method` for more information.
#> Warning: GRanges object contains 23 out-of-bound ranges located on sequences chr1,
#>   chr2, chr3, chr4, chr5, chr6, chr7, chr8, chr9, chr10, chr11, chr12, chr13,
#>   chr14, chr15, chr16, chr17, chr18, chr19, chr20, chr21, chr22, and chrX. Note
#>   that ranges located on a sequence whose length is unknown (NA) or on a
#>   circular sequence are not considered out-of-bound (use seqlengths() and
#>   isCircular() to get the lengths and circularity flags of the underlying
#>   sequences). You can use trim() to trim these ranges. See
#>   ?`trim,GenomicRanges-method` for more information.
#> Finding overlaps.
#> Injecting eigenvalues.
#> Contiguous runs. Finding inflections.
```
