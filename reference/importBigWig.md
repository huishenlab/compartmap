# Import and optionally summarize a bigwig at a given resolution

Import and optionally summarize a bigwig at a given resolution

## Usage

``` r
importBigWig(
  bw,
  bins = NULL,
  summarize = FALSE,
  genome = c("hg19", "hg38", "mm9", "mm10")
)
```

## Arguments

- bw:

  Path a bigwig file

- bins:

  Optional set of bins as a GRanges to summarize the bigwig to

- summarize:

  Whether to perform mean summarization

- genome:

  Which genome is the bigwig from ("hg19", "hg38", "mm9", "mm10")

## Value

SummarizedExperiment object with rowRanges corresponding to summarized
features
