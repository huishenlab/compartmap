# Condense a RaggedExperiment to a list of SummarizedExperiments

Condense a RaggedExperiment to a list of SummarizedExperiments

## Usage

``` r
condenseRE(obj)
```

## Arguments

- obj:

  Input RaggedExperiment

## Value

A list of SummarizedExperiments corresponding to the assays in the input

## Examples

``` r
grl <- GRangesList(
  GRanges(c("A:1-5", "A:4-6", "A:10-15"), score = 1:3),
  GRanges(c("A:1-5", "B:1-3"), score = 4:5)
)
names(grl) <- c("A", "B")
x <- RaggedExperiment(grl)
x.condense <- condenseRE(x)
```
