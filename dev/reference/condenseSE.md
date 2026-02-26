# Condense the output of condenseRE to reconstruct per-sample GRanges objects to plot

Condense the output of condenseRE to reconstruct per-sample GRanges
objects to plot

## Usage

``` r
condenseSE(obj, sample.name = NULL)
```

## Arguments

- obj:

  Output of condenseRE or can be a RaggedExperiment

- sample.name:

  Vector of samples/cells to extract

## Value

GRanges or list of per-sample GRanges to pass to plotAB or export

## Examples

``` r
grl <- GRangesList(
  GRanges(c("A:1-5", "A:4-6", "A:10-15"), score = 1:3),
  GRanges(c("A:1-5", "B:1-3"), score = 4:5)
)
names(grl) <- c("A", "B")
x <- RaggedExperiment(grl)
condense.x <- condenseSE(x, sample.name = "A")
```
