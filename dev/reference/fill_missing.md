# Fill missing genomic bins in `CompartmentCalls` using a reference GRanges

Compartmap may drop genomic bins with insufficient data and the
resulting GRanges object may not have all the bins of the region it was
run which means using the same region and resolution on different inputs
does not guarantee the same output bins. Having different set of bins
between calls despite the same input regions and resolution prevents the
creation of `MultiCompartmapCall` objcets that expect the same GRanges
across all input `CompartmentCall` objects. `fill_missing()` adds the
missing bins according to a larger reference set of bins, filling
missing data with NA. All bins in `x` must be present in the reference
bins.

## Usage

``` r
fill_missing(x, ref.gr)
```

## Arguments

- x:

  A `CompartmentCall` object

- ref.gr:

  The `GRanges` object to use as the full reference set of regions
