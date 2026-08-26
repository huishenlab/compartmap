# Find overlaps between `CompartmentCall` objects

Find overlaps between `CompartmentCall` objects

## Usage

``` r
## S7 method for classes <compartmap::CompartmentCall>, <compartmap::CompartmentCall>
findOverlaps(
  query,
  subject,
  maxgap = -1L,
  minoverlap = 0L,
  type = c("any", "start", "end", "within", "equal"),
  select = c("all", "first", "last", "arbitrary"),
  ignore.strand = FALSE
)
```

## Arguments

- query:

  A string vector of chromosomes to subset to

- subject:

  A `CompartmentCall` object
