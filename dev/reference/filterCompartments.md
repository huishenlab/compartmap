# Filter compartments using confidence estimates and eigenvalue thresholds

Filter compartments using confidence estimates and eigenvalue thresholds

## Usage

``` r
filterCompartments(obj, min.conf = 0.7, min.eigen = 0.02)
```

## Arguments

- obj:

  Output of condenseSE or fixCompartments

- min.conf:

  Minimum confidence estimate to use when filtering

- min.eigen:

  Minimum absolute eigenvalue to use when filtering

## Value

A filtered/subset of the input object/list
