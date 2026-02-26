# Helper to invert, or "fix", compartments that have a minimum confidence score (1-min.conf)

Helper to invert, or "fix", compartments that have a minimum confidence
score (1-min.conf)

## Usage

``` r
flipper(input_obj, min.conf)
```

## Arguments

- input_obj:

  Input RaggedExperiment or output of condenseSE

- min.conf:

  Minimum confidence score to use

## Value

A "fixed" set of compartments
