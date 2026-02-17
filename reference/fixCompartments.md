# Invert, or "fix", compartments that have a minimum confidence score (1-min.conf)

Invert, or "fix", compartments that have a minimum confidence score
(1-min.conf)

## Usage

``` r
fixCompartments(obj, min.conf = 0.8, parallel = FALSE, cores = 1)
```

## Arguments

- obj:

  Input RaggedExperiment or output of condenseSE

- min.conf:

  Minimum confidence score to use

- parallel:

  Whether to run in parallel

- cores:

  How many cores to use if running in parallel

## Value

A "fixed" set of compartments
