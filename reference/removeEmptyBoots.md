# Remove bootstrap estimates that failed

This can happen if the correlation between the bins and eigenvector
fails theoretically we can recover these but need an additional utility
to find consensus

## Usage

``` r
removeEmptyBoots(obj)
```

## Arguments

- obj:

  Input list object with elements 'pc' and 'gr'

## Value

A filtered list object
