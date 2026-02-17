# Generate function to filter rows/columns with NAs exceeding a threshold

Generate function to filter rows/columns with NAs exceeding a threshold

## Usage

``` r
cleanAssay(by = c("row", "col"))
```

## Arguments

- by:

  Whether to filter by rows or columns

## Value

A function to filter assay rows/columns

## Details

Since removing NAs from rows vs columns only differs by whether rowMeans
or colMeans is used, and by where the comma goes in the subset
operation, code repetition can be avoided by consolidating these
operations. This `cleanAssay` function can generate two functions to
remove NA's from rows and columns using the `by` argument based on which
it selects the appropriate 'mean' and subset functions. This maintains
the clarity of having the operation in the function name when used:
`cleanAssayRows` and `cleanAssayCols`.
