# Get the start and end of sequences within a vector of integers

It returns a data.table where rows may be duplicated to allow for
additional row-specific data if necessary

## Usage

``` r
get_sequential_idx(v)
```

## Arguments

- v:

  Indices of significant Mahalanobis distances

## Examples

``` r
v <- c(1:4, 5, 7:9)
get_sequential_idx(v)
#>    start_idx end_idx dc_id
#>        <num>   <num> <int>
#> 1:         1       5     1
#> 2:         1       5     1
#> 3:         1       5     1
#> 4:         1       5     1
#> 5:         1       5     1
#> 6:         7       9     2
#> 7:         7       9     2
#> 8:         7       9     2
```
