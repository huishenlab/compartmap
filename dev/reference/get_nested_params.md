# Set outer and inner params for nester parallelization The outer param is across the input samples/columns and the second is for bootstrapping. If `boot.parallel` is FALSE, the inner param is set to `SerialParam`.

Set outer and inner params for nester parallelization The outer param is
across the input samples/columns and the second is for bootstrapping. If
`boot.parallel` is FALSE, the inner param is set to `SerialParam`.

## Usage

``` r
get_nested_params(BPPARAM, boot.parallel)
```
