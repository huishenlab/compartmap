# Transform/normalize counts using TF-IDF

Transform/normalize counts using TF-IDF

## Usage

``` r
transformTFIDF(mat, scale.factor = 100000, binarize = FALSE)
```

## Arguments

- mat:

  n x p input matrix (n = samples/cells; p = rna counts)

- scale.factor:

  Scaling factor for the term-frequency (TF)

- binarize:

  Whether to binarize the input matrix: any value \> 0 is set to 1

## Value

A TF-IDF transformed matrix of the same dimensions as the input

## Details

This function and its helpers were modeled after or taken from:

- http://andrewjohnhill.com/images/posts/2019-5-6-dimensionality-reduction-for-scatac-data/analysis.html

- https://divingintogeneticsandgenomics.rbind.io/post/clustering-scatacseq-data-the-tf-idf-way/

## Examples

``` r
m <- 1000
n <- 100
mat <- round(matrix(runif(m * n), m, n))
# Input needs to be a tall matrix
tfidf <- transformTFIDF(mat)
```
