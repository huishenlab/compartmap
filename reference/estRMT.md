# Denoising of Covariance matrix using Random Matrix Theory

Denoising of Covariance matrix using Random Matrix Theory

## Usage

``` r
estRMT(
  R,
  Q = NA,
  cutoff = c("max", "each"),
  eigenTreat = c("average", "delete"),
  numEig = 1
)
```

## Arguments

- R:

  input matrix

- Q:

  ratio of rows/size. Can be supplied externally or fit using data

- cutoff:

  takes two values max/each. If cutoff is max, Q is fitted and cutoff
  for eigenvalues is calculated. If cutoff is each, Q is set to
  row/size. Individual cutoff for each eigenvalue is calculated and used
  for filteration.

- eigenTreat:

  takes 2 values, average/delete. If average then the noisy eigenvalues
  are averged and each value is replaced by average. If delete then
  noisy eigenvalues are ignored and the diagonal entries of the
  correlation matrix are replaced with 1 to make the matrix psd.

- numEig:

  number of eigenvalues that are known for variance calculation. Default
  is set to 1. If numEig = 0 then variance is assumed to be 1.

## Value

A denoised RMT object

## Details

This method takes in data as a matrix object. It then fits a marchenko
pastur density to eigenvalues of the correlation matrix. All eigenvalues
above the cutoff are retained and ones below the cutoff are replaced
such that the trace of the correlation matrix is 1 or non-significant
eigenvalues are deleted and diagonal of correlation matrix is changed
to 1. Finally, correlation matrix is converted to covariance matrix.
This function was taken and modified from the covmat package
(https://github.com/cran/covmat) which has since been deprecated on
CRAN.

## Author

Rohit Arora

## Examples

``` r
rand_cor_mat <- cor(matrix(rnorm(100), nrow = 10))
denoised_rand_cor_mat <- estRMT(rand_cor_mat)$cov
       
```
