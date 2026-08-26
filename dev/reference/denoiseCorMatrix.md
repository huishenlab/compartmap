# Denoising of Covariance matrix using Random Matrix Theory

Denoising of Covariance matrix using Random Matrix Theory

## Usage

``` r
denoiseCorMatrix(cormat, rows, cols)
```

## Details

Given a random matrix \\X\_{m\times n}\\ with independently and
identically distributed values, the eigendecomposition of its covariance
matrix \\C\\ is

\$\$C = VDV^T\$\$

where \\D\\ is a diagonal matrix of its eigenvalues and \$V\$ contains
corresponding eigenvectors of \\C\\. The Marchenko-Pastur distribution
describes the distribution of random matrix eigenvalues, specifying that
they fall with the bounds

\$\$lambda\_{\pm }=\sigma ^{2}\left(1\pm
{\sqrt{\frac{m}{n}}}\right)^{2}\$\$

where \\\sigma^2\\ is the variance, equal to 1 in the case of
correlation matrices and \\m/n\\ the aspect ratio of the matrix. As
\\\lambda\_+\\ is the maximum expected eigenvalue in the
Marchenko-Pastur distribution and eigenvalues greater than
\\\lambda\_+\\ are unlikely to derive from noise, we interpret these as
signal. Given \\k\\ eigenvalues greater than \\\lambda\_+\\, we can
reconstruct a \`denoised' correlation matrix

\$\$C' = V_kD_kV_k^T\$\$

where \\V_k\\ and \\D_k\\ denote the eigenvectors and eigenvalues
respectively corresponding to the \\k\\ largest eigenvalues. In effect,
we start by assuming that the correlation matrix is random noise to find
the upper bound of its Marchenko-Pastur distribution. Eigenvalues above
this bound violate this assumption and allow us to reconstruct the
correlation matrix excluding the elements representing random noise.
