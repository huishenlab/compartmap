  
![](reference/figures/compartmap_logo.png)

*Chromatin compartment inference from single-cell RNA- and ATAC-seq*

Compartmap extends methods proposed by [Fortin and Hansen 2015, Genome
Biology](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0741-y)
to perform direct inference of higher-order chromatin in *single cells*
from scRNA-seq and scATAC-seq. Originally, Fortin and Hansen
demonstrated that chromatin conformation could be inferred from
(sc)ATAC-seq, bisulfite sequencing, DNase-seq and methylation arrays,
similar to the results provided by HiC at the group level. Thus, in
addition to the base information provided by the aforementioned assays,
chromatin state could also be inferred.

Here, we propose a method to infer both group and single-cell level
higher-order chromatin states from scRNA-seq and scATAC-seq. To
accomplish this, we employ a James-Stein estimator (JSE) towards a
global or targeted mean, using either chromsome or genome-wide
information from scRNA-seq and scATAC-seq. Additionally, due to the
sparsity of single-cell data, we employ a bootstrap procedure to
quantify the uncertainty associated with the state and boundaries of
inferred compartments. The output from compartmap can then be visualized
directly, compared with orthogonal assay types, and/or embedded with
something like UMAP or t-SNE. Further, to explore the higher-order
interacting domains inferred from compartmap, we use a Random Matrix
Theory (RMT) approach to resolve the “plaid-like” patterning, similar to
what is observed in Hi-C and scHi-C.

## Installation

### Bioconductor

    # requires BiocManager
    install.packages("BiocManager")

    # Release
    BiocManager::install("compartmap")

    # Development
    install.packages("BiocManager")
    BiocManager::install("huishenlab/compartmap")

### GitHub

You can install the development version from Github by cloning the repo
and running

``` bash
git clone https://github.com/huishenlab/compartmap
R CMD INSTALL compartmap
```

You can also use
[`BiocManager`](https://bioconductor.github.io/BiocManager/),
[`devtools`](https://devtools.r-lib.org/) or
[`pak`](https://pak.r-lib.org/):

``` r
BiocManager::install("huishenlab/compartmap")

devtools::install_github("huishenlab/compartmap")

pak::pkg_install("huishenlab/compartmap")
```

### Usage

A user guide is available on the [package
website](https://huishenlab.github.io/compartmap). Bug reports may be
submitted through GitHub issues.
