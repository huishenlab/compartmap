# Package index

## Compartment mapping

- [`scCompartments()`](https://huishenlab.github.io/compartmap/reference/scCompartments.md)
  : Estimate A/B compartments from single-cell RNA or ATAC sequencing
  data
- [`arrayCompartments()`](https://huishenlab.github.io/compartmap/reference/arrayCompartments.md)
  : Estimate A/B compartments from methylation array data

## Compartment helper

- [`filterCompartments()`](https://huishenlab.github.io/compartmap/reference/filterCompartments.md)
  : Filter compartments using confidence estimates and eigenvalue
  thresholds
- [`fixCompartments()`](https://huishenlab.github.io/compartmap/reference/fixCompartments.md)
  : Invert, or "fix", compartments that have a minimum confidence score
  (1-min.conf)
- [`extractOpenClosed()`](https://huishenlab.github.io/compartmap/reference/extractOpenClosed.md)
  : Get the open and closed compartment calls based on sign of singular
  values
- [`getABSignal()`](https://huishenlab.github.io/compartmap/reference/getABSignal.md)
  : Calculate Pearson correlations of smoothed eigenvectors

## Correlation matrix

- [`estRMT()`](https://huishenlab.github.io/compartmap/reference/estRMT.md)
  : Denoising of Covariance matrix using Random Matrix Theory
- [`getCorMatrix()`](https://huishenlab.github.io/compartmap/reference/getCorMatrix.md)
  : Calculate Pearson correlations of a binned matrix
- [`getDenoisedCorMatrix()`](https://huishenlab.github.io/compartmap/reference/getDenoisedMatrix.md)
  : Wrapper to denoise a correlation matrix using a Random Matrix Theory
  approach

## Data transformation

- [`transformTFIDF()`](https://huishenlab.github.io/compartmap/reference/transformTFIDF.md)
  : Transform/normalize counts using TF-IDF
- [`hdf5TFIDF()`](https://huishenlab.github.io/compartmap/reference/hdf5TFIDF.md)
  : Transform/normalize compartment calls using TF-IDF on HDF5-backed
  objects

## Bootstrapping

- [`precomputeBootstrapMeans()`](https://huishenlab.github.io/compartmap/reference/precomputeBootstrapMeans.md)
  : Pre-compute the global means for bootstrapping compartments
- [`bootstrapCompartments()`](https://huishenlab.github.io/compartmap/reference/bootstrapCompartments.md)
  : Non-parametric bootstrapping of compartments and summarization of
  bootstraps/compute confidence intervals
- [`summarizeBootstraps()`](https://huishenlab.github.io/compartmap/reference/summarizeBootstraps.md)
  : Summarize the bootstrap compartment estimates and compute
  Agresti-Coull confidence intervals

## Data import

- [`importBigWig()`](https://huishenlab.github.io/compartmap/reference/importBigWig.md)
  : Import and optionally summarize a bigwig at a given resolution
- [`filterOpenSea()`](https://huishenlab.github.io/compartmap/reference/filterOpenSea.md)
  : Filter to open sea CpG loci
- [`getOpenSeas()`](https://huishenlab.github.io/compartmap/reference/getOpenSeas.md)
  : Gather open sea CpG from a GRanges of CpG islands
- [`preprocessArrays()`](https://huishenlab.github.io/compartmap/reference/preprocessArrays.md)
  : Preprocess arrays for compartment inference

## Bioconductor helpers

- [`condenseRE()`](https://huishenlab.github.io/compartmap/reference/condenseRE.md)
  : Condense a RaggedExperiment to a list of SummarizedExperiments
- [`condenseSE()`](https://huishenlab.github.io/compartmap/reference/condenseSE.md)
  : Condense the output of condenseRE to reconstruct per-sample GRanges
  objects to plot
- [`getAssayNames()`](https://huishenlab.github.io/compartmap/reference/getAssayNames.md)
  : Get the assay names from a SummarizedExperiment object
- [`getSeqLengths()`](https://huishenlab.github.io/compartmap/reference/getSeqLengths.md)
  : Get the seqlengths of a chromosome from a given genome's GRanges
- [`getChrs()`](https://huishenlab.github.io/compartmap/reference/getChrs.md)
  : Get the chromosomes from an object

## Stat functions

- [`agrestiCoullCI()`](https://huishenlab.github.io/compartmap/reference/agrestiCoullCI.md)
  : Agresti-Coull confidence interval for a binomial proportion
- [`fexpit()`](https://huishenlab.github.io/compartmap/reference/fexpit.md)
  : Helper function: expanded expit
- [`fisherZ()`](https://huishenlab.github.io/compartmap/reference/fisherZ.md)
  : Fisher's Z transformation
- [`flogit()`](https://huishenlab.github.io/compartmap/reference/flogit.md)
  : Helper function: squeezed logit
- [`getBinMatrix()`](https://huishenlab.github.io/compartmap/reference/getBinMatrix.md)
  : Generate bins for A/B compartment estimation
- [`getDomainInflections()`](https://huishenlab.github.io/compartmap/reference/getDomainInflections.md)
  : A wrapper function to generate a GRanges object of chromatin domain
  inflection points
- [`getGlobalMeans()`](https://huishenlab.github.io/compartmap/reference/getGlobalMeans.md)
  : Get the global means of a matrix
- [`getMatrixBlocks()`](https://huishenlab.github.io/compartmap/reference/getMatrixBlocks.md)
  : Get chunked sets of row-wise or column-wise indices of a matrix
- [`getSVD()`](https://huishenlab.github.io/compartmap/reference/getSVD.md)
  : Compute the SVD of a matrix using irlba
- [`getShrinkageTargets()`](https://huishenlab.github.io/compartmap/reference/getShrinkageTargets.md)
  : Get the specified samples to shrink towards instead of the global
  mean
- [`ifisherZ()`](https://huishenlab.github.io/compartmap/reference/ifisherZ.md)
  : Inverse Fisher's Z transformation
- [`imputeKNN()`](https://huishenlab.github.io/compartmap/reference/imputeKNN.md)
  : Impute missing values/NAs with KNN
- [`meanSmoother()`](https://huishenlab.github.io/compartmap/reference/meanSmoother.md)
  : Windowed mean smoother
- [`precomputeBootstrapMeans()`](https://huishenlab.github.io/compartmap/reference/precomputeBootstrapMeans.md)
  : Pre-compute the global means for bootstrapping compartments
- [`preprocessArrays()`](https://huishenlab.github.io/compartmap/reference/preprocessArrays.md)
  : Preprocess arrays for compartment inference
- [`removeEmptyBoots()`](https://huishenlab.github.io/compartmap/reference/removeEmptyBoots.md)
  : Remove bootstrap estimates that failed
- [`scCompartments()`](https://huishenlab.github.io/compartmap/reference/scCompartments.md)
  : Estimate A/B compartments from single-cell RNA or ATAC sequencing
  data
- [`shrinkBins()`](https://huishenlab.github.io/compartmap/reference/shrinkBins.md)
  : Employ an eBayes shrinkage approach for bin-level estimates for A/B
  inference
- [`sparseToDenseMatrix()`](https://huishenlab.github.io/compartmap/reference/sparseToDenseMatrix.md)
  : Convert a sparse matrix to a dense matrix in a block-wise fashion

## Plotting

- [`plotAB()`](https://huishenlab.github.io/compartmap/reference/plotAB.md)
  : Plots A/B compartment estimates on a per chromosome basis
- [`plotCorMatrix()`](https://huishenlab.github.io/compartmap/reference/plotCorMatrix.md)
  : Plot a denoised correlation matrix

## Helpers

- [`cleanAssayCols()`](https://huishenlab.github.io/compartmap/reference/cleanAssayCols.md)
  :

  Remove columns/cells/samples with NAs exceeding a threshold. See
  [`cleanAssay()`](https://huishenlab.github.io/compartmap/reference/cleanAssay.md)

- [`cleanAssayRows()`](https://huishenlab.github.io/compartmap/reference/cleanAssayRows.md)
  :

  Remove rows with NAs exceeding a threshold. See
  [`cleanAssay()`](https://huishenlab.github.io/compartmap/reference/cleanAssay.md)

- [`.n_approx()`](https://huishenlab.github.io/compartmap/reference/dot-n_approx.md)
  : n_tilde in AC

- [`.p_approx()`](https://huishenlab.github.io/compartmap/reference/dot-p_approx.md)
  : p_tilde in AC

- [`.z()`](https://huishenlab.github.io/compartmap/reference/dot-z.md) :
  Normal alpha/2 quantile

## Datasets

- [`hg19.gr`](https://huishenlab.github.io/compartmap/reference/hg19.gr.md)
  : hg19 seqlengths as a GRanges object
- [`hg38.gr`](https://huishenlab.github.io/compartmap/reference/hg38.gr.md)
  : hg38 seqlengths as a GRanges object
- [`mm10.gr`](https://huishenlab.github.io/compartmap/reference/mm10.gr.md)
  : mm10 seqlengths as a GRanges object
- [`mm9.gr`](https://huishenlab.github.io/compartmap/reference/mm9.gr.md)
  : mm9 seqlengths as a GRanges object
- [`openSeas.hg19`](https://huishenlab.github.io/compartmap/reference/openSeas.hg19.md)
  : hg19 open sea CpG as a GRanges object
- [`openSeas.hg38`](https://huishenlab.github.io/compartmap/reference/openSeas.hg38.md)
  : hg38 open sea CpG as a GRanges object
- [`openSeas.mm10`](https://huishenlab.github.io/compartmap/reference/openSeas.mm10.md)
  : mm10 open sea CpG as a GRanges object
- [`openSeas.mm9`](https://huishenlab.github.io/compartmap/reference/openSeas.mm9.md)
  : mm9 open sea CpG as a GRanges object
- [`array_data_chr14`](https://huishenlab.github.io/compartmap/reference/array.data.chr14.md)
  : Example Illumina 450k methylation array data for compartmap
- [`k562_scatac_chr14`](https://huishenlab.github.io/compartmap/reference/k562_scatac_chr14.md)
  : Example scATAC-seq data for compartmap
- [`k562_scrna_chr14`](https://huishenlab.github.io/compartmap/reference/k562_scrna_chr14.md)
  : Example scRNA-seq data for compartmap
- [`k562_scrna_raw`](https://huishenlab.github.io/compartmap/reference/k562_scrna_se_chr14.md)
  : Example scRNA-seq data for compartmap
- [`ss3_umi_sce`](https://huishenlab.github.io/compartmap/reference/ss3_umi_sce.md)
  : Example SMART-seq3 scRNA-seq data for compartmap
