### verify errors {{{
err.verifySE <- "Input needs to be a SummarizedExperiment"
err.verifyCoords <- paste(
  "The SummarizedExperiment you have provided has no coordinates.\n",
  "Compartment extraction will fail.\n",
  "Please provide rowRanges with genomic coordinates for the object."
)
err.verifyAssayNames.beta <- paste("The 'assays' slot should contain", shQuote('Beta'), "for array data")
err.verifyAssayNames.atac_counts <- paste("The 'assays' slot should contain", shQuote('counts'), "for atac data")
err.verifyAssayNames.rna_counts <- paste("The 'assays' slot should contain", shQuote('counts'), "for rna data")
# }}}
