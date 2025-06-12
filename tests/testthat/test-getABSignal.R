# .isCompartmentOpen
test_that(".isCompartmentOpen", {
  .isCompartmentOpen <- compartmap:::.isCompartmentOpen
  expect_true(
    .isCompartmentOpen(is.atac_or_rna = TRUE, eigen = 1, cutoff = 0),
    .isCompartmentOpen(is.atac_or_rna = FALSE, eigen = 0, cutoff = 0),
    .isCompartmentOpen(is.atac_or_rna = FALSE, eigen = -1, cutoff = 0)
  )
  expect_false(
    .isCompartmentOpen(is.atac_or_rna = TRUE, eigen = 0, cutoff = 0),
    .isCompartmentOpen(is.atac_or_rna = TRUE, eigen = -1, cutoff = 0),
    .isCompartmentOpen(is.atac_or_rna = FALSE, eigen = 1, cutoff = 0)
  )
})
