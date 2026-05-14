test_that("load_omics_inputs preserves sample and feature identifiers", {
  inputs <- load_omics_inputs(
    matrix_path = test_path("fixtures", "clean_assay.tsv"),
    metadata_path = test_path("fixtures", "clean_metadata.tsv")
  )

  expect_equal(inputs$sample_ids, c("sample_1", "sample_2", "sample_3"))
  expect_equal(inputs$feature_ids, c("gene_a", "gene_b", "gene_c", "gene_d"))
  expect_equal(dim(inputs$assay), c(4, 3))
})
