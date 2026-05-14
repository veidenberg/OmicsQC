test_that("sample filters remove samples with excessive missingness or zeros", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "clean_assay.tsv"),
    include_default_rules = FALSE,
    sample_max_missing_fraction = 0.15,
    sample_max_zero_fraction = 0.7
  )

  expect_true("sample_4" %in% result$removed$samples$sample_id)
  expect_true(!("sample_4" %in% rownames(result$filtered$assay)))
})

test_that("feature filters remove features with excessive zeros", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "clean_assay.tsv"),
    include_default_rules = FALSE,
    feature_max_zero_fraction = 0.8
  )

  expect_true("gene_e" %in% result$removed$features$feature_id)
  expect_true(!("gene_e" %in% colnames(result$filtered$assay)))
})
