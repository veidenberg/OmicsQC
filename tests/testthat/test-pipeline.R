test_that("clean fixture passes end-to-end with default rules", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "clean_assay.tsv")
  )

  expect_identical(result$status, "pass")
  expect_equal(nrow(result$sample_qc), 4)
  expect_equal(nrow(result$feature_qc), 6)
  expect_true(nrow(result$dataset_qc) > 0)
  expect_true(all(result$validate$results$status == "pass"))
})

test_that("feature linking harmonizes the filtered assay output", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "clean_assay.tsv"),
    include_default_rules = FALSE,
    feature_linking_paths = test_path("fixtures", "feature_links.tsv")
  )

  expect_true(any(result$harmonization$summary$mapping_status == "mapped"))
  expect_true("GENE:B" %in% rownames(result$filtered$assay))
})

test_that("write_qc_outputs writes summary artifacts", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "clean_assay.tsv"),
    include_default_rules = FALSE
  )

  output_dir <- tempfile("omicsqc-output-")
  outputs <- write_qc_outputs(result, output_dir = output_dir)

  expect_true(file.exists(outputs$built_in_checks))
  expect_true(file.exists(outputs$sample_qc))
  expect_true(file.exists(outputs$feature_qc))
  expect_true(file.exists(outputs$filtered_assay))
  expect_true(file.exists(outputs$html_report))
  expect_true(file.exists(outputs$summary))
})
