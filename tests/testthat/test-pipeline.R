test_that("clean fixture passes end-to-end with default rules", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "clean_assay.tsv"),
    metadata_path = test_path("fixtures", "clean_metadata.tsv")
  )

  expect_identical(result$status, "pass")
  expect_equal(nrow(result$sample_qc), 3)
  expect_true(all(result$built_in_checks$status == "pass"))
  expect_true(all(result$validate$results$status == "pass"))
})

test_that("write_qc_outputs writes summary artifacts", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "clean_assay.tsv"),
    metadata_path = test_path("fixtures", "clean_metadata.tsv"),
    include_default_rules = FALSE
  )

  output_dir <- tempfile("omicsqc-output-")
  outputs <- write_qc_outputs(result, output_dir = output_dir)

  expect_true(file.exists(outputs$built_in_checks))
  expect_true(file.exists(outputs$sample_qc))
  expect_true(file.exists(outputs$summary))
})
