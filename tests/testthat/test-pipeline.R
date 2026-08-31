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
  expect_true("GENE:B" %in% colnames(result$filtered$assay))
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
  expect_true(file.exists(outputs$pca_scores))
  expect_true(file.exists(outputs$pca_variance))
  expect_true(file.exists(outputs$html_report))
  expect_true(file.exists(outputs$summary))

  report_html <- paste(readLines(outputs$html_report, warn = FALSE), collapse = "\n")
  expect_match(report_html, "PCA Outlier Analysis", fixed = TRUE)
  expect_match(report_html, "PCA Sample Diagnostics", fixed = TRUE)
})

test_that("PCA diagnostics and outputs reflect filtered assay membership", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "pca_outlier_assay.tsv"),
    include_default_rules = FALSE,
    remove_pca_outliers = TRUE
  )
  output_dir <- tempfile("omicsqc-pca-output-")
  outputs <- write_qc_outputs(result, output_dir = output_dir)
  exported_scores <- utils::read.delim(outputs$pca_scores, check.names = FALSE)

  expect_true(exported_scores$outlier[exported_scores$sample_id == "outlier_sample"])
  expect_false("outlier_sample" %in% rownames(result$filtered$assay))
  expect_identical(result$settings$pca_variance_target, 0.8)
  expect_identical(result$settings$pca_cutoff_probability, 0.99)
})
