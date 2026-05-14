test_that("custom validate rules can fail sample summaries", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "clean_assay.tsv"),
    metadata_path = test_path("fixtures", "clean_metadata.tsv"),
    include_default_rules = FALSE,
    rule_files = test_path("fixtures", "strict_rules.yaml")
  )

  expect_identical(result$status, "fail")
  expect_equal(nrow(result$validate$config_errors), 0)
  expect_true(any(result$validate$results$status == "fail"))
  expect_true(any(result$validate$results$rule == "detected_features_floor"))
})
