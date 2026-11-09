test_that("default rules and the synthetic example are installed", {
  rules <- system.file("validate", "base_rules.yaml", package = "OmicsQC")
  example <- system.file("extdata", "example_assay.tsv", package = "OmicsQC")

  expect_true(nzchar(rules))
  expect_true(file.exists(rules))
  expect_true(nzchar(example))
  expect_true(file.exists(example))

  result <- run_qc_pipeline(example)
  expect_identical(result$status, "pass")
  expect_gt(nrow(result$validate$results), 0)
  expect_equal(nrow(result$validate$config_errors), 0)
})
