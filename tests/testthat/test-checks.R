test_that("duplicate metadata IDs fail built-in checks", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "clean_assay.tsv"),
    metadata_path = test_path("fixtures", "duplicate_metadata.tsv"),
    include_default_rules = FALSE
  )

  duplicate_row <- result$built_in_checks[
    result$built_in_checks$check == "duplicate_metadata_sample_ids",
    ,
    drop = FALSE
  ]

  expect_identical(result$status, "fail")
  expect_equal(duplicate_row$status, "fail")
  expect_match(duplicate_row$items, "sample_2")
  expect_null(result$sample_qc)
})

test_that("sample mismatches are reported separately from duplicates", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "clean_assay.tsv"),
    metadata_path = test_path("fixtures", "mismatch_metadata.tsv"),
    include_default_rules = FALSE
  )

  mismatch_row <- result$built_in_checks[
    result$built_in_checks$check == "sample_id_mismatch",
    ,
    drop = FALSE
  ]

  expect_identical(result$status, "fail")
  expect_equal(mismatch_row$status, "fail")
  expect_match(mismatch_row$detail, "Missing in metadata")
  expect_match(mismatch_row$detail, "sample_3")
  expect_match(mismatch_row$detail, "sample_4")
})
