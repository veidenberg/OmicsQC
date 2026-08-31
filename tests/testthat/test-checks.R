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

test_that("PCA filtering removes a robust score outlier and retains incomplete samples", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "pca_outlier_assay.tsv"),
    include_default_rules = FALSE,
    remove_pca_outliers = TRUE
  )

  expect_identical(result$pca$status, "completed")
  expect_identical(result$pca$n_retained_pcs, 1L)
  expect_equal(result$pca$cutoff, stats::qchisq(0.99, df = 1))
  expect_true(result$pca$variance$retained[[1]])
  expect_true(result$pca$variance$cumulative_variance[[1]] >= 0.8)
  expect_true("outlier_sample" %in% result$removed$samples$sample_id)
  expect_false("outlier_sample" %in% rownames(result$filtered$assay))
  expect_true("incomplete_sample" %in% rownames(result$filtered$assay))

  incomplete_diagnostic <- result$pca$scores[result$pca$scores$sample_id == "incomplete_sample", ]
  expect_identical(incomplete_diagnostic$analysis_status, "incomplete")
  expect_false(incomplete_diagnostic$outlier)
})

test_that("PCA removals compose with scalar sample filters", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "pca_outlier_assay.tsv"),
    include_default_rules = FALSE,
    sample_max_missing_fraction = 0.1,
    remove_pca_outliers = TRUE
  )

  expect_true(all(c("incomplete_sample", "outlier_sample") %in% result$removed$samples$sample_id))
  incomplete_diagnostic <- result$pca$scores[result$pca$scores$sample_id == "incomplete_sample", ]
  expect_identical(incomplete_diagnostic$analysis_status, "scalar_filtered")
})

test_that("PCA filtering is disabled by default", {
  result <- run_qc_pipeline(
    matrix_path = test_path("fixtures", "pca_outlier_assay.tsv"),
    include_default_rules = FALSE
  )

  expect_identical(result$pca$status, "disabled")
  expect_true("outlier_sample" %in% rownames(result$filtered$assay))
  expect_false(any(result$pca$scores$outlier))
})

test_that("PCA filtering skips data with no usable robust score scale", {
  assay <- rbind(
    sample_1 = c(0, 0),
    sample_2 = c(0, 0),
    sample_3 = c(1, 1)
  )

  result <- OmicsQC:::evaluate_pca_outliers(assay, enabled = TRUE)

  expect_identical(result$status, "skipped")
  expect_match(result$reason, "robust scale")
  expect_true(all(result$keep_mask))
})

test_that("PCA settings reject invalid values", {
  expect_error(
    run_qc_pipeline(
      test_path("fixtures", "clean_assay.tsv"),
      remove_pca_outliers = NA
    ),
    "TRUE or FALSE"
  )
  expect_error(
    run_qc_pipeline(
      test_path("fixtures", "clean_assay.tsv"),
      pca_variance_target = 1
    ),
    "strictly between 0 and 1"
  )
})
