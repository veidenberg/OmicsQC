test_that("load_omics_matrix preserves sample and feature identifiers", {
  inputs <- load_omics_matrix(
    matrix_path = test_path("fixtures", "clean_assay.tsv")
  )

  expect_equal(inputs$sample_ids, c("sample_1", "sample_2", "sample_3", "sample_4"))
  expect_equal(inputs$feature_ids, c("gene_a", "gene_b", "gene_c", "gene_d", "gene_e", "gene_f"))
  expect_equal(dim(inputs$assay), c(6, 4))
})

test_that("duplicate sample identifiers are surfaced as input issues", {
  inputs <- load_omics_matrix(
    matrix_path = test_path("fixtures", "duplicate_sample_assay.tsv")
  )

  duplicate_row <- inputs$input_issues[
    inputs$input_issues$check == "duplicate_sample_ids",
    ,
    drop = FALSE
  ]

  expect_equal(duplicate_row$status, "fail")
  expect_match(duplicate_row$detail, "must be unique")
})
