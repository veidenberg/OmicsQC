run_qc_pipeline <- function(
    matrix_path,
    metadata_path,
    sample_id_col = "sample_id",
    feature_id_col = NULL,
    matrix_sep = NULL,
    metadata_sep = NULL,
    rule_files = NULL,
    include_default_rules = TRUE) {
  inputs <- load_omics_inputs(
    matrix_path = matrix_path,
    metadata_path = metadata_path,
    sample_id_col = sample_id_col,
    feature_id_col = feature_id_col,
    matrix_sep = matrix_sep,
    metadata_sep = metadata_sep
  )

  built_in_checks <- run_builtin_checks(inputs)
  sample_qc <- NULL
  validate_result <- empty_validate_result()

  if (!has_fatal_builtin_failures(built_in_checks)) {
    sample_qc <- build_sample_qc_table(inputs)
    validate_result <- run_validate_checks(
      sample_qc = sample_qc,
      rule_files = rule_files,
      include_default_rules = include_default_rules
    )
  }

  overall_fail <- any(built_in_checks$status == "fail") || has_validate_failures(validate_result)

  structure(
    list(
      status = if (overall_fail) "fail" else "pass",
      inputs = inputs$paths,
      settings = list(
        sample_id_col = sample_id_col,
        feature_id_col = inputs$feature_id_col,
        include_default_rules = include_default_rules,
        rule_files = rule_files
      ),
      built_in_checks = built_in_checks,
      sample_qc = sample_qc,
      validate = validate_result
    ),
    class = "OmicsQCResult"
  )
}
