print.OmicsQCResult <- function(x, ...) {
  built_in_failures <- sum(x$built_in_checks$status == "fail")
  validate_failures <- if (nrow(x$validate$results)) {
    sum(x$validate$results$status == "fail")
  } else {
    0L
  }

  cat("OmicsQC result:", toupper(x$status), "\n")
  cat("Assay file:", x$inputs$matrix, "\n")
  cat("Metadata file:", x$inputs$metadata, "\n")
  cat("Built-in checks:", nrow(x$built_in_checks), "total,", built_in_failures, "failing\n")
  cat("Validate checks:", nrow(x$validate$results), "total,", validate_failures, "failing\n")
  cat("Validate config errors:", nrow(x$validate$config_errors), "\n")

  if (!is.null(x$sample_qc)) {
    cat("Sample summaries:", nrow(x$sample_qc), "samples\n")
  }

  if (built_in_failures) {
    cat("\nFailing built-in checks:\n")
    print(x$built_in_checks[x$built_in_checks$status == "fail", , drop = FALSE], row.names = FALSE)
  }

  if (validate_failures) {
    cat("\nFailing validate checks:\n")
    print(x$validate$results[x$validate$results$status == "fail", , drop = FALSE], row.names = FALSE)
  }

  if (nrow(x$validate$config_errors)) {
    cat("\nValidate configuration errors:\n")
    print(x$validate$config_errors, row.names = FALSE)
  }

  invisible(x)
}

write_tsv <- function(data, path) {
  utils::write.table(
    data,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
}

write_qc_outputs <- function(result, output_dir) {
  if (!inherits(result, "OmicsQCResult")) {
    stop("result must inherit from 'OmicsQCResult'.", call. = FALSE)
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  built_in_checks_path <- file.path(output_dir, "built_in_checks.tsv")
  validate_results_path <- file.path(output_dir, "validate_results.tsv")
  validate_errors_path <- file.path(output_dir, "validate_config_errors.tsv")
  summary_path <- file.path(output_dir, "summary.txt")
  sample_qc_path <- file.path(output_dir, "sample_qc.tsv")

  write_tsv(result$built_in_checks, built_in_checks_path)
  write_tsv(result$validate$results, validate_results_path)
  write_tsv(result$validate$config_errors, validate_errors_path)

  if (!is.null(result$sample_qc)) {
    write_tsv(result$sample_qc, sample_qc_path)
  }

  writeLines(capture.output(print(result)), con = summary_path)

  invisible(list(
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
    sample_qc = if (file.exists(sample_qc_path)) normalizePath(sample_qc_path, winslash = "/", mustWork = TRUE) else NA_character_,
    built_in_checks = normalizePath(built_in_checks_path, winslash = "/", mustWork = TRUE),
    validate_results = normalizePath(validate_results_path, winslash = "/", mustWork = TRUE),
    validate_config_errors = normalizePath(validate_errors_path, winslash = "/", mustWork = TRUE),
    summary = normalizePath(summary_path, winslash = "/", mustWork = TRUE)
  ))
}
