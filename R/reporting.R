html_escape <- function(text) {
  escaped <- gsub("&", "&amp;", as.character(text), fixed = TRUE)
  escaped <- gsub("<", "&lt;", escaped, fixed = TRUE)
  escaped <- gsub(">", "&gt;", escaped, fixed = TRUE)
  escaped
}

render_html_table <- function(data) {
  if (is.null(data) || !nrow(data)) {
    return("<p>None</p>")
  }

  header <- paste0(
    "<tr>",
    paste(sprintf("<th>%s</th>", html_escape(names(data))), collapse = ""),
    "</tr>"
  )
  rows <- vapply(seq_len(nrow(data)), function(index) {
    paste0(
      "<tr>",
      paste(sprintf("<td>%s</td>", html_escape(data[index, , drop = TRUE])), collapse = ""),
      "</tr>"
    )
  }, character(1))

  paste0("<table>", header, paste(rows, collapse = ""), "</table>")
}

print.OmicsQCResult <- function(x, ...) {
  built_in_failures <- sum(x$built_in_checks$status == "fail")
  validate_failures <- if (nrow(x$validate$results)) {
    sum(x$validate$results$status == "fail")
  } else {
    0L
  }

  cat("OmicsQC result:", toupper(x$status), "\n")
  cat("Assay file:", x$inputs$matrix, "\n")
  cat("Built-in checks:", nrow(x$built_in_checks), "total,", built_in_failures, "failing\n")
  cat("Validate checks:", nrow(x$validate$results), "total,", validate_failures, "failing\n")
  cat("Validate config errors:", nrow(x$validate$config_errors), "\n")
  cat("Dataset metrics:", nrow(x$dataset_qc), "\n")

  if (!is.null(x$sample_qc)) {
    cat("Sample summaries:", nrow(x$sample_qc), "samples\n")
  }

  if (!is.null(x$feature_qc)) {
    cat("Feature summaries:", nrow(x$feature_qc), "features\n")
  }

  if (!is.null(x$filtered$assay)) {
    cat("Filtered assay dimensions:", nrow(x$filtered$assay), "features x", ncol(x$filtered$assay), "samples\n")
  }

  if (nrow(x$removed$samples)) {
    cat("Removed samples:", nrow(x$removed$samples), "\n")
  }

  if (nrow(x$removed$features)) {
    cat("Removed features:", nrow(x$removed$features), "\n")
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
  dataset_qc_path <- file.path(output_dir, "dataset_qc.tsv")
  validate_results_path <- file.path(output_dir, "validate_results.tsv")
  validate_errors_path <- file.path(output_dir, "validate_config_errors.tsv")
  summary_path <- file.path(output_dir, "summary.txt")
  sample_qc_path <- file.path(output_dir, "sample_qc.tsv")
  feature_qc_path <- file.path(output_dir, "feature_qc.tsv")
  removed_samples_path <- file.path(output_dir, "removed_samples.tsv")
  removed_features_path <- file.path(output_dir, "removed_features.tsv")
  harmonization_path <- file.path(output_dir, "harmonization.tsv")
  filtered_assay_path <- file.path(output_dir, "assay_qc_harmonized.tsv")
  html_report_path <- file.path(output_dir, "qc_report.html")

  write_tsv(result$built_in_checks, built_in_checks_path)
  write_tsv(result$dataset_qc, dataset_qc_path)
  write_tsv(result$validate$results, validate_results_path)
  write_tsv(result$validate$config_errors, validate_errors_path)

  if (!is.null(result$sample_qc)) {
    write_tsv(result$sample_qc, sample_qc_path)
  }

  if (!is.null(result$feature_qc)) {
    write_tsv(result$feature_qc, feature_qc_path)
  }

  write_tsv(result$removed$samples, removed_samples_path)
  write_tsv(result$removed$features, removed_features_path)

  if (!is.null(result$harmonization$summary)) {
    write_tsv(result$harmonization$summary, harmonization_path)
  }

  if (!is.null(result$filtered$assay)) {
    filtered_export <- data.frame(
      feature_id = rownames(result$filtered$assay),
      result$filtered$assay,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    write_tsv(filtered_export, filtered_assay_path)
  }

  writeLines(capture.output(print(result)), con = summary_path)

  html_lines <- c(
    "<html><head><meta charset='utf-8'><title>OmicsQC Report</title>",
    "<style>body{font-family:Helvetica,Arial,sans-serif;margin:24px;}table{border-collapse:collapse;width:100%;margin:12px 0;}th,td{border:1px solid #d7d7d7;padding:6px 8px;text-align:left;}h1,h2{margin-bottom:8px;} .metric-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:12px;} .metric-card{border:1px solid #d7d7d7;padding:12px;border-radius:8px;background:#fafafa;}</style>",
    "</head><body>",
    "<h1>OmicsQC Report</h1>",
    sprintf("<p><strong>Status:</strong> %s</p>", html_escape(toupper(result$status))),
    "<h2>Dataset Summary</h2>",
    "<div class='metric-grid'>",
    vapply(seq_len(nrow(result$dataset_qc)), function(index) {
      sprintf(
        "<div class='metric-card'><strong>%s</strong><br/>%s</div>",
        html_escape(result$dataset_qc$metric[[index]]),
        html_escape(format(result$dataset_qc$value[[index]], digits = 4))
      )
    }, character(1)),
    "</div>",
    "<h2>Built-in Checks</h2>",
    render_html_table(result$built_in_checks),
    "<h2>Removed Samples</h2>",
    render_html_table(result$removed$samples),
    "<h2>Removed Features</h2>",
    render_html_table(result$removed$features),
    "<h2>Harmonization</h2>",
    render_html_table(result$harmonization$summary),
    "<h2>Validate Results</h2>",
    render_html_table(result$validate$results),
    "</body></html>"
  )
  writeLines(html_lines, con = html_report_path)

  invisible(list(
    output_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
    dataset_qc = normalizePath(dataset_qc_path, winslash = "/", mustWork = TRUE),
    sample_qc = if (file.exists(sample_qc_path)) normalizePath(sample_qc_path, winslash = "/", mustWork = TRUE) else NA_character_,
    feature_qc = if (file.exists(feature_qc_path)) normalizePath(feature_qc_path, winslash = "/", mustWork = TRUE) else NA_character_,
    built_in_checks = normalizePath(built_in_checks_path, winslash = "/", mustWork = TRUE),
    removed_samples = normalizePath(removed_samples_path, winslash = "/", mustWork = TRUE),
    removed_features = normalizePath(removed_features_path, winslash = "/", mustWork = TRUE),
    harmonization = if (file.exists(harmonization_path)) normalizePath(harmonization_path, winslash = "/", mustWork = TRUE) else NA_character_,
    filtered_assay = if (file.exists(filtered_assay_path)) normalizePath(filtered_assay_path, winslash = "/", mustWork = TRUE) else NA_character_,
    validate_results = normalizePath(validate_results_path, winslash = "/", mustWork = TRUE),
    validate_config_errors = normalizePath(validate_errors_path, winslash = "/", mustWork = TRUE),
    summary = normalizePath(summary_path, winslash = "/", mustWork = TRUE),
    html_report = normalizePath(html_report_path, winslash = "/", mustWork = TRUE)
  ))
}
