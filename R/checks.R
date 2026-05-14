empty_check_frame <- function() {
  data.frame(
    scope = character(),
    check = character(),
    severity = character(),
    status = character(),
    item_count = integer(),
    items = character(),
    detail = character(),
    stringsAsFactors = FALSE
  )
}

new_check_record <- function(scope, check, severity, status, items = character(), detail) {
  cleaned_items <- unique(items[!is.na(items) & nzchar(items)])

  data.frame(
    scope = scope,
    check = check,
    severity = severity,
    status = status,
    item_count = as.integer(length(cleaned_items)),
    items = paste(cleaned_items, collapse = ";"),
    detail = detail,
    stringsAsFactors = FALSE
  )
}

combine_check_frames <- function(frames) {
  frames <- Filter(Negate(is.null), frames)
  if (!length(frames)) {
    return(empty_check_frame())
  }

  combined <- do.call(rbind, frames)
  rownames(combined) <- NULL
  combined
}

compute_zero_fraction <- function(assay_matrix) {
  apply(assay_matrix, 2, function(values) {
    non_missing <- values[!is.na(values)]
    if (!length(non_missing)) {
      return(NA_real_)
    }

    sum(non_missing == 0) / length(non_missing)
  })
}

build_sample_qc_table <- function(inputs) {
  metadata_order <- match(inputs$sample_ids, inputs$metadata[[inputs$sample_id_col]])
  ordered_metadata <- inputs$metadata[metadata_order, , drop = FALSE]
  extra_metadata_cols <- setdiff(names(ordered_metadata), inputs$sample_id_col)

  sample_qc <- data.frame(
    sample_id = inputs$sample_ids,
    total_signal = colSums(inputs$assay, na.rm = TRUE),
    detected_features = colSums(!is.na(inputs$assay) & inputs$assay > 0),
    missing_fraction = colMeans(is.na(inputs$assay)),
    zero_fraction = compute_zero_fraction(inputs$assay),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (length(extra_metadata_cols)) {
    sample_qc[extra_metadata_cols] <- ordered_metadata[extra_metadata_cols]
  }

  sample_qc
}

run_builtin_checks <- function(inputs) {
  metadata_ids <- inputs$metadata[[inputs$sample_id_col]]
  matrix_ids <- inputs$sample_ids
  assay_matrix <- inputs$assay

  duplicate_metadata_ids <- unique(metadata_ids[duplicated(metadata_ids) & nzchar(metadata_ids)])
  duplicate_matrix_ids <- unique(matrix_ids[duplicated(matrix_ids) & nzchar(matrix_ids)])
  duplicate_feature_ids <- unique(inputs$feature_ids[duplicated(inputs$feature_ids) & nzchar(inputs$feature_ids)])
  missing_in_metadata <- setdiff(matrix_ids, metadata_ids)
  missing_in_matrix <- setdiff(metadata_ids, matrix_ids)

  non_numeric_items <- character()
  if (nrow(inputs$non_numeric_entries) > 0) {
    non_numeric_items <- paste0(
      colnames(assay_matrix)[inputs$non_numeric_entries[, "col"]],
      "/",
      rownames(assay_matrix)[inputs$non_numeric_entries[, "row"]]
    )
  }

  all_zero_samples <- colnames(assay_matrix)[apply(assay_matrix, 2, function(values) {
    non_missing <- values[!is.na(values)]
    length(non_missing) > 0 && all(non_missing == 0)
  })]

  all_missing_samples <- colnames(assay_matrix)[apply(assay_matrix, 2, function(values) {
    all(is.na(values))
  })]

  sample_mismatch_detail <- if (length(missing_in_metadata) || length(missing_in_matrix)) {
    paste(
      c(
        if (length(missing_in_metadata)) {
          paste("Missing in metadata:", paste(missing_in_metadata, collapse = ", "))
        },
        if (length(missing_in_matrix)) {
          paste("Missing in matrix:", paste(missing_in_matrix, collapse = ", "))
        }
      ),
      collapse = " | "
    )
  } else {
    "Matrix samples and metadata samples align."
  }

  combine_check_frames(list(
    new_check_record(
      scope = "dataset",
      check = "empty_assay_matrix",
      severity = "error",
      status = if (nrow(assay_matrix) == 0 || ncol(assay_matrix) == 0) "fail" else "pass",
      detail = if (nrow(assay_matrix) == 0 || ncol(assay_matrix) == 0) {
        "The assay matrix has no features or no samples."
      } else {
        sprintf("Loaded %d features across %d samples.", nrow(assay_matrix), ncol(assay_matrix))
      }
    ),
    new_check_record(
      scope = "dataset",
      check = "empty_metadata",
      severity = "error",
      status = if (nrow(inputs$metadata) == 0) "fail" else "pass",
      detail = if (nrow(inputs$metadata) == 0) {
        "Metadata contains no rows."
      } else {
        sprintf("Loaded %d metadata rows.", nrow(inputs$metadata))
      }
    ),
    new_check_record(
      scope = "matrix",
      check = "duplicate_matrix_sample_ids",
      severity = "error",
      status = if (length(duplicate_matrix_ids)) "fail" else "pass",
      items = duplicate_matrix_ids,
      detail = if (length(duplicate_matrix_ids)) {
        "Sample column names in the assay matrix must be unique."
      } else {
        "Assay matrix sample IDs are unique."
      }
    ),
    new_check_record(
      scope = "metadata",
      check = "duplicate_metadata_sample_ids",
      severity = "error",
      status = if (length(duplicate_metadata_ids)) "fail" else "pass",
      items = duplicate_metadata_ids,
      detail = if (length(duplicate_metadata_ids)) {
        "Metadata sample IDs must be unique."
      } else {
        "Metadata sample IDs are unique."
      }
    ),
    new_check_record(
      scope = "matrix",
      check = "duplicate_feature_ids",
      severity = "warning",
      status = if (length(duplicate_feature_ids)) "fail" else "pass",
      items = duplicate_feature_ids,
      detail = if (length(duplicate_feature_ids)) {
        "Feature identifiers are duplicated in the assay matrix."
      } else {
        "Feature identifiers are unique."
      }
    ),
    new_check_record(
      scope = "dataset",
      check = "sample_id_mismatch",
      severity = "error",
      status = if (length(missing_in_metadata) || length(missing_in_matrix)) "fail" else "pass",
      items = c(missing_in_metadata, missing_in_matrix),
      detail = sample_mismatch_detail
    ),
    new_check_record(
      scope = "matrix",
      check = "non_numeric_assay_values",
      severity = "warning",
      status = if (length(non_numeric_items)) "fail" else "pass",
      items = non_numeric_items,
      detail = if (length(non_numeric_items)) {
        "Non-numeric assay values were coerced to NA."
      } else {
        "All assay values are numeric."
      }
    ),
    new_check_record(
      scope = "sample",
      check = "all_zero_samples",
      severity = "warning",
      status = if (length(all_zero_samples)) "fail" else "pass",
      items = all_zero_samples,
      detail = if (length(all_zero_samples)) {
        "One or more samples contain only zeros across non-missing features."
      } else {
        "No all-zero samples were detected."
      }
    ),
    new_check_record(
      scope = "sample",
      check = "all_missing_samples",
      severity = "warning",
      status = if (length(all_missing_samples)) "fail" else "pass",
      items = all_missing_samples,
      detail = if (length(all_missing_samples)) {
        "One or more samples contain only missing values."
      } else {
        "No all-missing samples were detected."
      }
    )
  ))
}

has_fatal_builtin_failures <- function(checks) {
  any(checks$severity == "error" & checks$status == "fail")
}
