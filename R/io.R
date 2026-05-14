detect_delimiter <- function(path, sep = NULL) {
  if (!is.null(sep)) {
    return(sep)
  }

  extension <- tolower(tools::file_ext(path))
  if (identical(extension, "csv")) "," else "\t"
}

read_delimited_file <- function(path, sep = NULL) {
  if (!file.exists(path)) {
    stop("File does not exist: ", path, call. = FALSE)
  }

  utils::read.delim(
    file = path,
    sep = detect_delimiter(path, sep),
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

normalize_ids <- function(values) {
  trimws(as.character(values))
}

coerce_assay_to_numeric <- function(raw_assay) {
  value_matrix <- as.matrix(raw_assay)
  numeric_matrix <- suppressWarnings(
    matrix(
      as.numeric(value_matrix),
      nrow = nrow(value_matrix),
      ncol = ncol(value_matrix),
      dimnames = dimnames(value_matrix)
    )
  )

  non_numeric <- which(
    is.na(numeric_matrix) & !is.na(value_matrix) & trimws(value_matrix) != "",
    arr.ind = TRUE
  )

  list(matrix = numeric_matrix, non_numeric = non_numeric)
}

build_input_issues <- function(assay_matrix, feature_ids, sample_ids, non_numeric_entries) {
  duplicate_sample_ids <- unique(sample_ids[duplicated(sample_ids) & nzchar(sample_ids)])
  duplicate_feature_ids <- unique(feature_ids[duplicated(feature_ids) & nzchar(feature_ids)])
  missing_sample_ids <- sample_ids[!nzchar(sample_ids)]
  missing_feature_ids <- feature_ids[!nzchar(feature_ids)]

  non_numeric_items <- character()
  if (nrow(non_numeric_entries) > 0) {
    non_numeric_items <- paste0(
      sample_ids[non_numeric_entries[, "col"]],
      "/",
      feature_ids[non_numeric_entries[, "row"]]
    )
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
      scope = "matrix",
      check = "missing_sample_ids",
      severity = "error",
      status = if (length(missing_sample_ids)) "fail" else "pass",
      items = missing_sample_ids,
      detail = if (length(missing_sample_ids)) {
        "Every sample column must have a non-empty identifier."
      } else {
        "All sample identifiers are present."
      }
    ),
    new_check_record(
      scope = "matrix",
      check = "missing_feature_ids",
      severity = "error",
      status = if (length(missing_feature_ids)) "fail" else "pass",
      items = missing_feature_ids,
      detail = if (length(missing_feature_ids)) {
        "Every feature row must have a non-empty identifier."
      } else {
        "All feature identifiers are present."
      }
    ),
    new_check_record(
      scope = "matrix",
      check = "duplicate_sample_ids",
      severity = "error",
      status = if (length(duplicate_sample_ids)) "fail" else "pass",
      items = duplicate_sample_ids,
      detail = if (length(duplicate_sample_ids)) {
        "Sample identifiers in the matrix header must be unique."
      } else {
        "Sample identifiers are unique."
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
    )
  ))
}

load_omics_matrix <- function(
    matrix_path,
    feature_id_col = NULL,
    matrix_sep = "\t") {
  assay_raw <- read_delimited_file(matrix_path, sep = matrix_sep)

  if (ncol(assay_raw) < 2) {
    stop(
      "The assay matrix must contain a feature column and at least one sample column.",
      call. = FALSE
    )
  }

  feature_column <- if (is.null(feature_id_col)) names(assay_raw)[1] else feature_id_col
  if (!(feature_column %in% names(assay_raw))) {
    stop("Assay matrix is missing the feature identifier column: ", feature_column, call. = FALSE)
  }

  feature_column_index <- match(feature_column, names(assay_raw))
  sample_column_indices <- setdiff(seq_along(assay_raw), feature_column_index)
  if (!length(sample_column_indices)) {
    stop("No sample columns were found in the assay matrix.", call. = FALSE)
  }

  sample_columns <- names(assay_raw)[sample_column_indices]

  feature_ids <- normalize_ids(assay_raw[[feature_column]])
  sample_ids <- normalize_ids(sample_columns)

  assay_values <- assay_raw[sample_column_indices]
  colnames(assay_values) <- sample_ids

  coerced <- coerce_assay_to_numeric(assay_values)
  assay_matrix <- coerced$matrix
  rownames(assay_matrix) <- feature_ids
  colnames(assay_matrix) <- sample_ids

  structure(
    list(
      assay = assay_matrix,
      feature_ids = feature_ids,
      sample_ids = sample_ids,
      feature_id_col = feature_column,
      non_numeric_entries = coerced$non_numeric,
      paths = list(
        matrix = normalizePath(matrix_path, winslash = "/", mustWork = TRUE)
      ),
      input_issues = build_input_issues(
        assay_matrix = assay_matrix,
        feature_ids = feature_ids,
        sample_ids = sample_ids,
        non_numeric_entries = coerced$non_numeric
      )
    ),
    class = "OmicsQCInputs"
  )
}

load_omics_inputs <- function(
    matrix_path,
    metadata_path = NULL,
    sample_id_col = "sample_id",
    feature_id_col = NULL,
    matrix_sep = "\t",
    metadata_sep = NULL) {
  if (!is.null(metadata_path)) {
    warning("metadata_path is ignored; OmicsQC now expects a single assay matrix input.", call. = FALSE)
  }

  load_omics_matrix(
    matrix_path = matrix_path,
    feature_id_col = feature_id_col,
    matrix_sep = matrix_sep
  )
}
