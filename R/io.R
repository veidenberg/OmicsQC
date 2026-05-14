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

load_omics_inputs <- function(
    matrix_path,
    metadata_path,
    sample_id_col = "sample_id",
    feature_id_col = NULL,
    matrix_sep = NULL,
    metadata_sep = NULL) {
  assay_raw <- read_delimited_file(matrix_path, sep = matrix_sep)
  metadata <- read_delimited_file(metadata_path, sep = metadata_sep)

  if (ncol(assay_raw) < 2) {
    stop(
      "The assay matrix must contain a feature column and at least one sample column.",
      call. = FALSE
    )
  }

  if (!(sample_id_col %in% names(metadata))) {
    stop("Metadata is missing the sample identifier column: ", sample_id_col, call. = FALSE)
  }

  feature_column <- if (is.null(feature_id_col)) names(assay_raw)[1] else feature_id_col
  if (!(feature_column %in% names(assay_raw))) {
    stop("Assay matrix is missing the feature identifier column: ", feature_column, call. = FALSE)
  }

  sample_columns <- setdiff(names(assay_raw), feature_column)
  if (!length(sample_columns)) {
    stop("No sample columns were found in the assay matrix.", call. = FALSE)
  }

  feature_ids <- normalize_ids(assay_raw[[feature_column]])
  sample_ids <- normalize_ids(sample_columns)

  assay_values <- assay_raw[sample_columns]
  colnames(assay_values) <- sample_ids

  coerced <- coerce_assay_to_numeric(assay_values)
  assay_matrix <- coerced$matrix
  rownames(assay_matrix) <- feature_ids
  colnames(assay_matrix) <- sample_ids

  metadata[[sample_id_col]] <- normalize_ids(metadata[[sample_id_col]])

  structure(
    list(
      assay = assay_matrix,
      metadata = metadata,
      feature_ids = feature_ids,
      sample_ids = sample_ids,
      feature_id_col = feature_column,
      sample_id_col = sample_id_col,
      non_numeric_entries = coerced$non_numeric,
      paths = list(
        matrix = normalizePath(matrix_path, winslash = "/", mustWork = TRUE),
        metadata = normalizePath(metadata_path, winslash = "/", mustWork = TRUE)
      )
    ),
    class = "OmicsQCInputs"
  )
}
