load_feature_links <- function(feature_linking_paths = NULL) {
  if (!length(feature_linking_paths)) {
    return(data.frame(
      source_feature_id = character(),
      target_feature_id = character(),
      stringsAsFactors = FALSE
    ))
  }

  link_frames <- lapply(feature_linking_paths, function(path) {
    link_frame <- read_delimited_file(path, sep = "\t")
    required_columns <- c("source_feature_id", "target_feature_id")
    if (!all(required_columns %in% names(link_frame))) {
      stop(
        "Feature linking files must contain source_feature_id and target_feature_id columns: ",
        path,
        call. = FALSE
      )
    }

    link_frame <- link_frame[required_columns]
    link_frame$source_feature_id <- normalize_ids(link_frame$source_feature_id)
    link_frame$target_feature_id <- normalize_ids(link_frame$target_feature_id)
    link_frame
  })

  combined <- do.call(rbind, link_frames)
  rownames(combined) <- NULL
  combined
}

harmonize_feature_ids <- function(filtered_assay, feature_links) {
  if (!nrow(filtered_assay)) {
    return(list(
      assay = filtered_assay,
      summary = data.frame(
        original_feature_id = character(),
        harmonized_feature_id = character(),
        mapping_status = character(),
        stringsAsFactors = FALSE
      )
    ))
  }

  if (!nrow(feature_links)) {
    return(list(
      assay = filtered_assay,
      summary = data.frame(
        original_feature_id = rownames(filtered_assay),
        harmonized_feature_id = rownames(filtered_assay),
        mapping_status = rep("unmapped", nrow(filtered_assay)),
        stringsAsFactors = FALSE
      )
    ))
  }

  split_links <- split(feature_links$target_feature_id, feature_links$source_feature_id)
  harmonized_ids <- vapply(rownames(filtered_assay), function(feature_id) {
    mapped <- unique(split_links[[feature_id]])
    if (!length(mapped) || length(mapped) > 1) {
      return(feature_id)
    }

    mapped[[1]]
  }, character(1))
  mapping_status <- vapply(rownames(filtered_assay), function(feature_id) {
    mapped <- unique(split_links[[feature_id]])
    if (!length(mapped)) {
      return("unmapped")
    }

    if (length(mapped) > 1) {
      return("ambiguous")
    }

    "mapped"
  }, character(1))

  rownames(filtered_assay) <- harmonized_ids

  list(
    assay = filtered_assay,
    summary = data.frame(
      original_feature_id = names(mapping_status),
      harmonized_feature_id = harmonized_ids,
      mapping_status = mapping_status,
      stringsAsFactors = FALSE
    )
  )
}

run_qc_pipeline <- function(
    matrix_path,
    metadata_path = NULL,
    feature_id_col = NULL,
    matrix_sep = "\t",
    rule_files = NULL,
    include_default_rules = TRUE,
    feature_linking_paths = NULL,
    sample_max_missing_fraction = 0.2,
    sample_max_zero_fraction = 0.8,
    sample_min_detected_features = 1,
    sample_min_total_signal = 0,
    feature_max_missing_fraction = 0.2,
    feature_max_zero_fraction = 0.8,
    feature_max_outlier_fraction = 0.2,
    enforce_normal_distribution = FALSE,
    normality_alpha = 0.05) {
  inputs <- load_omics_matrix(
    matrix_path = matrix_path,
    feature_id_col = feature_id_col,
    matrix_sep = matrix_sep
  )

  dataset_qc <- build_dataset_qc(inputs)
  sample_qc <- build_sample_qc_table(inputs)
  feature_qc <- build_feature_qc_table(inputs, normality_enforcement = enforce_normal_distribution)
  built_in_checks <- run_builtin_checks(inputs, sample_qc = sample_qc, feature_qc = feature_qc)

  validate_result <- empty_validate_result()
  filtered <- list(assay = NULL, sample_qc = NULL, feature_qc = NULL)
  removed <- list(
    samples = data.frame(sample_id = character(), reasons = character(), stringsAsFactors = FALSE),
    features = data.frame(feature_id = character(), reasons = character(), stringsAsFactors = FALSE)
  )
  harmonization <- list(
    summary = data.frame(
      original_feature_id = character(),
      harmonized_feature_id = character(),
      mapping_status = character(),
      stringsAsFactors = FALSE
    ),
    feature_links = data.frame(
      source_feature_id = character(),
      target_feature_id = character(),
      stringsAsFactors = FALSE
    )
  )

  if (!has_fatal_builtin_failures(built_in_checks)) {
    sample_filters <- evaluate_sample_filters(
      sample_qc = sample_qc,
      max_missing_fraction = sample_max_missing_fraction,
      max_zero_fraction = sample_max_zero_fraction,
      min_detected_features = sample_min_detected_features,
      min_total_signal = sample_min_total_signal
    )
    feature_filters <- evaluate_feature_filters(
      feature_qc = feature_qc,
      max_missing_fraction = feature_max_missing_fraction,
      max_zero_fraction = feature_max_zero_fraction,
      max_outlier_fraction = feature_max_outlier_fraction,
      enforce_normal_distribution = enforce_normal_distribution,
      normality_alpha = normality_alpha
    )

    removed <- list(samples = sample_filters$removed, features = feature_filters$removed)
    filtered <- filter_assay_matrix(inputs, sample_filters = sample_filters, feature_filters = feature_filters)
    filtered$sample_qc <- sample_qc[sample_filters$keep_mask, , drop = FALSE]
    filtered$feature_qc <- feature_qc[feature_filters$keep_mask, , drop = FALSE]

    feature_links <- load_feature_links(feature_linking_paths = feature_linking_paths)
    harmonized <- harmonize_feature_ids(filtered$assay, feature_links = feature_links)
    filtered$assay <- harmonized$assay
    harmonization <- list(summary = harmonized$summary, feature_links = feature_links)

    validate_result <- run_validate_checks(
      sample_qc = filtered$sample_qc,
      rule_files = rule_files,
      include_default_rules = include_default_rules
    )
  }

  overall_fail <- has_fatal_builtin_failures(built_in_checks) || has_validate_failures(validate_result)

  structure(
    list(
      status = if (overall_fail) "fail" else "pass",
      inputs = inputs$paths,
      settings = list(
        feature_id_col = inputs$feature_id_col,
        include_default_rules = include_default_rules,
        rule_files = rule_files,
        feature_linking_paths = feature_linking_paths,
        sample_max_missing_fraction = sample_max_missing_fraction,
        sample_max_zero_fraction = sample_max_zero_fraction,
        sample_min_detected_features = sample_min_detected_features,
        sample_min_total_signal = sample_min_total_signal,
        feature_max_missing_fraction = feature_max_missing_fraction,
        feature_max_zero_fraction = feature_max_zero_fraction,
        feature_max_outlier_fraction = feature_max_outlier_fraction,
        enforce_normal_distribution = enforce_normal_distribution,
        normality_alpha = normality_alpha,
        preprocessing_assumption = "Input matrix is assumed to be already normalized, transformed, and corrected."
      ),
      built_in_checks = built_in_checks,
      dataset_qc = dataset_qc,
      sample_qc = sample_qc,
      feature_qc = feature_qc,
      filtered = filtered,
      removed = removed,
      harmonization = harmonization,
      validate = validate_result
    ),
    class = "OmicsQCResult"
  )
}
