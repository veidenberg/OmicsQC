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
  if (!ncol(filtered_assay)) {
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
    original_feature_ids <- colnames(filtered_assay)
    return(list(
      assay = filtered_assay,
      summary = data.frame(
        original_feature_id = original_feature_ids,
        harmonized_feature_id = original_feature_ids,
        mapping_status = rep("unmapped", ncol(filtered_assay)),
        stringsAsFactors = FALSE
      )
    ))
  }

  split_links <- split(feature_links$target_feature_id, feature_links$source_feature_id)
  original_feature_ids <- colnames(filtered_assay)
  harmonized_ids <- vapply(original_feature_ids, function(feature_id) {
    mapped <- unique(split_links[[feature_id]])
    if (!length(mapped) || length(mapped) > 1) {
      return(feature_id)
    }

    mapped[[1]]
  }, character(1))
  mapping_status <- vapply(original_feature_ids, function(feature_id) {
    mapped <- unique(split_links[[feature_id]])
    if (!length(mapped)) {
      return("unmapped")
    }

    if (length(mapped) > 1) {
      return("ambiguous")
    }

    "mapped"
  }, character(1))

  colnames(filtered_assay) <- harmonized_ids

  list(
    assay = filtered_assay,
    summary = data.frame(
      original_feature_id = original_feature_ids,
      harmonized_feature_id = harmonized_ids,
      mapping_status = mapping_status,
      stringsAsFactors = FALSE
    )
  )
}

validate_pca_settings <- function(remove_pca_outliers, variance_target, cutoff_probability) {
  if (!is.logical(remove_pca_outliers) || length(remove_pca_outliers) != 1 || is.na(remove_pca_outliers)) {
    stop("remove_pca_outliers must be TRUE or FALSE.", call. = FALSE)
  }

  probability_settings <- list(
    pca_variance_target = variance_target,
    pca_cutoff_probability = cutoff_probability
  )
  for (setting_name in names(probability_settings)) {
    value <- probability_settings[[setting_name]]
    if (!is.numeric(value) || length(value) != 1 || !is.finite(value) || value <= 0 || value >= 1) {
      stop(setting_name, " must be a finite number strictly between 0 and 1.", call. = FALSE)
    }
  }
}

run_qc_pipeline <- function(
    matrix_path,
    metadata_path = NULL,
    sample_id_col = NULL,
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
    normality_alpha = 0.05,
    remove_pca_outliers = FALSE,
    pca_variance_target = 0.8,
    pca_cutoff_probability = 0.99) {
  validate_pca_settings(
    remove_pca_outliers = remove_pca_outliers,
    variance_target = pca_variance_target,
    cutoff_probability = pca_cutoff_probability
  )

  inputs <- load_omics_matrix(
    matrix_path = matrix_path,
    sample_id_col = sample_id_col,
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
  pca <- empty_pca_result(
    inputs$sample_ids,
    status = if (remove_pca_outliers) "skipped" else "disabled",
    reason = if (remove_pca_outliers) "PCA was not run because built-in input checks failed." else ""
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

    pca <- evaluate_pca_outliers(
      assay_matrix = inputs$assay[, feature_filters$keep_mask, drop = FALSE],
      enabled = remove_pca_outliers,
      variance_target = pca_variance_target,
      cutoff_probability = pca_cutoff_probability,
      candidate_mask = sample_filters$keep_mask
    )
    sample_filters$keep_mask <- sample_filters$keep_mask & pca$keep_mask
    removed_samples <- rbind(sample_filters$removed, pca$removed)
    rownames(removed_samples) <- NULL
    removed <- list(samples = removed_samples, features = feature_filters$removed)
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
        sample_id_col = inputs$sample_id_col,
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
        remove_pca_outliers = remove_pca_outliers,
        pca_variance_target = pca_variance_target,
        pca_cutoff_probability = pca_cutoff_probability,
        preprocessing_assumption = "Input matrix is assumed to be already normalized, transformed, and corrected."
      ),
      built_in_checks = built_in_checks,
      dataset_qc = dataset_qc,
      sample_qc = sample_qc,
      feature_qc = feature_qc,
      filtered = filtered,
      removed = removed,
      pca = pca,
      harmonization = harmonization,
      validate = validate_result
    ),
    class = "OmicsQCResult"
  )
}
