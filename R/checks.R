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
  apply(assay_matrix, 1, function(values) {
    non_missing <- values[!is.na(values)]
    if (!length(non_missing)) {
      return(NA_real_)
    }

    sum(non_missing == 0) / length(non_missing)
  })
}

compute_feature_zero_fraction <- function(assay_matrix) {
  apply(assay_matrix, 2, function(values) {
    non_missing <- values[!is.na(values)]
    if (!length(non_missing)) {
      return(NA_real_)
    }

    sum(non_missing == 0) / length(non_missing)
  })
}

safe_stat <- function(values, fun) {
  observed <- values[!is.na(values)]
  if (!length(observed)) {
    return(NA_real_)
  }

  fun(observed)
}

compute_feature_outlier_fraction <- function(assay_matrix, multiplier = 1.5) {
  apply(assay_matrix, 2, function(values) {
    observed <- values[!is.na(values)]
    if (length(observed) < 4) {
      return(NA_real_)
    }

    quartiles <- stats::quantile(observed, probs = c(0.25, 0.75), names = FALSE)
    iqr_value <- quartiles[[2]] - quartiles[[1]]
    if (is.na(iqr_value) || iqr_value == 0) {
      return(0)
    }

    lower_bound <- quartiles[[1]] - multiplier * iqr_value
    upper_bound <- quartiles[[2]] + multiplier * iqr_value
    sum(observed < lower_bound | observed > upper_bound) / length(observed)
  })
}

compute_normality_pvalue <- function(assay_matrix) {
  apply(assay_matrix, 2, function(values) {
    observed <- values[!is.na(values)]
    if (length(observed) < 3 || length(observed) > 5000 || length(unique(observed)) < 3) {
      return(NA_real_)
    }

    stats::shapiro.test(observed)$p.value
  })
}

build_dataset_qc <- function(inputs) {
  assay_matrix <- inputs$assay
  observed <- assay_matrix[!is.na(assay_matrix)]

  data.frame(
    metric = c(
      "n_features",
      "n_samples",
      "missing_fraction_overall",
      "zero_fraction_overall",
      "sample_missing_fraction_median",
      "feature_missing_fraction_median",
      "sample_zero_fraction_median",
      "feature_zero_fraction_median"
    ),
    value = c(
      ncol(assay_matrix),
      nrow(assay_matrix),
      mean(is.na(assay_matrix)),
      if (!length(observed)) NA_real_ else sum(observed == 0) / length(observed),
      stats::median(rowMeans(is.na(assay_matrix))),
      stats::median(colMeans(is.na(assay_matrix))),
      stats::median(compute_zero_fraction(assay_matrix), na.rm = TRUE),
      stats::median(compute_feature_zero_fraction(assay_matrix), na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
}

build_sample_qc_table <- function(inputs) {
  data.frame(
    sample_id = inputs$sample_ids,
    total_signal = rowSums(inputs$assay, na.rm = TRUE),
    detected_features = rowSums(!is.na(inputs$assay) & inputs$assay > 0),
    missing_fraction = rowMeans(is.na(inputs$assay)),
    zero_fraction = compute_zero_fraction(inputs$assay),
    mean_signal = apply(inputs$assay, 1, safe_stat, fun = mean),
    median_signal = apply(inputs$assay, 1, safe_stat, fun = stats::median),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

build_feature_qc_table <- function(inputs, normality_enforcement = FALSE) {
  assay_matrix <- inputs$assay
  normality_pvalue <- if (isTRUE(normality_enforcement)) {
    compute_normality_pvalue(assay_matrix)
  } else {
    rep(NA_real_, ncol(assay_matrix))
  }

  data.frame(
    feature_id = inputs$feature_ids,
    missing_fraction = colMeans(is.na(assay_matrix)),
    zero_fraction = compute_feature_zero_fraction(assay_matrix),
    mean_signal = apply(assay_matrix, 2, safe_stat, fun = mean),
    median_signal = apply(assay_matrix, 2, safe_stat, fun = stats::median),
    sd_signal = apply(assay_matrix, 2, safe_stat, fun = stats::sd),
    outlier_fraction = compute_feature_outlier_fraction(assay_matrix),
    normality_pvalue = normality_pvalue,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

run_builtin_checks <- function(inputs, sample_qc, feature_qc) {
  all_zero_samples <- sample_qc$sample_id[!is.na(sample_qc$zero_fraction) & sample_qc$zero_fraction == 1]
  all_missing_samples <- sample_qc$sample_id[!is.na(sample_qc$missing_fraction) & sample_qc$missing_fraction == 1]
  all_zero_features <- feature_qc$feature_id[!is.na(feature_qc$zero_fraction) & feature_qc$zero_fraction == 1]
  all_missing_features <- feature_qc$feature_id[!is.na(feature_qc$missing_fraction) & feature_qc$missing_fraction == 1]

  combine_check_frames(list(
    inputs$input_issues,
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
    ),
    new_check_record(
      scope = "feature",
      check = "all_zero_features",
      severity = "warning",
      status = if (length(all_zero_features)) "fail" else "pass",
      items = all_zero_features,
      detail = if (length(all_zero_features)) {
        "One or more features are zero across all observed samples."
      } else {
        "No all-zero features were detected."
      }
    ),
    new_check_record(
      scope = "feature",
      check = "all_missing_features",
      severity = "warning",
      status = if (length(all_missing_features)) "fail" else "pass",
      items = all_missing_features,
      detail = if (length(all_missing_features)) {
        "One or more features are missing across all samples."
      } else {
        "No all-missing features were detected."
      }
    )
  ))
}

has_fatal_builtin_failures <- function(checks) {
  any(checks$severity == "error" & checks$status == "fail")
}

empty_pca_result <- function(sample_ids, status, reason = "") {
  list(
    status = status,
    reason = reason,
    n_complete_samples = 0L,
    n_variable_features = 0L,
    n_retained_pcs = 0L,
    cutoff = NA_real_,
    variance = data.frame(
      component = character(),
      explained_variance = numeric(),
      cumulative_variance = numeric(),
      retained = logical(),
      used_for_distance = logical(),
      stringsAsFactors = FALSE
    ),
    scores = data.frame(
      sample_id = sample_ids,
      analysis_status = rep(if (identical(status, "disabled")) "disabled" else "not_analyzed", length(sample_ids)),
      robust_distance_squared = rep(NA_real_, length(sample_ids)),
      cutoff = rep(NA_real_, length(sample_ids)),
      outlier = rep(FALSE, length(sample_ids)),
      stringsAsFactors = FALSE
    ),
    keep_mask = rep(TRUE, length(sample_ids)),
    removed = data.frame(sample_id = character(), reasons = character(), stringsAsFactors = FALSE)
  )
}

evaluate_pca_outliers <- function(
    assay_matrix,
    enabled = FALSE,
    variance_target = 0.8,
    cutoff_probability = 0.99,
    candidate_mask = rep(TRUE, nrow(assay_matrix))) {
  sample_ids <- rownames(assay_matrix)
  if (!isTRUE(enabled)) {
    result <- empty_pca_result(sample_ids, status = "disabled")
    result$scores$analysis_status[!candidate_mask] <- "scalar_filtered"
    return(result)
  }

  result <- empty_pca_result(sample_ids, status = "skipped")
  complete_case_mask <- if (ncol(assay_matrix)) {
    stats::complete.cases(assay_matrix)
  } else {
    rep(TRUE, nrow(assay_matrix))
  }
  complete_mask <- candidate_mask & complete_case_mask
  result$n_complete_samples <- sum(complete_mask)
  result$scores$analysis_status[!candidate_mask] <- "scalar_filtered"
  result$scores$analysis_status[candidate_mask & !complete_mask] <- "incomplete"
  result$scores$analysis_status[complete_mask] <- "eligible"

  if (result$n_complete_samples < 3) {
    result$reason <- "PCA requires at least three complete samples."
    return(result)
  }

  complete_assay <- assay_matrix[complete_mask, , drop = FALSE]
  feature_sd <- apply(complete_assay, 2, stats::sd)
  variable_feature_mask <- is.finite(feature_sd) & feature_sd > 0
  result$n_variable_features <- sum(variable_feature_mask)
  if (result$n_variable_features < 2) {
    result$reason <- "PCA requires at least two variable features."
    return(result)
  }

  pca <- stats::prcomp(
    complete_assay[, variable_feature_mask, drop = FALSE],
    center = TRUE,
    scale. = TRUE
  )
  component_variance <- pca$sdev^2
  positive_component_mask <- is.finite(component_variance) & component_variance > 0
  if (!any(positive_component_mask)) {
    result$reason <- "PCA produced no components with positive variance."
    return(result)
  }

  component_variance <- component_variance[positive_component_mask]
  component_scores <- pca$x[, positive_component_mask, drop = FALSE]
  explained_variance <- component_variance / sum(component_variance)
  cumulative_variance <- cumsum(explained_variance)
  retained_count <- which(cumulative_variance >= variance_target)[1]
  retained_mask <- seq_along(component_variance) <= retained_count
  retained_scores <- component_scores[, retained_mask, drop = FALSE]

  score_centers <- apply(retained_scores, 2, stats::median)
  score_scales <- apply(retained_scores, 2, stats::mad)
  usable_component_mask <- is.finite(score_scales) & score_scales > 0

  result$variance <- data.frame(
    component = colnames(component_scores),
    explained_variance = explained_variance,
    cumulative_variance = cumulative_variance,
    retained = retained_mask,
    used_for_distance = retained_mask & c(usable_component_mask, rep(FALSE, sum(!retained_mask))),
    stringsAsFactors = FALSE
  )
  result$n_retained_pcs <- retained_count

  score_frame <- as.data.frame(retained_scores, stringsAsFactors = FALSE)
  score_rows <- match(sample_ids[complete_mask], result$scores$sample_id)
  for (component in names(score_frame)) {
    result$scores[[component]] <- NA_real_
    result$scores[[component]][score_rows] <- score_frame[[component]]
  }

  if (!any(usable_component_mask)) {
    result$reason <- "Retained PCA components have no usable robust scale."
    return(result)
  }

  standardized_scores <- sweep(
    retained_scores[, usable_component_mask, drop = FALSE],
    2,
    score_centers[usable_component_mask],
    "-"
  )
  standardized_scores <- sweep(
    standardized_scores,
    2,
    score_scales[usable_component_mask],
    "/"
  )
  robust_distance_squared <- rowSums(standardized_scores^2)
  cutoff <- stats::qchisq(cutoff_probability, df = sum(usable_component_mask))
  outlier <- robust_distance_squared > cutoff

  result$status <- "completed"
  result$reason <- ""
  result$cutoff <- cutoff
  result$scores$analysis_status[score_rows] <- "analyzed"
  result$scores$robust_distance_squared[score_rows] <- robust_distance_squared
  result$scores$cutoff[score_rows] <- cutoff
  result$scores$outlier[score_rows] <- outlier
  result$keep_mask <- !result$scores$outlier
  result$removed <- data.frame(
    sample_id = result$scores$sample_id[result$scores$outlier],
    reasons = rep("pca_robust_distance>cutoff", sum(result$scores$outlier)),
    stringsAsFactors = FALSE
  )
  result
}

evaluate_sample_filters <- function(
    sample_qc,
    max_missing_fraction = 0.2,
    max_zero_fraction = 0.8,
    min_detected_features = 1,
    min_total_signal = 0) {
  reason_map <- vector("list", nrow(sample_qc))
  names(reason_map) <- sample_qc$sample_id

  register_reason <- function(ids, reason) {
    for (sample_id in ids) {
      reason_map[[sample_id]] <<- unique(c(reason_map[[sample_id]], reason))
    }
  }

  register_reason(
    sample_qc$sample_id[!is.na(sample_qc$missing_fraction) & sample_qc$missing_fraction > max_missing_fraction],
    sprintf("missing_fraction>%.3f", max_missing_fraction)
  )
  register_reason(
    sample_qc$sample_id[!is.na(sample_qc$zero_fraction) & sample_qc$zero_fraction > max_zero_fraction],
    sprintf("zero_fraction>%.3f", max_zero_fraction)
  )
  register_reason(
    sample_qc$sample_id[!is.na(sample_qc$detected_features) & sample_qc$detected_features < min_detected_features],
    sprintf("detected_features<%d", as.integer(min_detected_features))
  )
  register_reason(
    sample_qc$sample_id[!is.na(sample_qc$total_signal) & sample_qc$total_signal <= min_total_signal],
    sprintf("total_signal<=%.3f", min_total_signal)
  )

  removed_ids <- names(reason_map)[vapply(reason_map, length, integer(1)) > 0]
  audit <- data.frame(
    sample_id = sample_qc$sample_id,
    keep = !(sample_qc$sample_id %in% removed_ids),
    reasons = vapply(reason_map[sample_qc$sample_id], paste, collapse = ";", character(1)),
    stringsAsFactors = FALSE
  )

  list(
    keep_mask = audit$keep,
    removed = audit[!audit$keep, c("sample_id", "reasons"), drop = FALSE],
    audit = audit
  )
}

evaluate_feature_filters <- function(
    feature_qc,
    max_missing_fraction = 0.2,
    max_zero_fraction = 0.8,
    max_outlier_fraction = 0.2,
    enforce_normal_distribution = FALSE,
    normality_alpha = 0.05) {
  reason_map <- vector("list", nrow(feature_qc))
  names(reason_map) <- feature_qc$feature_id

  register_reason <- function(ids, reason) {
    for (feature_id in ids) {
      reason_map[[feature_id]] <<- unique(c(reason_map[[feature_id]], reason))
    }
  }

  register_reason(
    feature_qc$feature_id[!is.na(feature_qc$missing_fraction) & feature_qc$missing_fraction > max_missing_fraction],
    sprintf("missing_fraction>%.3f", max_missing_fraction)
  )
  register_reason(
    feature_qc$feature_id[!is.na(feature_qc$zero_fraction) & feature_qc$zero_fraction > max_zero_fraction],
    sprintf("zero_fraction>%.3f", max_zero_fraction)
  )
  register_reason(
    feature_qc$feature_id[!is.na(feature_qc$outlier_fraction) & feature_qc$outlier_fraction > max_outlier_fraction],
    sprintf("outlier_fraction>%.3f", max_outlier_fraction)
  )

  if (isTRUE(enforce_normal_distribution)) {
    register_reason(
      feature_qc$feature_id[!is.na(feature_qc$normality_pvalue) & feature_qc$normality_pvalue < normality_alpha],
      sprintf("normality_pvalue<%.3f", normality_alpha)
    )
  }

  removed_ids <- names(reason_map)[vapply(reason_map, length, integer(1)) > 0]
  audit <- data.frame(
    feature_id = feature_qc$feature_id,
    keep = !(feature_qc$feature_id %in% removed_ids),
    reasons = vapply(reason_map[feature_qc$feature_id], paste, collapse = ";", character(1)),
    stringsAsFactors = FALSE
  )

  list(
    keep_mask = audit$keep,
    removed = audit[!audit$keep, c("feature_id", "reasons"), drop = FALSE],
    audit = audit
  )
}

filter_assay_matrix <- function(inputs, sample_filters, feature_filters) {
  assay_matrix <- inputs$assay[sample_filters$keep_mask, feature_filters$keep_mask, drop = FALSE]

  list(
    assay = assay_matrix,
    sample_ids = rownames(assay_matrix),
    feature_ids = colnames(assay_matrix)
  )
}
