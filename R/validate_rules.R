empty_validate_result <- function() {
  list(
    rule_files = character(),
    results = data.frame(
      source = character(),
      rule = character(),
      expression = character(),
      passes = integer(),
      fails = integer(),
      n_na = integer(),
      error = logical(),
      warning = logical(),
      status = character(),
      stringsAsFactors = FALSE
    ),
    config_errors = data.frame(
      source = character(),
      stage = character(),
      message = character(),
      stringsAsFactors = FALSE
    )
  )
}

find_package_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current, "DESCRIPTION"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      return("")
    }

    current <- parent
  }
}

default_rule_file_path <- function() {
  installed_path <- system.file("validate", "base_rules.yaml", package = "OmicsQC")
  if (nzchar(installed_path)) {
    return(installed_path)
  }

  candidate_paths <- character()

  package_root_option <- getOption("OmicsQC.package_root", "")
  if (nzchar(package_root_option)) {
    candidate_paths <- c(
      candidate_paths,
      file.path(package_root_option, "inst", "validate", "base_rules.yaml")
    )
  }

  detected_root <- find_package_root(getwd())
  if (nzchar(detected_root)) {
    candidate_paths <- c(
      candidate_paths,
      file.path(detected_root, "inst", "validate", "base_rules.yaml")
    )
  }

  existing_paths <- unique(candidate_paths[file.exists(candidate_paths)])
  if (!length(existing_paths)) {
    stop("Could not locate the default validate rule file.", call. = FALSE)
  }

  normalizePath(existing_paths[[1]], winslash = "/", mustWork = TRUE)
}

resolve_rule_files <- function(rule_files = NULL, include_default_rules = TRUE) {
  resolved <- character()

  if (isTRUE(include_default_rules)) {
    resolved <- c(resolved, default_rule_file_path())
  }

  if (length(rule_files)) {
    user_paths <- vapply(rule_files, function(path) {
      if (!file.exists(path)) {
        stop("Rule file does not exist: ", path, call. = FALSE)
      }

      normalizePath(path, winslash = "/", mustWork = TRUE)
    }, character(1))

    resolved <- c(resolved, user_paths)
  }

  unique(resolved)
}

summarize_validate_confrontation <- function(summary_frame, source) {
  rule_column <- if ("rule" %in% names(summary_frame)) "rule" else "name"
  expression_column <- if ("expression" %in% names(summary_frame)) "expression" else "expr"

  summarized <- data.frame(
    source = source,
    rule = as.character(summary_frame[[rule_column]]),
    expression = as.character(summary_frame[[expression_column]]),
    passes = as.integer(summary_frame$passes),
    fails = as.integer(summary_frame$fails),
    n_na = as.integer(summary_frame$nNA),
    error = as.logical(summary_frame$error),
    warning = as.logical(summary_frame$warning),
    status = ifelse(summary_frame$error | summary_frame$fails > 0, "fail", "pass"),
    stringsAsFactors = FALSE
  )

  rownames(summarized) <- NULL
  summarized
}

confrontation_summary_frame <- function(confrontation) {
  summary_method <- methods::selectMethod("summary", signature = "validation")
  summary_method(confrontation)
}

run_validate_checks <- function(sample_qc, rule_files = NULL, include_default_rules = TRUE) {
  validate_result <- empty_validate_result()

  resolved_files <- tryCatch(
    resolve_rule_files(rule_files = rule_files, include_default_rules = include_default_rules),
    error = function(error) error
  )

  if (inherits(resolved_files, "error")) {
    validate_result$config_errors <- data.frame(
      source = NA_character_,
      stage = "resolve",
      message = conditionMessage(resolved_files),
      stringsAsFactors = FALSE
    )
    return(validate_result)
  }

  validate_result$rule_files <- resolved_files
  if (!length(resolved_files)) {
    return(validate_result)
  }

  result_frames <- list()
  config_errors <- list()

  for (rule_file in resolved_files) {
    validator_object <- tryCatch(
      validate::validator(.file = rule_file),
      error = function(error) error
    )

    if (inherits(validator_object, "error")) {
      config_errors[[length(config_errors) + 1]] <- data.frame(
        source = rule_file,
        stage = "load",
        message = conditionMessage(validator_object),
        stringsAsFactors = FALSE
      )
      next
    }

    confrontation <- tryCatch(
      validate::confront(sample_qc, validator_object),
      error = function(error) error
    )

    if (inherits(confrontation, "error")) {
      config_errors[[length(config_errors) + 1]] <- data.frame(
        source = rule_file,
        stage = "confront",
        message = conditionMessage(confrontation),
        stringsAsFactors = FALSE
      )
      next
    }

    result_frames[[length(result_frames) + 1]] <- summarize_validate_confrontation(
      confrontation_summary_frame(confrontation),
      source = rule_file
    )
  }

  if (length(result_frames)) {
    validate_result$results <- do.call(rbind, result_frames)
    rownames(validate_result$results) <- NULL
  }

  if (length(config_errors)) {
    validate_result$config_errors <- do.call(rbind, config_errors)
    rownames(validate_result$config_errors) <- NULL
  }

  validate_result
}

has_validate_failures <- function(validate_result) {
  nrow(validate_result$config_errors) > 0 || any(validate_result$results$status == "fail")
}
