write_manifest <- function(manifest, run_dir) {
  temporary <- file.path(run_dir, "manifest.tmp")
  destination <- file.path(run_dir, "manifest.rds")
  saveRDS(manifest, temporary)
  if (!file.rename(temporary, destination)) {
    stop("Could not save the QC result.", call. = FALSE)
  }
}

run_worker <- function(matrix_path, run_dir) {
  tryCatch(
    {
      result <- OmicsQC::run_qc_pipeline(matrix_path = matrix_path)
      result$inputs$matrix <- "uploaded-assay.tsv"
      result$validate$results$source <- basename(result$validate$results$source)
      result$validate$config_errors$source <- basename(result$validate$config_errors$source)

      output_dir <- file.path(run_dir, "outputs")
      outputs <- OmicsQC::write_qc_outputs(result, output_dir)
      output_files <- unlist(outputs[names(outputs) != "output_dir"], use.names = FALSE)
      output_files <- basename(output_files[!is.na(output_files)])
      archive_path <- file.path(run_dir, "omicsqc_outputs.zip")
      zip::zipr(archive_path, output_files, root = output_dir)
      if (!file.exists(archive_path)) {
        stop("The output archive could not be created.", call. = FALSE)
      }

      write_manifest(
        list(
          state = "completed",
          status = result$status,
          dataset_qc = result$dataset_qc,
          built_in_checks = result$built_in_checks,
          validate_results = result$validate$results,
          validate_errors = result$validate$config_errors,
          removed_samples = result$removed$samples,
          removed_features = result$removed$features,
          files = output_files
        ),
        run_dir
      )
      0L
    },
    error = function(error) {
      message <- conditionMessage(error)
      for (path in c(matrix_path, run_dir)) {
        message <- gsub(path, "[temporary file]", message, fixed = TRUE)
      }
      write_manifest(list(state = "error", message = message), run_dir)
      1L
    }
  )
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 2L) {
    stop("Expected an uploaded assay path and run directory.", call. = FALSE)
  }
  quit(status = run_worker(args[[1]], args[[2]]))
}
