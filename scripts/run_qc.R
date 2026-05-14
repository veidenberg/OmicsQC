#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat(
    "Usage: Rscript scripts/run_qc.R <assay.tsv|csv> <metadata.tsv|csv> <output_dir> [rule.yaml ...]\n",
    file = stderr()
  )
  quit(status = 1)
}

all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
script_path <- if (length(file_arg)) {
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = TRUE)
} else {
  normalizePath("scripts/run_qc.R", winslash = "/", mustWork = TRUE)
}

project_root <- dirname(dirname(script_path))
options(OmicsQC.package_root = project_root)

r_files <- sort(list.files(file.path(project_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE))
invisible(lapply(r_files, source))

result <- run_qc_pipeline(
  matrix_path = args[[1]],
  metadata_path = args[[2]],
  rule_files = if (length(args) > 3) args[4:length(args)] else NULL
)

print(result)
outputs <- write_qc_outputs(result, output_dir = args[[3]])
cat("Outputs written to:", outputs$output_dir, "\n")

quit(status = if (identical(result$status, "pass")) 0 else 1)
