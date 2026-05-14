#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat(
    "Usage: Rscript scripts/run_qc.R <assay.tsv> <output_dir> [rule.yaml ...]\n",
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
  rule_files = if (length(args) > 2) args[3:length(args)] else NULL
)

print(result)
outputs <- write_qc_outputs(result, output_dir = args[[2]])
cat("Outputs written to:", outputs$output_dir, "\n")

quit(status = if (identical(result$status, "pass")) 0 else 1)
