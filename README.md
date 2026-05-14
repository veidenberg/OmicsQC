# CVDLINK Omics QC

OmicsQC is a small R-based omics QC pipeline for a generic feature-by-sample
matrix plus sample metadata table.

The first version focuses on:

- loading a delimited assay matrix and metadata table
- running built-in sanity checks such as duplicate IDs and sample mismatches
- computing simple sample-level QC summaries
- applying additional checks from validate YAML rule files

## Input contract

The package expects two files:

1. An assay matrix in TSV or CSV format where the first column contains feature
   identifiers and each remaining column is a sample.
2. A metadata table in TSV or CSV format with a sample identifier column named
   `sample_id` by default.

## Quick start

```r
library(OmicsQC)

result <- run_qc_pipeline(
  matrix_path = system.file("extdata", "example_assay.tsv", package = "OmicsQC"),
  metadata_path = system.file("extdata", "example_metadata.tsv", package = "OmicsQC")
)

print(result)
write_qc_outputs(result, output_dir = "qc-output")
```

## Custom checks with validate

OmicsQC calculates a canonical sample-level QC table with these columns:

- `sample_id`
- `total_signal`
- `detected_features`
- `missing_fraction`
- `zero_fraction`

You can add project-specific checks by supplying one or more validate YAML rule
files:

```r
result <- run_qc_pipeline(
  matrix_path = "assay.tsv",
  metadata_path = "metadata.tsv",
  rule_files = "custom_rules.yaml"
)
```

Example custom rule file:

```yaml
rules:
  - expr: missing_fraction <= 0.10
    name: missing_fraction_cap
    label: Missing fraction cap
    description: |
      Samples should have at most 10% missing values.
    meta:
      severity: warning
```

## Runner script

You can also run the pipeline directly:

```sh
Rscript scripts/run_qc.R inst/extdata/example_assay.tsv inst/extdata/example_metadata.tsv qc-output
```

## Contact

Urmo Võsa (urmo.vosa at gmail.com), Andres Veidenberg (andres.veidenberg at gmail.com)
