# CVDLINK Omics QC

OmicsQC is a generic R-based omics QC pipeline for different omic layers,
including transcriptomics, metabolomics, proteomics, and related data types.
It assumes that the input matrix is already optimally preprocessed,
normalized, transformed, and corrected according to best practices for the
specific omics layer. OmicsQC does not perform preprocessing itself.

The current implementation focuses on:

- validating a single tab-separated sample-by-feature matrix
- computing dataset-level, sample-level, and feature-level QC summaries
- applying configurable sample and feature removal filters
- optionally removing multivariate sample outliers with PCA
- integrating validate YAML rule files for custom checks
- harmonizing feature IDs with optional linking files
- writing a QCed harmonized dataset and a compact HTML QC report

## Input contract

The package expects one required tab-separated assay matrix where:

1. The first column contains sample identifiers.
2. Each remaining column is a feature identifier in the header row.
3. Every remaining cell contains a numeric measurement or a missing value.
4. A value of `0` is treated as below detection limit by default.

OmicsQC checks the file format automatically and records violations such as
missing identifiers, duplicate sample IDs, duplicate feature IDs, empty
matrices, and non-numeric values.

Optional feature linking files can be supplied to harmonize dataset-specific
feature IDs to a common identifier space. Each linking file must contain at
least `source_feature_id` and `target_feature_id` columns.

## Quick start

```r
library(OmicsQC)

result <- run_qc_pipeline(
  matrix_path = system.file("tmp", "extdata", "example_assay.tsv", package = "OmicsQC")
)

print(result)
write_qc_outputs(result, output_dir = "qc-output")
```

## Built-in diagnostics and filters

OmicsQC reports:

- dataset dimensions and global missing-value and zero-value rates
- sample-level missingness, below-detection-limit fraction, and signal summaries
- feature-level missingness, below-detection-limit fraction, distribution summaries, and outlier burden

It can automatically remove samples and features using adjustable thresholds,
including:

- samples with many missing values
- samples with many values below detection limit
- samples with extreme robust PCA score distances
- features with many missing values
- features with many values below detection limit
- features that fail optional normality enforcement

The exported QCed dataset contains only retained samples and retained features.
Audit tables record which samples and features were removed and why.

### PCA sample outliers

PCA-based sample removal is opt-in and runs after the scalar sample and feature
filters. PCA is fitted once using centered and scaled retained features. Only
candidate samples with complete measurements across those features participate;
incomplete samples are recorded as ineligible for PCA and are not removed by
this stage.

The pipeline retains the fewest principal components that explain at least 80%
of variance by default. Scores on retained components are standardized using
their median and median absolute deviation. A sample is removed when its
squared robust score distance exceeds the 99th percentile of a chi-square
distribution with degrees of freedom equal to the usable component count.

```r
result <- run_qc_pipeline(
  matrix_path = "assay.tsv",
  remove_pca_outliers = TRUE,
  pca_variance_target = 0.80,
  pca_cutoff_probability = 0.99
)
```

If fewer than three complete candidate samples, fewer than two variable
features, or no retained component with a usable robust scale remains, PCA is
skipped without failing the pipeline. The reason is available in
`result$pca$reason`. Per-sample eligibility, component scores, robust distances,
cutoffs, and flags are stored in `result$pca$scores`; explained variance is in
`result$pca$variance`.

## Custom checks with validate

OmicsQC computes metric tables that can be checked with additional validate
YAML rules. The default rule file targets sample-level metrics, and
user-supplied rules can be layered on top.

```r
result <- run_qc_pipeline(
  matrix_path = "assay.tsv",
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

## Harmonization

You can harmonize feature IDs across datasets with one or more linking files:

```r
result <- run_qc_pipeline(
  matrix_path = "assay.tsv",
  feature_linking_paths = "feature_links.tsv"
)
```

The resulting QCed assay export uses harmonized IDs where a unique mapping is
available and records mapped, unmapped, and ambiguous features in audit outputs.

## Reporting and outputs

`write_qc_outputs()` writes:

- the QCed harmonized assay matrix
- dataset, sample, and feature QC tables
- removed sample and feature audit tables
- PCA score and explained-variance diagnostics
- harmonization summary tables
- validate results
- a compact standardized HTML report suitable for cross-dataset comparison

## Runner script

You can also run the pipeline directly:

```sh
Rscript scripts/run_qc.R tmp/extdata/example_assay.tsv qc-output [rule.yaml ...]
```

## Contact

Urmo Võsa (urmo.vosa at gmail.com), Andres Veidenberg (andres.veidenberg at gmail.com)
