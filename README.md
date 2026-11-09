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
  matrix_path = system.file("extdata", "example_assay.tsv", package = "OmicsQC")
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
their median and median absolute deviation. By default, a sample is removed
when its squared robust score distance exceeds the mean squared robust score
distance plus two standard deviations. Set `pca_cutoff_method = "percentile"`
to use a chi-square percentile instead.

For a one-PC analysis, a two-sided $k$ SD score cutoff corresponds to the
chi-square percentile `pchisq(k^2, df = 1)`: 2 SD is about 95.45%, 2.58 SD is
about 99%, and 3 SD is about 99.73%. With multiple PCs, the percentile cutoff
uses the usable component count as the chi-square degrees of freedom, so there
is no single SD equivalent.

```r
result <- run_qc_pipeline(
  matrix_path = "assay.tsv",
  remove_pca_outliers = TRUE,
  pca_variance_target = 0.80,
  pca_cutoff_sd = 2
)

percentile_result <- run_qc_pipeline(
  matrix_path = "assay.tsv",
  remove_pca_outliers = TRUE,
  pca_cutoff_method = "percentile",
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
Rscript scripts/run_qc.R inst/extdata/example_assay.tsv qc-output [rule.yaml ...]
```

## CVDLink integration test

The `shiny/` application is a local integration-test interface for the
existing default QC pipeline. It is not a clinical decision system. Use only
synthetic or appropriately de-identified test data; the CERTH test stack uses
HTTP and does not provide production access controls.

The test interface accepts one tab-separated assay of at most 500 MiB
(524,288,000 bytes), runs one background QC job at a time, and offers the HTML
report, the QCed assay when available, and a ZIP of the available diagnostics.
Input-format failures and validate failures are displayed; a failed QC run can
still have diagnostic downloads. A thrown error does not reuse an earlier run.
Advanced settings, feature linking and custom rule uploads are not exposed.

The following assumes OmicsQC and `CVDLink` are sibling directories and that the
OmicsQC proxy configuration and tool registration have been added to the CERTH
test stack:

```sh
# From the OmicsQC repository
docker network inspect certh-cvdlink-test-net >/dev/null 2>&1 || \
  docker network create certh-cvdlink-test-net

cd ../CVDLink/certh-cvdlink
docker compose up --build -d

cd ../../OmicsQC
docker compose -f compose.cvdlink.yaml up --build -d
docker exec certh-cvdlink-rproxy-integration-test sh -c 'nginx -t && nginx -s reload'
```

Create a test account and log in at `http://localhost`, then open
`http://localhost/omicsqc/` manually. Tool cards are not implemented in the
test frontend. `http://localhost/api/v1/external-tools/all` should list
`omicsqc` in its `data` array. The app provides a small synthetic example for
testing.

The app's internal authentication endpoint is configured by `OMICSQC_AUTH_URL`
in `compose.cvdlink.yaml`; its default targets the CERTH test backend. The app
checks the captured browser cookie at Shiny session startup and redirects
unauthenticated sessions to `/login`. The test backend does not revoke a JWT
when logging out, and an already-open Shiny session is not immediately
revoked or continuously reauthenticated. Downloads additionally require a
session-specific browser cookie and a current successful authentication check;
a copied download URL alone is insufficient. Registry roles are not enforced.
Closing and reopening the app starts a new authentication check. Do not use
this PoC for sensitive or production data.

Uploads and outputs use a private directory per Shiny session. Closing the
session cancels its worker and removes session files; the container's temporary
storage is a 4 GiB tmpfs that is cleared when the container stops. This tmpfs
consumes host memory as used, in addition to the memory required by R. The
500 MiB input cap is not a promise that a large or wide dataset will fit in
available memory or run quickly. If the single worker is busy, another run
must be retried after it finishes.

To stop only the OmicsQC service, from the OmicsQC repository:

```sh
docker compose -f compose.cvdlink.yaml down
```

Do not remove the shared external network or reset the CERTH database as part
of normal OmicsQC cleanup. See the CERTH integration guide for the proxy and
registry formats and how to operate the rest of the test stack.

### Development checks

The package tests use `testthat`. The Shiny integration tests additionally
need `shiny`, `httr2`, `processx`, `zip`, `openssl`, and Node.js for the local
authentication stub. With the OmicsQC package installed, run:

```sh
OMICSQC_SOURCE="$PWD" Rscript -e \
  'testthat::test_file("shiny/tests/test-integration.R", reporter = "summary", stop_on_failure = TRUE)'
node --test shiny/tests/test-upload-limit.cjs
```

The container installs its runtime dependencies from the dated CRAN snapshot
specified in the Dockerfile. The Shiny tests exercise the real background
worker with synthetic fixtures, upload-size boundaries, authentication,
download isolation, error handling and session cleanup.

## Contact

Urmo Võsa (urmo.vosa at gmail.com), Andres Veidenberg (andres.veidenberg at gmail.com)
