library(shiny)

source("helpers.R", local = TRUE)
source("app-server.R", local = TRUE)

options(shiny.maxRequestSize = 501 * 1024^2)
auth_url <- Sys.getenv("OMICSQC_AUTH_URL", unset = "")
if (!nzchar(auth_url)) {
  stop("OMICSQC_AUTH_URL must be configured.", call. = FALSE)
}

ui <- fluidPage(
  class = "omics-app",
  title = "OmicsQC",
  tags$head(
    tags$script(src = "embedding.js"),
    tags$link(rel = "stylesheet", type = "text/css", href = "omicsqc.css"),
    tags$script(src = "upload-limit.js")
  ),
  tags$header(
    class = "omics-nav",
    tags$a(class = "omics-brand", href = "/data-qc", "CVDLINK"),
    tags$span(class = "omics-nav-section", "OmicsQC"),
    tags$a(class = "omics-back-link", href = "/data-qc", "Back to Data QC")
  ),
  tags$main(
    class = "omics-main",
    tags$div(
      class = "omics-shell",
      tags$header(
        class = "omics-header",
        tags$div(
          tags$p(class = "omics-eyebrow", "CVDLINK / QUALITY CONTROL"),
          tags$h1("OmicsQC"),
          tags$p(
            class = "omics-subtitle",
            "Run the default quality-control pipeline and review the results in one workspace."
          )
        ),
        tags$span(class = "omics-environment", "Quality control")
      ),
      tags$div(
        class = "omics-content",
        tags$div(
          class = "omics-notice",
          tags$strong("Test data only"),
          tags$span("Use synthetic or de-identified data. Results are not clinical decisions.")
        ),
        tags$div(
          class = "omics-workspace",
          tags$section(
            class = "omics-card omics-input-card",
            `aria-labelledby` = "assay-heading",
            tags$p(class = "omics-eyebrow", "GET STARTED"),
            tags$h2(id = "assay-heading", "Upload an assay"),
            tags$p(class = "omics-card-intro", "Choose a tab-separated assay to run default QC."),
            tags$div(
              class = "omics-upload",
              fileInput("assay", "Assay TSV", accept = c(".tsv", ".txt")),
              tags$p(class = "omics-help", "TSV or TXT, up to 500 MiB.")
            ),
            tags$div(role = "status", `aria-live` = "polite", textOutput("upload_status")),
            tags$div(
              class = "omics-actions",
              actionButton("run", "Run default QC", disabled = "disabled"),
              uiOutput("example_link", class = "omics-example")
            ),
            tags$div(role = "status", `aria-live` = "polite", textOutput("run_status"))
          ),
          tags$aside(
            class = "omics-card omics-guidance",
            `aria-labelledby` = "requirements-heading",
            tags$p(class = "omics-eyebrow", "INPUT GUIDANCE"),
            tags$h2(id = "requirements-heading", "Before you begin"),
            tags$ul(
              tags$li("Use a sample-by-feature matrix with sample identifiers in the first column."),
              tags$li("All remaining columns must contain numeric feature measurements."),
              tags$li("Inputs must already be normalized, transformed and corrected."),
              tags$li("Zero is treated as below detection by default.")
            ),
            tags$p(class = "omics-guidance-note", "Only one QC run can execute at a time.")
          )
        ),
        uiOutput("result_ui")
      )
    )
  )
)

shinyApp(
  ui = ui,
  server = make_omicsqc_server(
    auth_url = auth_url,
    worker_script = normalizePath("worker.R", mustWork = TRUE),
    manager = new_job_manager()
  )
)
