library(shiny)

source("helpers.R", local = TRUE)
source("app-server.R", local = TRUE)

options(shiny.maxRequestSize = 501 * 1024^2)
auth_url <- Sys.getenv("OMICSQC_AUTH_URL", unset = "")
if (!nzchar(auth_url)) {
  stop("OMICSQC_AUTH_URL must be configured.", call. = FALSE)
}

ui <- fluidPage(
  tags$head(tags$script(src = "upload-limit.js")),
  titlePanel("OmicsQC integration test"),
  p("For synthetic or de-identified test data only. Results are not clinical decisions."),
  p(paste(
    "Upload a tab-separated sample-by-feature matrix. The first column must contain",
    "sample identifiers; remaining columns must contain numeric feature measurements.",
    "Inputs must already be normalized, transformed and corrected. Zero is treated",
    "as below detection by default."
  )),
  p("Maximum assay size: 500 MiB. Only one QC run can execute at a time."),
  uiOutput("example_link"),
  fileInput("assay", "Assay TSV", accept = c(".tsv", ".txt")),
  textOutput("upload_status"),
  actionButton("run", "Run default QC", disabled = "disabled"),
  p(textOutput("run_status")),
  uiOutput("result_ui")
)

shinyApp(
  ui = ui,
  server = make_omicsqc_server(
    auth_url = auth_url,
    worker_script = normalizePath("worker.R", mustWork = TRUE),
    manager = new_job_manager()
  )
)
