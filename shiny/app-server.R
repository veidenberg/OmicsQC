make_omicsqc_server <- function(
    auth_url, worker_script, manager, authenticate = check_auth,
    register_download = function(session, name, filter) session$registerDataObj(name, NULL, filter)) {
  function(input, output, session) {
    if (!authenticate(session$request$HTTP_COOKIE, auth_url)) {
      session$sendCustomMessage("omicsqc-redirect", "/login")
      return(invisible())
    }

    session_dir <- tempfile("omicsqc-session-", tmpdir = tempdir())
    if (!dir.create(session_dir, mode = "0700")) {
      stop("Could not create temporary session storage.", call. = FALSE)
    }
    state <- shiny::reactiveValues(
      upload_ready = FALSE,
      upload_nonce = NULL,
      upload_message = "Select a TSV assay to begin.",
      worker = NULL,
      run_dir = NULL,
      result = NULL,
      error = NULL
    )
    cookie_name <- paste0("omicsqc_", session$token)
    cookie_secret <- paste(sprintf("%02x", as.integer(openssl::rand_bytes(32))), collapse = "")
    session$sendCustomMessage("omicsqc-session-cookie", list(name = cookie_name, value = cookie_secret))

    clear_result <- function() {
      state$result <- NULL
      state$error <- NULL
      if (!is.null(state$run_dir)) {
        remove_session_files(state$run_dir)
        state$run_dir <- NULL
      }
    }

    shiny::observeEvent(input$assay_selection, {
      if (!is.null(state$worker)) {
        return()
      }
      clear_result()
      state$upload_ready <- FALSE
      state$upload_nonce <- input$assay_selection$nonce
      if (isTRUE(input$assay_selection$oversize)) {
        state$upload_message <- "The assay exceeds the 500 MiB upload limit."
      } else {
        state$upload_message <- "Uploading the assay. If the upload fails, select the file again."
      }
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$assay, {
      if (!is.null(state$worker)) {
        return()
      }
      clear_result()
      tryCatch(
        {
          check_assay_upload(input$assay)
          state$upload_ready <- TRUE
          state$upload_nonce <- input$assay_selection$nonce
          state$upload_message <- "Assay uploaded. Ready to run default QC."
        },
        error = function(error) {
          state$upload_ready <- FALSE
          state$upload_message <- conditionMessage(error)
        }
      )
    }, ignoreNULL = TRUE)

    shiny::observeEvent(input$run, {
      if (!is.null(state$worker)) {
        return()
      }
      clear_result()
      if (!state$upload_ready) {
        state$error <- "Select and upload an assay before running QC."
        return()
      }
      if (!is.null(manager$owner)) {
        state$error <- "Another QC job is running. Please retry when it finishes."
        return()
      }

      upload_path <- tryCatch(
        check_assay_upload(input$assay),
        error = function(error) {
          state$error <- conditionMessage(error)
          NULL
        }
      )
      if (is.null(upload_path)) {
        return()
      }
      run_dir <- tempfile("run-", tmpdir = session_dir)
      if (!dir.create(run_dir, mode = "0700")) {
        state$error <- "Could not create temporary storage for this run."
        return()
      }
      state$run_dir <- run_dir
      manager$owner <- session$token
      worker <- tryCatch(
        processx::process$new(
          file.path(R.home("bin"), "Rscript"),
          c("--vanilla", worker_script, upload_path, run_dir),
          stdout = nullfile(),
          stderr = nullfile(),
          cleanup = TRUE
        ),
        error = function(error) {
          message("Could not start the OmicsQC worker: ", conditionMessage(error))
          NULL
        }
      )
      if (is.null(worker)) {
        manager$owner <- NULL
        state$error <- "Could not start the QC worker."
        return()
      }
      state$worker <- worker
      state$upload_message <- "QC is running in the background."
      session$sendCustomMessage("omicsqc-disable-upload", TRUE)
    })

    shiny::observe({
      shiny::invalidateLater(500, session)
      worker <- state$worker
      if (is.null(worker) || worker$is_alive()) {
        return()
      }
      exit_status <- worker$get_exit_status()
      manifest_path <- file.path(state$run_dir, "manifest.rds")
      manifest <- tryCatch(
        if (file.exists(manifest_path)) readRDS(manifest_path) else NULL,
        error = function(error) {
          message("Could not read the QC worker result: ", conditionMessage(error))
          NULL
        }
      )
      state$worker <- NULL
      if (identical(manager$owner, session$token)) {
        manager$owner <- NULL
      }
      session$sendCustomMessage("omicsqc-disable-upload", FALSE)
      state$upload_message <- "Assay uploaded. Ready to run default QC."
      if (is.null(manifest)) {
        state$error <- paste0("The QC worker exited without a result (exit status ", exit_status, ").")
      } else if (identical(manifest$state, "error")) {
        state$error <- paste("QC could not be completed:", manifest$message)
      } else if (!identical(exit_status, 0L) ||
                 !identical(manifest$state, "completed") ||
                 !(manifest$status %in% c("pass", "fail"))) {
        state$error <- "The QC worker returned an invalid result."
      } else {
        state$result <- manifest
      }
    })

    shiny::observe({
      session$sendCustomMessage(
        "omicsqc-run-button",
        list(enabled = state$upload_ready && is.null(state$worker), nonce = state$upload_nonce)
      )
    })
    output$upload_status <- shiny::renderText(state$upload_message)
    output$run_status <- shiny::renderText({
      if (!is.null(state$error)) {
        state$error
      } else if (!is.null(state$result)) {
        paste("QC status:", toupper(state$result$status))
      } else if (!is.null(state$worker)) {
        "Running QC. Please keep this page open."
      } else {
        ""
      }
    })

    output$result_ui <- shiny::renderUI({
      result <- state$result
      shiny::req(result)
      shiny::tagList(
        shiny::h3("Dataset summary"),
        shiny::tableOutput("dataset_qc"),
        shiny::h3("Built-in checks and warnings"),
        shiny::tableOutput("built_in_checks"),
        shiny::h3("Validation results"),
        shiny::tableOutput("validate_results"),
        if (nrow(result$validate_errors)) {
          shiny::tagList(shiny::h3("Validation configuration errors"), shiny::tableOutput("validate_errors"))
        },
        shiny::h3(paste("Removed samples:", nrow(result$removed_samples))),
        shiny::tableOutput("removed_samples"),
        shiny::h3(paste("Removed features:", nrow(result$removed_features))),
        shiny::tableOutput("removed_features"),
        shiny::p("Tables show up to 100 rows; complete diagnostics are in the output archive."),
        download_link(report_url, "Download HTML report", "qc_report.html"),
        if ("assay_qc_harmonized.tsv" %in% result$files) {
          download_link(assay_url, "Download QCed assay", "assay_qc_harmonized.tsv")
        },
        download_link(archive_url, "Download all available outputs", "omicsqc_outputs.zip")
      )
    })
    output$dataset_qc <- shiny::renderTable(utils::head(state$result$dataset_qc, 100))
    output$built_in_checks <- shiny::renderTable(utils::head(state$result$built_in_checks, 100))
    output$validate_results <- shiny::renderTable(utils::head(state$result$validate_results, 100))
    output$validate_errors <- shiny::renderTable(utils::head(state$result$validate_errors, 100))
    output$removed_samples <- shiny::renderTable(utils::head(state$result$removed_samples, 100))
    output$removed_features <- shiny::renderTable(utils::head(state$result$removed_features, 100))

    output_file <- function(name) {
      if (is.null(state$result) || !(name %in% state$result$files)) {
        stop("The requested output is not available.", call. = FALSE)
      }
      path <- file.path(state$run_dir, "outputs", name)
      if (!file.exists(path)) {
        stop("The requested output is no longer available.", call. = FALSE)
      }
      path
    }
    download_link <- function(url, label, filename) {
      shiny::tags$a(
        href = url, class = "btn btn-default", role = "button",
        download = filename, label
      )
    }
    make_download <- function(name, filename, content_type) {
      register_download(session, name, function(data, request) {
        if (!session_cookie_matches(request$HTTP_COOKIE, cookie_name, cookie_secret)) {
          return(shiny::httpResponse(403L, "text/plain", "This download belongs to another session."))
        }
        if (!authenticate(request$HTTP_COOKIE, auth_url)) {
          return(shiny::httpResponse(401L, "text/plain", "Authentication is required."))
        }
        tryCatch(
          shiny::isolate({
            owned <- FALSE
            if (identical(name, "example")) {
              path <- system.file("extdata", "example_assay.tsv", package = "OmicsQC")
              if (!nzchar(path) || !file.exists(path)) {
                stop("The synthetic example is unavailable.", call. = FALSE)
              }
            } else if (identical(name, "archive")) {
              if (is.null(state$result)) {
                stop("The requested outputs are not available.", call. = FALSE)
              }
              path <- file.path(state$run_dir, "omicsqc_outputs.zip")
              if (!file.exists(path)) {
                stop("The output archive is no longer available.", call. = FALSE)
              }
            } else {
              path <- output_file(filename)
            }
            shiny::httpResponse(
              200L, content_type, list(file = path, owned = owned),
              headers = list(
                "Content-Disposition" = paste0('attachment; filename="', filename, '"'),
                "Cache-Control" = "no-store"
              )
            )
          }),
          error = function(error) {
            message(
              "OmicsQC download failed: ",
              gsub(session_dir, "[temporary directory]", conditionMessage(error), fixed = TRUE)
            )
            shiny::httpResponse(500L, "text/plain", "The requested download could not be created.")
          }
        )
      })
    }
    example_url <- make_download("example", "example_assay.tsv", "text/tab-separated-values")
    report_url <- make_download("report", "qc_report.html", "text/html")
    assay_url <- make_download("assay", "assay_qc_harmonized.tsv", "text/tab-separated-values")
    archive_url <- make_download("archive", "omicsqc_outputs.zip", "application/zip")
    output$example_link <- shiny::renderUI(download_link(example_url, "Download synthetic example", "example_assay.tsv"))

    session$onSessionEnded(function() {
      worker <- shiny::isolate(state$worker)
      if (!is.null(worker) && worker$is_alive()) {
        worker$kill()
        worker$wait(5000)
        if (worker$is_alive()) {
          message("Could not stop the OmicsQC worker; temporary files were retained.")
          return()
        }
      }
      if (identical(manager$owner, session$token)) {
        manager$owner <- NULL
      }
      remove_session_files(session_dir)
    })
  }
}
