library(testthat)
library(shiny)
library(OmicsQC)

root <- normalizePath(Sys.getenv("OMICSQC_SOURCE"), mustWork = TRUE)
source(file.path(root, "shiny", "helpers.R"))
source(file.path(root, "shiny", "app-server.R"))
source(file.path(root, "shiny", "worker.R"))

fixture <- function(name) file.path(root, "tests", "testthat", "fixtures", name)
test_downloads <- new.env(parent = emptyenv())
register_test_download <- function(session, name, filter) {
  test_downloads[[name]] <- filter
  paste0("session/test/dataobj/", name)
}
uploaded_file <- function(path) {
  data.frame(
    name = basename(path),
    size = file.info(path)$size,
    type = "text/tab-separated-values",
    datapath = path,
    stringsAsFactors = FALSE
  )
}

wait_for_worker <- function(session, state) {
  for (attempt in seq_len(100)) {
    Sys.sleep(0.1)
    session$elapse(500)
    if (is.null(state$worker)) {
      return(invisible())
    }
  }
  stop("The QC worker did not finish during the test.")
}

test_that("the Shiny app directory serves the OmicsQC UI", {
  previous <- Sys.getenv("OMICSQC_AUTH_URL", unset = NA_character_)
  on.exit(if (is.na(previous)) Sys.unsetenv("OMICSQC_AUTH_URL") else Sys.setenv(OMICSQC_AUTH_URL = previous))
  Sys.setenv(OMICSQC_AUTH_URL = "http://127.0.0.1:12345/api/v1/auth/me")
  app <- shiny::shinyAppDir(file.path(root, "shiny"))
  response <- app$httpHandler(list(PATH_INFO = "/", REQUEST_METHOD = "GET"))
  expect_equal(response$status, 200)
  expect_match(response$content, "OmicsQC integration test", fixed = TRUE)
  expect_false(grepl("No UI defined", response$content, fixed = TRUE))
})

test_that("authentication forwards cookies and fails closed", {
  port <- httpuv::randomPort()
  script <- paste0(
    "require('http').createServer((req,res)=>{",
    "res.statusCode=req.headers.cookie==='token=valid'?200:401;",
    "res.end('{}');}).listen(", port, ", '127.0.0.1',()=>console.log('ready'));"
  )
  backend <- processx::process$new("node", c("-e", script), stdout = "|", stderr = "|", cleanup = TRUE)
  on.exit(backend$kill())
  for (attempt in seq_len(100)) {
    if (length(backend$read_output_lines())) break
    Sys.sleep(0.05)
  }
  expect_true(backend$is_alive())
  url <- paste0("http://127.0.0.1:", port, "/api/v1/auth/me")
  expect_true(check_auth("token=valid", url))
  expect_false(check_auth("token=invalid", url))
  expect_false(check_auth(NULL, url))
  backend$kill()
  expect_false(check_auth("token=valid", url))
})

test_that("the worker distinguishes QC failure from execution error", {
  run_dir <- tempfile("worker-test-")
  dir.create(run_dir)
  on.exit(unlink(run_dir, recursive = TRUE))

  expect_identical(run_worker(fixture("duplicate_sample_assay.tsv"), run_dir), 0L)
  result <- readRDS(file.path(run_dir, "manifest.rds"))
  expect_identical(result$state, "completed")
  expect_identical(result$status, "fail")
  expect_true("built_in_checks.tsv" %in% result$files)
  expect_false("assay_qc_harmonized.tsv" %in% result$files)

  bad_input <- file.path(run_dir, "invalid.tsv")
  writeLines("sample_id", bad_input)
  expect_identical(run_worker(bad_input, run_dir), 1L)
  error <- readRDS(file.path(run_dir, "manifest.rds"))
  expect_identical(error$state, "error")
  expect_match(error$message, "at least one feature")
  expect_false(grepl(run_dir, error$message, fixed = TRUE))
})

test_that("the exact 500 MiB server-side file limit is enforced", {
  path <- tempfile("upload-size-")
  on.exit(unlink(path))
  connection <- file(path, "wb")
  seek(connection, where = MAX_ASSAY_BYTES - 1, origin = "start")
  writeBin(as.raw(0), connection)
  close(connection)

  expect_equal(file.info(path)$size, MAX_ASSAY_BYTES)
  expect_invisible(check_assay_upload(uploaded_file(path)))
  connection <- file(path, "ab")
  writeBin(as.raw(0), connection)
  close(connection)
  expect_equal(file.info(path)$size, MAX_ASSAY_BYTES + 1)
  expect_error(check_assay_upload(uploaded_file(path)), "500 MiB")
})

test_that("download cookies must match the current session exactly", {
  expect_true(session_cookie_matches("token=valid; omicsqc_id=secret", "omicsqc_id", "secret"))
  expect_false(session_cookie_matches("omicsqc_id=wrong", "omicsqc_id", "secret"))
  expect_false(session_cookie_matches("omicsqc_id=secret; omicsqc_id=secret", "omicsqc_id", "secret"))
  expect_false(session_cookie_matches(NULL, "omicsqc_id", "secret"))
})

test_that("a Shiny session runs QC, downloads results and clears stale state", {
  manager <- new_job_manager()
  auth_state <- new.env()
  auth_state$ok <- TRUE
  server <- make_omicsqc_server(
    "http://unused",
    file.path(root, "shiny", "worker.R"),
    manager,
    authenticate = function(cookie, auth_url) auth_state$ok,
    register_download = register_test_download
  )
  testServer(server, {
    session$setInputs(assay = uploaded_file(fixture("clean_assay.tsv")))
    expect_true(state$upload_ready)
    session$setInputs(run = 1)
    expect_false(is.null(state$worker))
    wait_for_worker(session, state)
    expect_identical(state$result$status, "pass")
    expect_null(state$error)
    expect_null(manager$owner)

    cookie <- paste0(cookie_name, "=", cookie_secret)
    request <- list(HTTP_COOKIE = cookie)
    expect_equal(test_downloads$report(NULL, list(HTTP_COOKIE = NULL))$status, 403)
    auth_state$ok <- FALSE
    expect_equal(test_downloads$report(NULL, request)$status, 401)
    auth_state$ok <- TRUE
    report_response <- test_downloads$report(NULL, request)
    expect_equal(report_response$status, 200)
    report <- report_response$content$file
    expect_true(file.exists(report))
    expect_match(paste(readLines(report, warn = FALSE), collapse = "\n"), "OmicsQC Report")
    archive_response <- test_downloads$archive(NULL, request)
    expect_equal(archive_response$status, 200)
    archive <- archive_response$content$file
    expect_true(file.exists(archive))
    archive_files <- zip::zip_list(archive)$filename
    expect_true(all(c("qc_report.html", "assay_qc_harmonized.tsv", "validate_config_errors.tsv") %in% archive_files))

    session$setInputs(assay_selection = list(oversize = TRUE, nonce = 1))
    expect_null(state$result)
    expect_false(state$upload_ready)
    expect_match(state$upload_message, "500 MiB")

    session$setInputs(assay = uploaded_file(fixture("duplicate_sample_assay.tsv")))
    session$setInputs(run = 2)
    wait_for_worker(session, state)
    expect_identical(state$result$status, "fail")
    expect_false("assay_qc_harmonized.tsv" %in% state$result$files)
    archive <- test_downloads$archive(NULL, request)$content$file
    expect_true("built_in_checks.tsv" %in% zip::zip_list(archive)$filename)

    invalid <- tempfile("invalid-assay-")
    writeLines("sample_id", invalid)
    on.exit(unlink(invalid), add = TRUE)
    session$setInputs(assay = uploaded_file(invalid))
    session$setInputs(run = 3)
    wait_for_worker(session, state)
    expect_null(state$result)
    expect_match(state$error, "at least one feature")
  })
})

test_that("busy and unauthenticated sessions cannot start a job", {
  manager <- new_job_manager()
  manager$owner <- "another-session"
  server <- make_omicsqc_server(
    "http://unused", file.path(root, "shiny", "worker.R"), manager,
    authenticate = function(cookie, auth_url) TRUE,
    register_download = register_test_download
  )
  testServer(server, {
    session$setInputs(assay = uploaded_file(fixture("clean_assay.tsv")), run = 1)
    expect_null(state$worker)
    expect_match(state$error, "Another QC job")
    expect_identical(manager$owner, "another-session")
  })

  manager <- new_job_manager()
  server <- make_omicsqc_server(
    "http://unused", file.path(root, "shiny", "worker.R"), manager,
    authenticate = function(cookie, auth_url) FALSE,
    register_download = register_test_download
  )
  testServer(server, {
    session$setInputs(assay = uploaded_file(fixture("clean_assay.tsv")), run = 1)
    expect_null(manager$owner)
  })
})

test_that("disconnect cancels the worker and removes session files", {
  script <- tempfile("slow-worker-", fileext = ".R")
  writeLines("Sys.sleep(30)", script)
  on.exit(unlink(script))
  manager <- new_job_manager()
  server <- make_omicsqc_server(
    "http://unused", script, manager,
    authenticate = function(cookie, auth_url) TRUE,
    register_download = register_test_download
  )
  testServer(server, {
    session$setInputs(assay = uploaded_file(fixture("clean_assay.tsv")), run = 1)
    worker <- state$worker
    path <- session_dir
    expect_true(worker$is_alive())
    expect_true(dir.exists(path))
    session$close()
    expect_false(worker$is_alive())
    expect_false(dir.exists(path))
    expect_null(manager$owner)
  })
})
