MAX_ASSAY_BYTES <- 500 * 1024^2

check_auth <- function(cookie, auth_url) {
  if (is.null(cookie) || !nzchar(cookie)) {
    return(FALSE)
  }

  response <- tryCatch(
    httr2::request(auth_url) |>
      httr2::req_headers(Cookie = cookie) |>
      httr2::req_error(is_error = function(response) FALSE) |>
      httr2::req_timeout(5) |>
      httr2::req_perform(),
    error = function(error) {
      message("Authentication service unavailable: ", conditionMessage(error))
      NULL
    }
  )
  if (is.null(response)) {
    return(FALSE)
  }
  if (httr2::resp_status(response) != 200L) {
    message("Authentication rejected with HTTP ", httr2::resp_status(response))
    return(FALSE)
  }
  TRUE
}

new_job_manager <- function() {
  manager <- new.env(parent = emptyenv())
  manager$owner <- NULL
  manager
}

session_cookie_matches <- function(header, name, secret) {
  if (is.null(header) || !nzchar(header)) {
    return(FALSE)
  }
  cookies <- trimws(strsplit(header, ";", fixed = TRUE)[[1]])
  prefix <- paste0(name, "=")
  matches <- cookies[startsWith(cookies, prefix)]
  length(matches) == 1L && identical(substring(matches[[1]], nchar(prefix) + 1L), secret)
}

check_assay_upload <- function(upload) {
  if (is.null(upload) || nrow(upload) != 1L || !file.exists(upload$datapath[[1]])) {
    stop("Select an assay file before running QC.", call. = FALSE)
  }
  size <- file.info(upload$datapath[[1]])$size
  if (is.na(size)) {
    stop("The uploaded assay could not be read.", call. = FALSE)
  }
  if (size > MAX_ASSAY_BYTES) {
    stop("The assay exceeds the 500 MiB upload limit.", call. = FALSE)
  }
  invisible(upload$datapath[[1]])
}

remove_session_files <- function(path) {
  if (dir.exists(path)) {
    unlink(path, recursive = TRUE)
    if (dir.exists(path)) {
      message("Could not remove temporary OmicsQC session files.")
    }
  }
}
