# Rscript --vanilla tests/action-errors.R
# Exercises the shared guards and actual button/upload/download registrations.
suppressPackageStartupMessages({ library(shiny); library(dplyr); library(tidyr); library(writexl) })
is_assignment <- function(x, name) is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(name))
server_code <- Filter(function(x) is_assignment(x, "server"), as.list(parse("server.R")))[[1]]
statements <- as.list(server_code[[3]][[3]])[-1]
assignment <- function(name) Filter(function(x) is_assignment(x, name), statements)[[1]]
event <- function(name) Filter(function(x) is.call(x) && identical(x[[1]], as.name("observe_action")) &&
  identical(x[[2]], as.call(list(as.name("$"), as.name("input"), as.name(name)))), statements)[[1]]
guard_names <- c("action_label", "report_action_error", "with_action_error", "observe_action", "observe_updates")
walk <- function(expr, predicate) {
  if (missing(expr)) return(list())
  result <- if (predicate(expr)) list(expr) else list()
  if (is.call(expr) || is.expression(expr)) for (child in as.list(expr)[-1]) result <- c(result, walk(child, predicate))
  result
}
bare_observer <- function(x) is.call(x) && is.symbol(x[[1]]) && as.character(x[[1]]) %in% c("observe", "observeEvent")
stopifnot(!length(walk(server_code, bare_observer)))
registered <- walk(server_code, function(x) is.call(x) && identical(x[[1]], as.name("observe_action")))
buttons <- c(walk(parse("ui.R"), function(x) is.call(x) && identical(x[[1]], as.name("actionButton"))),
             walk(server_code, function(x) is.call(x) && identical(x[[1]], as.name("actionButton"))))
button_ids <- unique(vapply(buttons, function(x) as.character(x[[2]]), character(1)))
event_ids <- unique(unlist(lapply(registered, function(x) {
  e <- x[[2]]
  if (is.call(e) && identical(e[[1]], as.name("$")) && identical(e[[2]], as.name("input"))) as.character(e[[3]])
})))
stopifnot(!length(setdiff(button_ids, event_ids)))

signals <- new.env(); signals$errors <- list(); signals$messages <- list(); signals$notifications <- list()
testServer(function(input, output, session) {
  for (name in guard_names) eval(assignment(name), envir = environment())
  original_report <- report_action_error
  report_action_error <- function(action, error) {
    signals$errors[[length(signals$errors) + 1L]] <- list(action = action, message = conditionMessage(error))
    original_report(action, error)
  }
  session$sendCustomMessage <- function(type, message) signals$messages[[length(signals$messages) + 1L]] <- list(type = type, message = message)
  session$sendNotification <- function(type, message) signals$notifications[[length(signals$notifications) + 1L]] <- list(type = type, message = message)
  state <- reactiveValues(value = "kept", count = 0L, cleanups = 0L, updates = 0L)
  observe_action(input$test_failure, {
    if (input$fail) stop("invalid file fixture")
    state$count <- state$count + 1L
    input$unrelated  # observeEvent must keep handler reads isolated.
  }, ignoreInit = TRUE, on_error = function() state$cleanups <- state$cleanups + 1L)
  observe_action(input$test_silent, { req(input$allowed); state$count <- state$count + 1L }, ignoreInit = TRUE)
  observe_action(input$test_validation, { shiny::validate(shiny::need(FALSE, "Select a file")) }, ignoreInit = TRUE)
  observe_action(input$test_return, {
    local_value <- "captured"
    output$captured <- renderText(local_value)
    if (input$early) return(NULL)
    state$count <- state$count + 1L
  }, ignoreInit = TRUE)
  observe_action(input$test_once, { state$count <- state$count + 1L }, once = TRUE, ignoreInit = TRUE)
  observe_updates({
    req(input$refresh)
    if (input$refresh == "fail") stop("refresh fixture")
    state$updates <- isolate(state$updates) + 1L
  })
}, {
  session$flushReact()
  zero <- structure(0L, class = "shinyActionButtonValue")
  session$setInputs(fail = TRUE, early = FALSE, allowed = FALSE, test_failure = zero,
                   test_silent = zero, test_return = zero, test_once = zero, test_validation = zero)
  session$setInputs(test_failure = 1)
  stopifnot(!session$isClosed(), state$value == "kept", state$count == 0L, state$cleanups == 1L,
            length(signals$errors) == 1L)
  notice <- tail(signals$notifications, 1)[[1]]$message
  stopifnot(notice$type == "error", is.null(notice$duration),
            grepl("invalid file fixture", notice$html), any(vapply(signals$messages, function(x) x$type == "waiter-hide", logical(1))))
  session$setInputs(fail = FALSE, test_failure = 2)
  stopifnot(!session$isClosed(), state$count == 1L)
  session$setInputs(unrelated = "changed")
  stopifnot(state$count == 1L)  # no accidental dependency on handler inputs
  before <- length(signals$errors)
  session$setInputs(allowed = FALSE, test_silent = 1, test_validation = 1)
  stopifnot(!session$isClosed(), length(signals$errors) == before)
  session$setInputs(allowed = TRUE, test_silent = 2)
  stopifnot(state$count == 2L)
  session$setInputs(early = TRUE, test_return = 1)
  stopifnot(output$captured == "captured", state$count == 2L)
  session$setInputs(early = FALSE, test_return = 2)
  stopifnot(state$count == 3L)
  session$setInputs(test_once = 1); session$setInputs(test_once = 2)
  stopifnot(state$count == 4L)
  session$setInputs(refresh = "fail")
  stopifnot(!session$isClosed(), tail(signals$errors, 1)[[1]]$action == "Refresh interface")
  session$setInputs(refresh = "ok")
  stopifnot(state$updates == 1L, !session$isClosed())
})

# Actual button bodies: metadata fails without replacing current data, then can
# be retried; batch processing and uploads fail without disconnecting the session.
work <- tempfile("action-errors-"); dir.create(file.path(work, "Data"), recursive = TRUE)
write.csv(data.frame(WrongHeader = "Run01", `1` = "Patient1", check.names = FALSE), file.path(work, "Data", "metadata.csv"), row.names = FALSE)
raw <- file.path(work, "fixture.mz5"); writeLines("synthetic raw input", raw)
signals$errors <- list()
testServer(function(input, output, session) {
  for (name in c(guard_names, "metadata_file_paths", "connect_peak_area_metadata", "apply_batch_correction", "new_run_timer", "run_new_analysis", "valid_mz5_names")) eval(assignment(name), envir = environment())
  original_report <- report_action_error
  report_action_error <- function(action, error) {
    signals$errors[[length(signals$errors) + 1L]] <- list(action = action, message = conditionMessage(error))
    original_report(action, error)
  }
  original_area <- data.frame(file.name = "Run01", peak.number = 1, `100_A` = 10, `200_B` = 20, check.names = FALSE)
  peak_areas_data <- reactiveVal(original_area)
  main_folder2 <- function() work
  batch_order <- function() data.frame(SCode = "ref", Batch = 1, Seq = 1)
  doBC <- function(...) stop("batch model fixture")
  adherence_patient_data <- reactiveVal("previous upload")
  manual_peak_integration <- function(...) stop("invalid integration fixture")
  regenerate_metabolite_peak_areas <- function() stop("must not continue after integration failure")
  analysis_mode <- function() "new"
  uploadedmz5 <- function() data.frame(FileName = "Run01.mz5", FilePath = raw, Source = "upload")
  create_results_folder <- function(...) work
  prepare_analysis_parameters <- function(...) list(mi_df = data.frame())
  prepare_results_folder <- function(...) invisible(NULL)
  massData <- refMassListData <- parametersData <- userMTI <- function() data.frame()
  process_mz5_file <- function(raw_path, filename, file_name, config, timing) {
    withProgress(message = "Reading fixture", value = 0, {
      timing$stage("Read mz5 into memory", filename)
      stop("raw reading fixture")
    })
  }
  for (id in c("connect_metadata", "batch_corr", "adherence_results", "confirm_peak_number", "initialize")) eval(event(id), envir = environment())
}, {
  session$setInputs(connect_metadata = 1)
  stopifnot(!session$isClosed(), identical(peak_areas_data(), original_area),
            tail(signals$errors, 1)[[1]]$action == "Connect Metadata to Results")
  write.csv(data.frame(`Data file` = "Run01", `1` = "Patient1", check.names = FALSE), file.path(work, "Data", "metadata.csv"), row.names = FALSE)
  session$setInputs(connect_metadata = 2)
  stopifnot(!session$isClosed(), peak_areas_data()$`PBM Sample ID` == "Patient1",
            file.exists(file.path(work, "merged_peak_area.csv")))
  session$setInputs(batch_corr = 1)
  stopifnot(!session$isClosed(), tail(signals$errors, 1)[[1]]$action == "Perform Batch Correction")
  session$setInputs(adherence_results = list(datapath = file.path(work, "missing.xlsx")))
  stopifnot(!session$isClosed(), adherence_patient_data() == "previous upload",
            tail(signals$errors, 1)[[1]]$action == "Load adherence results")
  session$setInputs(manual_peak_number = "bad", confirm_peak_number = 1)
  stopifnot(!session$isClosed(), tail(signals$errors, 1)[[1]]$message == "invalid integration fixture")
  session$setInputs(initialize = 1)
  stopifnot(!session$isClosed(), tail(signals$errors, 1)[[1]]$action == "Initialize Run")
  timings <- read.csv(file.path(work, "run_timings.csv"))
  stopifnot(tail(timings$status, 1) == "failed", any(timings$stage == "Read mz5 into memory" & timings$status == "failed"))
})

# Use the actual download content callback; a failed write must remove its
# partial temporary file rather than returning a corrupt download.
download <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(paste(deparse(x[[2]]), collapse = ""), "output$downloadData"), statements)[[1]][[3]]
testServer(function(input, output, session) {
  for (name in guard_names) eval(assignment(name), envir = environment())
  parametersData <- massData <- refMassListData <- function() data.frame(value = 1)
  fail_write <- reactiveVal(TRUE)
  write_xlsx <- function(x, path) {
    writeLines("partial", path)
    if (fail_write()) stop("disk write fixture")
    writexl::write_xlsx(x, path)
  }
  download_content <- eval(download$content, envir = environment())
}, {
  target <- file.path(work, "download.xlsx")
  download_content(target)
  stopifnot(!session$isClosed(), !file.exists(target))
  fail_write(FALSE)
  download_content(target)
  stopifnot(!session$isClosed(), file.exists(target), file.info(target)$size > 100)
})
unlink(work, recursive = TRUE)
cat("PASS: all registered buttons use action guards; ordinary errors show notifications and release waiters; silent validation/returns/isolation/once are preserved; real metadata, batch, upload, integration, raw-run timing and download failures keep the session alive and permit retry.\n")
