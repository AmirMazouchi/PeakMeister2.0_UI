# Rscript --vanilla tests/investigate.R
# Exact raw-window extraction, bounded preview, real Plotly objects and Shiny
# actions. The mzR interface is supplied by a deterministic scan fixture here;
# tests/investigate-mz5.R additionally exercises a genuine mz5 file with mzR.
suppressPackageStartupMessages({ library(shiny); library(plotly) })
is_assignment <- function(x, name) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(paste(deparse(x[[2]]), collapse = ""), name)
code <- readLines("server.R", warn = FALSE)
first <- which(grepl("^  #####Matabo3.6 Investigate raw", code))
last <- which(grepl("^  #####Matabo3.6 End Investigate", code))
stopifnot(length(first) == 1L, length(last) == 1L, first < last)
investigate_block <- parse(text = code[(first + 1L):(last - 1L)])
definition <- Filter(function(x) is_assignment(x, "server"), as.list(parse("server.R")))[[1]]
statements <- as.list(definition[[3]][[3]])[-1]
assignment <- function(name) Filter(function(x) is_assignment(x, name), statements)[[1]]
functions <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function")), as.list(investigate_block))
tasks <- new.env(parent = globalenv())
for (expression in functions) eval(expression, tasks)
same <- function(a, b) stopifnot(isTRUE(all.equal(a, b, check.attributes = TRUE)))
error_text <- function(expr) tryCatch({ force(expr); stop("Expected failure was not raised") },
                                     error = function(e) conditionMessage(e))

# Unsorted spectra; MS2 scans; missing windows; finite-value handling; exact
# range endpoints; more than eight scans to exercise bounded batch reading.
fixture_header <- data.frame(msLevel = c(1L, 2L, rep(1L, 10)),
  retentionTime = c(0, 30, seq(60, 600, by = 60)), lowMZ = 99.9, highMZ = 101)
fixture_scans <- lapply(seq_len(nrow(fixture_header)), function(i) {
  cbind(mz = c(101, 100.0008, 100, 99.9, NA, -1),
        intensity = c(i * 100, i * 4, i * 2, i, 999, 999))
})
fixture_scans[[4]] <- cbind(mz = c(99.9, 101), intensity = c(1, 20))
fixture_scans[[5]] <- cbind(mz = c(100, 100.0008), intensity = c(NA, Inf))
signals <- new.env(); signals$opens <- 0L; signals$closes <- 0L
signals$batches <- integer(); signals$errors <- list(); signals$input_updates <- list()
fixture_opener <- function(path) { signals$opens <- signals$opens + 1L; new.env() }
fixture_closer <- function(handle) { signals$closes <- signals$closes + 1L; invisible(NULL) }
fixture_reader <- function(handle, scans) {
  signals$batches <- c(signals$batches, length(scans))
  result <- fixture_scans[scans]
  if (length(scans) == 1L) result[[1]] else result
}
fixture_path <- tempfile(fileext = ".mz5"); writeLines("raw reader fixture", fixture_path)
open_fixture <- function(path = fixture_path, name = "Fixture.mz5", header = fixture_header) {
  tasks$investigate_open_raw(path, name, opener = fixture_opener,
    read_header = function(handle) header, read_peaks = fixture_reader, closer = fixture_closer)
}
raw <- open_fixture()
same(raw$scans, c(1L, 3:12)); same(raw$time_range, c(0, 10)); same(raw$mz_range, c(99.9, 101))
stopifnot(grepl("no MS1", error_text(open_fixture(header = transform(fixture_header, msLevel = 2L)))))
stopifnot(signals$closes == 1L)
fallback_header <- fixture_header[, c("msLevel", "retentionTime")]
fallback <- open_fixture(header = fallback_header)
same(fallback$mz_range, c(99.9, 101)); tasks$investigate_close_raw(fallback)
inconsistent <- open_fixture(header = transform(fixture_header, lowMZ = 500, highMZ = 1000))
same(inconsistent$mz_range, c(99.9, 101))
stopifnot(inconsistent$range_source == "sampled MS1 measurements")
tasks$investigate_close_raw(inconsistent)
same(tasks$investigate_window_limits(100, 10, "ppm"), c(99.999, 100.001))
same(tasks$investigate_window_limits(100, 0.001, "Da"), c(99.999, 100.001))
stopifnot(grepl("positive", error_text(tasks$investigate_window_limits(100, 0, "ppm"))))
slice <- tasks$investigate_extract_slice(raw, 100, 10, "ppm", c(0, 10))
same(slice$data$time, 0:10)
same(slice$data$intensity, c(3, 9, 0, 0, seq(18, 36, by = 3)))
summed <- tasks$investigate_extract_slice(raw, 100, 0.001, "Da", c(0, 1), "sum")
same(summed$data$intensity, c(6, 18))
maximum <- tasks$investigate_extract_slice(raw, 100, 10, "ppm", c(0, 1), "max")
same(maximum$data$intensity, c(4, 12))
missing <- tasks$investigate_extract_slice(raw, 105, 10, "ppm", c(0, 1))
same(missing$data$intensity, c(0, 0))
stopifnot(grepl("No MS1 scans", error_text(tasks$investigate_extract_slice(raw, 100, 10, "ppm", c(50, 51)))))
grid <- tasks$investigate_build_grid(raw, c(99.9, 101), c(0, 10), 100L, 80L)
same(dim(grid$z), c(80L, 100L))
stopifnot(max(signals$batches) <= 8L, grid$scans == 11L,
          max(grid$z) == 1200, grid$z[80, 100] == 1200, grid$z[1, 1] == 1)
stopifnot(grepl("No m/z measurements", error_text(tasks$investigate_build_grid(raw, c(200, 201), c(0, 10)))))
stopifnot(grepl("minimum < maximum", error_text(tasks$investigate_build_grid(raw, c(101, 99), c(0, 10)))))
negative <- raw
negative$read_peaks <- function(handle, scans) {
  result <- lapply(scans, function(index) matrix(c(100, -5), ncol = 2))
  if (length(result) == 1L) result[[1]] else result
}
negative_grid <- tasks$investigate_build_grid(negative, c(99.9, 101), c(0, 10), 100L, 80L)
stopifnot(min(negative_grid$z) == -5)

# Actual Plotly traces, orientation, range controls, logarithmic display and band.
grid$key <- "fixture"; grid$intensity_range <- c(0, 1500)
surface <- plotly_build(tasks$investigate_plot_grid(grid, c(99.999, 100.001)))
stopifnot(surface$x$data[[1]]$type == "surface", surface$x$data[[2]]$type == "mesh3d")
same(surface$x$layout$scene$xaxis$range, c(99.9, 101))
same(surface$x$layout$scene$yaxis$range, c(0, 10))
same(surface$x$layout$scene$zaxis$range, c(0, 1500))
same(dim(surface$x$data[[1]]$z), dim(grid$z))
same(as.numeric(surface$x$data[[1]]$z), as.numeric(grid$z))
logged <- plotly_build(tasks$investigate_plot_grid(grid, logarithmic = TRUE))
same(as.numeric(logged$x$data[[1]]$z), as.numeric(log10(1 + grid$z)))
slice$key <- "fixture"; slice$intensity_range <- c(0, 40)
trace <- plotly_build(tasks$investigate_plot_slice(slice))
stopifnot(trace$x$data[[1]]$type == "scattergl")
same(as.numeric(trace$x$data[[1]]$y), slice$data$intensity)
same(trace$x$layout$yaxis$range, c(0, 40))
tasks$investigate_close_raw(raw)

# Harness runs the whole actual Matabo section, including all registrations.
guard_names <- c("action_label", "report_action_error", "with_action_error", "observe_action", "observe_updates", "plot_event_data")
investigate_harness <- function(input, output, session) {
  uploadedmz5 <- reactiveVal(data.frame(FileName = "Fixture.mz5", FilePath = fixture_path, Source = "upload"))
  massData <- reactiveVal(data.frame(mz = 100)); refMassListData <- reactiveVal(data.frame(mz = 100))
  for (name in guard_names) eval(assignment(name), environment())
  report_action_error <- function(action, error) {
    signals$errors[[length(signals$errors) + 1L]] <- list(action = action, message = conditionMessage(error))
  }
  # Reproduce packages loaded after Shiny masking its progress updater. Calling
  # this S3 generic with a numeric fraction matches the reported runtime error.
  setProgress <- function(value, ...) UseMethod("setProgress")
  withProgress <- function(...) stop("A masked withProgress() was called.")
  eval(investigate_block, environment())
  investigate_open_raw <- open_fixture
  if (inherits(session, "MockShinySession")) {
    old_send <- session$sendInputMessage
    session$sendInputMessage <- function(id, message) {
      signals$input_updates[[id]] <- message
      old_send(id, message)
    }
  }
}
hash_before <- unname(tools::md5sum(fixture_path))
signals$errors <- list()
testServer(investigate_harness, {
  session$setInputs(investigate_source = "engine", investigate_engine_file = fixture_path,
    investigate_mz_min = 99.9, investigate_mz_max = 101,
    investigate_time_min = 0, investigate_time_max = 10,
    investigate_slice_time_min = 0, investigate_slice_time_max = 10,
    investigate_mz_bins = "100", investigate_time_bins = "80",
    investigate_center = 100, investigate_window = 10, investigate_unit = "ppm",
    investigate_auto_z = TRUE, investigate_auto_y = TRUE, investigate_log_z = FALSE,
    investigate_aggregation = "mean")
  session$setInputs(investigate_draw_3d = 1)
  stopifnot(is.null(investigate_grid()), tail(signals$errors, 1)[[1]]$action == "Investigate: draw 3D preview", !session$isClosed())
  session$setInputs(investigate_load = 1)
  stopifnot(investigate_raw()$name == "Fixture.mz5", is.null(investigate_grid()), is.null(investigate_trace()))
  session$setInputs(investigate_draw_3d = 2, investigate_draw_slice = 1)
  stopifnot(investigate_grid()$scans == 11L)
  same(investigate_trace()$data$intensity, c(3, 9, 0, 0, seq(18, 36, by = 3)))
  built <- jsonlite::fromJSON(output$investigate_3d)
  stopifnot(built$x$source == "metabomeister_investigate")
  session$flushReact()
  reads <- length(signals$batches)
  session$setInputs(investigate_auto_z = FALSE, investigate_z_min = 0, investigate_z_max = 500,
    investigate_auto_y = FALSE, investigate_y_min = 0, investigate_y_max = 50,
    investigate_draw_3d = 3, investigate_draw_slice = 2)
  stopifnot(length(signals$batches) == reads)
  same(investigate_grid()$intensity_range, c(0, 500))
  same(investigate_trace()$intensity_range, c(0, 50))
  session$setInputs(`plotly_click-metabomeister_investigate` = '{"x":100.25,"y":1,"z":10,"curveNumber":0}')
  same(as.numeric(signals$input_updates$investigate_center$value), 100.25)
  session$setInputs(investigate_mz_max = 90, investigate_draw_3d = 4)
  stopifnot(is.null(investigate_grid()), !session$isClosed(), grepl("minimum < maximum", tail(signals$errors, 1)[[1]]$message))
  session$setInputs(investigate_mz_max = 101, investigate_draw_3d = 5)
  stopifnot(!is.null(investigate_grid()), !session$isClosed())
  closed <- signals$closes
  session$setInputs(investigate_load = 2)
  stopifnot(signals$closes == closed + 1L, is.null(investigate_grid()), is.null(investigate_trace()))
  session$setInputs(investigate_clear = 1)
  stopifnot(is.null(investigate_raw()), signals$closes == closed + 2L)
  session$setInputs(investigate_load = 3)
  stopifnot(!is.null(investigate_raw()))
})
stopifnot(signals$opens == signals$closes)
same(unname(tools::md5sum(fixture_path)), hash_before)
cat("PASS: ppm/Da, inclusive boundaries, mean/sum/max, MS1-only and finite data; bounded batches/grid; real 3D/log/slice Plotly traces; actual GUI actions with masked progress functions, cache/range changes, click selection, recoverable errors, handle cleanup and unchanged raw input.\n")

# New UI IDs are unique, and the menu places Investigate between Engine and Visualization.
suppressPackageStartupMessages({ library(shinydashboard); library(DT); library(shinyjs); library(waiter) })
ui_code <- as.list(parse("ui.R")); title <- "PeakMeister v2.0"
actual_ui <- eval(Filter(function(x) is_assignment(x, "ui"), ui_code)[[1]])
menus <- Filter(function(tag) !is.null(tag$attribs[["data-value"]]),
                htmltools::tagQuery(actual_ui)$find("a")$selectedTags())
tabs <- vapply(menus, function(tag) tag$attribs[["data-value"]], character(1))
same(unname(tabs[match("engine", tabs) + 0:2]), c("engine", "investigate", "processing"))
walk_tags <- function(x) {
  if (inherits(x, "shiny.tag")) return(c(list(x), unlist(lapply(x$children, walk_tags), recursive = FALSE)))
  if (is.list(x)) return(unlist(lapply(x, walk_tags), recursive = FALSE))
  list()
}
all_ids <- vapply(Filter(function(tag) !is.null(tag$attribs$id), walk_tags(actual_ui)), function(tag) tag$attribs$id, character(1))
investigate_ids <- all_ids[startsWith(all_ids, "investigate_")]
stopifnot(length(investigate_ids) >= 30L, !anyDuplicated(investigate_ids))
cat("PASS: actual Investigate UI, sidebar order and unique control/output IDs.\n")
