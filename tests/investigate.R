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
open_fixture <- function(path = fixture_path, name = "Fixture.mz5", header = fixture_header, progress = function(fraction) NULL) {
  tasks$investigate_open_raw(path, name, opener = fixture_opener,
    read_header = function(handle) header, read_peaks = fixture_reader, closer = fixture_closer, progress = progress)
}
raw <- open_fixture()
same(raw$scans, c(1L, 3:12)); same(raw$time_range, c(0, 10)); same(raw$mz_range, c(99.9, 101))
stopifnot(grepl("no MS1", error_text(open_fixture(header = transform(fixture_header, msLevel = 2L)))))
stopifnot(signals$closes == 1L)
closed_before <- signals$closes
stopifnot(grepl("injected sampling failure", error_text(tasks$investigate_open_raw(fixture_path,
  "Fixture.mz5", opener = fixture_opener, read_header = function(handle) fixture_header,
  read_peaks = function(handle, scans) stop("injected sampling failure"), closer = fixture_closer))))
stopifnot(signals$closes == closed_before + 1L)
fallback_header <- fixture_header[, c("msLevel", "retentionTime")]
fallback <- open_fixture(header = fallback_header)
same(fallback$mz_range, c(99.9, 101)); tasks$investigate_close_raw(fallback)
inconsistent <- open_fixture(header = transform(fixture_header, lowMZ = 500, highMZ = 1000))
same(inconsistent$mz_range, c(99.9, 101))
stopifnot(inconsistent$range_source == "all valid MS1 measurements")
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
# Native sampling is measured across every valid MS1 spectrum, not inferred
# from header bounds or from the number of display bins.
stopifnot(!raw$sampling$shared, raw$sampling$time_count == 11L,
  raw$sampling$time_unique == 11L, raw$sampling$time_delta == 1)
same(raw$sampling$mz_count_range, c(0L, 4L))
same(raw$sampling$mz_delta, 0.0008)
stopifnot(grepl("positive whole grouping", error_text(tasks$investigate_group_factors(0, 2))),
  grepl("positive whole grouping", error_text(tasks$investigate_group_factors(2.5, 1))))

# Independent dense oracle on a specified alignment grid. Per-scan cell means
# precede grouping; missing positions and empty spectra contribute zeros.
width <- 0.1; edges <- seq(99.9, 101, by = width)
centers <- head(edges, -1) + diff(edges)/2
matrix_oracle <- matrix(0, nrow = 11L, ncol = length(centers))
for (row in seq_along(raw$scans)) {
  scan <- fixture_scans[[raw$scans[row]]]
  for (column in seq_along(centers)) {
    keep <- is.finite(scan[, 1]) & is.finite(scan[, 2]) & scan[, 1] >= edges[column] &
      if (column == length(centers)) scan[, 1] <= edges[column + 1] else scan[, 1] < edges[column + 1]
    if (any(keep)) matrix_oracle[row, column] <- mean(scan[keep, 2])
  }
}
grid <- tasks$investigate_build_grid(raw, c(99.9, 101), c(0, 10), 2, 3, width)
same(dim(grid$z), c(4L, 6L))
for (row in seq_len(nrow(grid$z))) for (column in seq_len(ncol(grid$z))) {
  rows <- seq.int((row - 1)*3 + 1, min(row*3, 11))
  columns <- seq.int((column - 1)*2 + 1, min(column*2, length(centers)))
  same(grid$z[row, column], mean(matrix_oracle[rows, columns]))
}
same(grid$time, c(1, 4, 7, 9.5))
stopifnot(max(signals$batches) <= 8L, grid$scans == 11L, !grid$native)
stopifnot(grepl("No m/z measurements", error_text(tasks$investigate_build_grid(raw, c(200, 201), c(0, 10), alignment_width = 0.1))))
stopifnot(grepl("minimum < maximum", error_text(tasks$investigate_build_grid(raw, c(101, 99), c(0, 10)))))
stopifnot(grepl("too large", error_text(tasks$investigate_build_grid(raw, c(99.9, 101), c(0, 10), alignment_width = 0.000001))))
stopifnot(grepl("positive alignment", error_text(tasks$investigate_mass_axis(raw, c(99.9, 101), 0))))

# Exact medians, including signed signals, odd/even scan counts, sparse zeros,
# and the disk-spill path. All paths agree with independent stats::median.
for (n in 1:12) for (trial in 1:30) {
  set.seed(n*100 + trial)
  observed <- sample(c(-9, -1, 0, 0, 0, 2, 8), n, replace = TRUE)
  same(tasks$investigate_zero_median(observed[observed != 0], n), median(observed))
}
median_ram <- tasks$investigate_median_spectrum(raw, c(99.9, 101), c(0, 10), width)
median_disk <- tasks$investigate_median_spectrum(raw, c(99.9, 101), c(0, 10), width, cache_bytes = 64)
same(median_ram$data$intensity, apply(matrix_oracle, 2, median))
same(median_disk$data, median_ram$data)
stopifnot(median_ram$method == "cached sparse values", median_disk$method == "partitioned sparse values")
median_crop <- tasks$investigate_median_spectrum(raw, c(99.9, 101), c(1, 2), width)
same(median_crop$data$intensity, apply(matrix_oracle[2:3, , drop = FALSE], 2, median))
stopifnot(median_crop$scans == 2L)

# Shared native coordinates need no alignment, including nonuniform axes,
# missing whole scans, duplicate raw points, negative values and short groups.
native_scans <- list(cbind(c(100, 100.1, 101, 103, 103), c(-4, 2, 6, 8, 12)),
  matrix(numeric(), ncol = 2), cbind(c(103, 101, 100.1, 100), c(18, 12, 4, -8)),
  cbind(c(100, 100.1, 101, 103), c(-2, 6, 18, 30)))
native <- raw; native$scans <- 1:4; native$times <- c(0, 0.5, 2, 5); native$time_range <- c(0, 5)
native$read_peaks <- function(handle, scans) { if (length(scans) == 1) native_scans[[scans]] else native_scans[scans] }
native <- tasks$investigate_measure_sampling(native)
stopifnot(native$sampling$shared, native$sampling$time_delta == 0.5)
same(native$sampling$mz_axis, c(100, 100.1, 101, 103)); same(native$sampling$mz_delta, 0.1)
native_matrix <- rbind(c(-4, 2, 6, 10), rep(0, 4), c(-8, 4, 12, 18), c(-2, 6, 18, 30))
native_grid <- tasks$investigate_build_grid(native, c(99, 104), c(0, 5), 3, 3, alignment_width = -1)
same(native_grid$z, rbind(c(mean(native_matrix[1:3, 1:3]), mean(native_matrix[1:3, 4])),
  c(mean(native_matrix[4, 1:3]), native_matrix[4, 4])))
same(native_grid$mz, c(mean(c(100, 100.1, 101)), 103))
same(native_grid$time, c(mean(c(0, 0.5, 2)), 5))
for (budget in c(64, 1e6)) {
  med <- tasks$investigate_median_spectrum(native, c(99, 104), c(0, 5), cache_bytes = budget)
  same(med$data$intensity, apply(native_matrix, 2, median)); stopifnot(med$native)
}
# More than sixteen spill partitions exercises connection eviction/reopening.
wide <- native; wide$scans <- 1:40; wide$times <- 0:39; wide$time_range <- c(0, 39)
wide$read_peaks <- function(handle, scans) {
  answer <- lapply(scans, function(i) cbind(100:129, (1:30)*i))
  if (length(scans) == 1) answer[[1]] else answer
}
wide <- tasks$investigate_measure_sampling(wide)
wide_disk <- tasks$investigate_median_spectrum(wide, c(100, 129), c(0, 39), cache_bytes = 64)
same(wide_disk$data$intensity, (1:30)*20.5)
# One coordinate in either dimension still produces a usable 3D line/marker.
single_grid <- tasks$investigate_build_grid(native, c(99, 104), c(0, 5), 10, 10)
single_grid$key <- "single"; single_grid$intensity_range <- c(-5, 30)
single_plot <- plotly::plotly_build(tasks$investigate_plot_grid(single_grid))
stopifnot(single_plot$x$data[[1]]$type == "scatter3d")
same(as.numeric(single_plot$x$data[[1]]$z), as.numeric(single_grid$z))
# A transient, sparse spectrum has exactly zero medians, without a second pass.
sparse <- native; sparse$read_peaks <- function(handle, scans) {
  answer <- lapply(scans, function(i) cbind(c(100, 100.1, 101, 103), if (i == 1) 1:4 else rep(0, 4)))
  if (length(scans) == 1) answer[[1]] else answer
}
zero_median <- tasks$investigate_median_spectrum(sparse, c(99, 104), c(0, 5), cache_bytes = 64)
stopifnot(all(zero_median$data$intensity == 0), zero_median$method == "implicit-zero counts")
# Temporary partitions clean up on both success and a failed second pass.
before_temp <- list.files(tempdir(), pattern = "^metabomeister-median-", full.names = TRUE)
read_count <- 0L; failing <- raw
failing$read_peaks <- function(handle, scans) {
  read_count <<- read_count + 1L
  if (read_count == 4L) stop("injected second-pass failure")
  fixture_reader(handle, scans)
}
stopifnot(grepl("injected", error_text(tasks$investigate_median_spectrum(failing,
  c(99.9, 101), c(0, 10), width, cache_bytes = 64))))
same(list.files(tempdir(), pattern = "^metabomeister-median-", full.names = TRUE), before_temp)

# Time-window means are calculated per scan before averaging, so differing
# point counts do not give a spectrum extra weight; absent bins contribute zero.
combined <- tasks$investigate_extract_time_slice(raw, 1, c(99.9, 101),
  mode = "combine", window_seconds = 60, mz_bin_width = 0.01)
stopifnot(combined$scan_count == 3L, combined$mode == "combine")
# Independently locate occupied bins from the reported grid centres.
mass_bin <- which.min(abs(combined$data$mz - 100.005))
same(combined$data$intensity[mass_bin], 4)
same(tail(combined$data$intensity, 1), 140)
overlaid <- tasks$investigate_extract_time_slice(raw, 1.2, c(99.9, 101),
  mode = "overlay", window_seconds = 120, n_spectra = 3)
same(vapply(overlaid$spectra, function(s) s$time, numeric(1)), c(1, 2, 0))
same(overlaid$spectra[[1]]$data$mz, sort(fixture_scans[[3]][1:4, 1]))
stopifnot(grepl("No MS1 scans", error_text(tasks$investigate_extract_time_slice(raw, 1.5,
  c(99.9, 101), window_seconds = 1))))
stopifnot(grepl("positive", error_text(tasks$investigate_extract_time_slice(raw, 1,
  c(99.9, 101), window_seconds = 0))))
stopifnot(grepl("whole spectrum count", error_text(tasks$investigate_extract_time_slice(raw, 1,
  c(99.9, 101), mode = "overlay", n_spectra = 4.5))))

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
stopifnot(surface$x$data[[1]]$lighting$ambient >= 0.8,
  any(grepl("dbeefe", unlist(surface$x$data[[1]]$colorscale), ignore.case = TRUE)))
slice$key <- "fixture"; slice$intensity_range <- c(0, 40)
trace <- plotly_build(tasks$investigate_plot_slice(slice))
stopifnot(trace$x$data[[1]]$type == "scattergl")
same(as.numeric(trace$x$data[[1]]$y), slice$data$intensity)
same(trace$x$layout$yaxis$range, c(0, 40))
combined$key <- "combined"; combined$intensity_range <- c(0, 150)
combined_plot <- plotly_build(tasks$investigate_plot_slice(combined))
stopifnot(combined_plot$x$layout$xaxis$title == "m/z")
overlaid$key <- "overlay"; overlaid$intensity_range <- c(0, 400)
overlay_plot <- plotly_build(tasks$investigate_plot_slice(overlaid))
stopifnot(length(overlay_plot$x$data) == 3L,
  length(unique(vapply(overlay_plot$x$data, function(s) s$line$color, character(1)))) == 3L)
same(as.numeric(overlay_plot$x$data[[1]]$x), overlaid$spectra[[1]]$data$mz)
same(as.numeric(overlay_plot$x$data[[1]]$y), overlaid$spectra[[1]]$data$intensity)
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
    investigate_mz_group = 2, investigate_time_group = 3, investigate_alignment_width = 0.1,
    investigate_slice_axis = "mz",
    investigate_center = 100, investigate_window = 10, investigate_unit = "ppm",
    investigate_auto_z = TRUE, investigate_auto_y = TRUE, investigate_log_z = FALSE,
    investigate_aggregation = "mean")
  session$setInputs(investigate_draw_3d = 1)
  stopifnot(is.null(investigate_grid()), tail(signals$errors, 1)[[1]]$action == "Investigate: draw 3D preview", !session$isClosed())
  session$setInputs(investigate_load = 1)
  stopifnot(investigate_raw()$name == "Fixture.mz5", is.null(investigate_grid()), is.null(investigate_trace()), is.null(investigate_median()))
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
  stopifnot(grepl("Mean intensity", output$investigate_grid_status),
    grepl("display positions", output$investigate_group_status),
    grepl("Minimum positive native Δm/z", output$investigate_status),
    grepl("minimum positive Δt", output$investigate_status))
  session$setInputs(investigate_draw_median = 1)
  same(investigate_median()$data$intensity, median_ram$data$intensity)
  median_json <- jsonlite::fromJSON(output$investigate_median_plot, simplifyVector = FALSE)
  same(as.numeric(unlist(median_json$x$data[[1]]$y)), median_ram$data$intensity)
  reads <- length(signals$batches)
  session$setInputs(investigate_mz_group = 3, investigate_draw_median = 2)
  stopifnot(length(signals$batches) == reads)
  session$setInputs(investigate_time_max = 2)
  stopifnot(grepl("Ranges have changed", output$investigate_median_status))
  session$setInputs(investigate_draw_median = 3)
  stopifnot(investigate_median()$scans == 3L)
  session$setInputs(investigate_time_max = 10)
  session$setInputs(`plotly_click-metabomeister_investigate` = '{"x":100.25,"y":1,"z":10,"curveNumber":0}')
  same(as.numeric(signals$input_updates$investigate_center$value), 100.25)
  session$setInputs(investigate_slice_axis = "time", investigate_time_mode = "combine",
    investigate_time_center = 1, investigate_time_window = 60,
    investigate_spectrum_bin_width = 0.01, investigate_spectrum_count = 5,
    investigate_spectrum_mz_min = 99.9, investigate_spectrum_mz_max = 101,
    investigate_draw_slice = 3)
  stopifnot(investigate_trace()$mode == "combine", investigate_trace()$scan_count == 3L)
  session$setInputs(investigate_time_mode = "overlay", investigate_time_center = 1.2,
    investigate_time_window = 120, investigate_spectrum_count = 3, investigate_draw_slice = 4)
  same(vapply(investigate_trace()$spectra, function(s) s$time, numeric(1)), c(1, 2, 0))
  reads <- length(signals$batches)
  session$setInputs(investigate_y_max = 400, investigate_draw_slice = 5)
  stopifnot(length(signals$batches) == reads)
  session$setInputs(`plotly_click-metabomeister_investigate` = '{"x":100,"y":2,"z":20,"curveNumber":0}')
  same(as.numeric(signals$input_updates$investigate_time_center$value), 2)
  session$setInputs(investigate_mz_max = 90, investigate_draw_3d = 4)
  stopifnot(is.null(investigate_grid()), !session$isClosed(), grepl("minimum < maximum", tail(signals$errors, 1)[[1]]$message))
  session$setInputs(investigate_mz_max = 101, investigate_draw_3d = 5)
  stopifnot(!is.null(investigate_grid()), !session$isClosed())
  closed <- signals$closes
  session$setInputs(investigate_load = 2)
  stopifnot(signals$closes == closed + 1L, is.null(investigate_grid()), is.null(investigate_trace()), is.null(investigate_median()))
  session$setInputs(investigate_clear = 1)
  stopifnot(is.null(investigate_raw()), signals$closes == closed + 2L)
  session$setInputs(investigate_load = 3)
  stopifnot(!is.null(investigate_raw()))
})
stopifnot(signals$opens == signals$closes)
same(unname(tools::md5sum(fixture_path)), hash_before)
cat("PASS: exact m/z traces; native sampling; adjacent grouped means; exact selected-time medians including sparse-zero and disk-spill paths; light palette; combined time-window means and original spectra overlays in nearest-first colour order; bounded batches/grid; actual guarded GUI actions, masked progress functions, cached intensity updates, axis-aware clicks, recovery and handle cleanup.\n")

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
for (pair in list(c("investigate_mz_group", "1"), c("investigate_time_group", "1"))) {
  control <- Filter(function(tag) identical(tag$attribs$id, pair[1]), walk_tags(actual_ui))[[1]]
  stopifnot(as.character(control$attribs$value) == pair[2], control$name == "input")
}
cat("PASS: actual Investigate UI, sidebar order and unique control/output IDs.\n")
