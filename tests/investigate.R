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
signals$read_scans <- integer(); signals$batches <- integer(); signals$errors <- list(); signals$input_updates <- list()
fixture_opener <- function(path) { signals$opens <- signals$opens + 1L; new.env() }
fixture_closer <- function(handle) { signals$closes <- signals$closes + 1L; invisible(NULL) }
fixture_reader <- function(handle, scans) {
  signals$batches <- c(signals$batches, length(scans))
  signals$read_scans <- c(signals$read_scans, scans)
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
occupied <- matrix(FALSE, nrow = 11L, ncol = length(centers))
for (row in seq_along(raw$scans)) {
  scan <- fixture_scans[[raw$scans[row]]]
  for (column in seq_along(centers)) {
    keep <- is.finite(scan[, 1]) & is.finite(scan[, 2]) & scan[, 1] >= edges[column] &
      if (column == length(centers)) scan[, 1] <= edges[column + 1] else scan[, 1] < edges[column + 1]
    occupied[row, column] <- any(keep)
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
stopifnot(grepl("too large", error_text(tasks$investigate_build_grid(raw, c(99.9, 111), c(0, 10), alignment_width = 0.0001))))
stopifnot(grepl("positive alignment", error_text(tasks$investigate_mass_axis(raw, c(99.9, 101), 0))))

# Independent per-scan interpolation oracle, retaining measured zeros and
# exterior missing cells. The 3D means above continue to use the unfilled data.
median_oracle <- matrix_oracle
for (row in seq_len(nrow(median_oracle))) {
  anchors <- which(occupied[row, ])
  if (length(anchors) < 2L) next
  for (pair in seq_len(length(anchors) - 1L)) {
    left <- anchors[pair]; right <- anchors[pair + 1L]
    if (right <= left + 1L || median_oracle[row, left] == 0 || median_oracle[row, right] == 0) next
    for (column in seq.int(left + 1L, right - 1L)) {
      fraction <- (centers[column] - centers[left])/(centers[right] - centers[left])
      median_oracle[row, column] <- median_oracle[row, left] + fraction *
        (median_oracle[row, right] - median_oracle[row, left])
    }
  }
}
same(raw$alignment_width, 0.01)
exact_width <- tasks$investigate_mass_axis(raw, c(921.9, 922.1), 0.01)
stopifnot(length(exact_width$mz) == 20L, all(diff(exact_width$edges) > 0))
same(tasks$investigate_alignment_value(0.01006), 0.0101)
stopifnot(grepl("at least", error_text(tasks$investigate_alignment_value(0.00001))))
axis_test <- list(mz = 0:8)
fill <- tasks$investigate_interpolate_empty(list(columns = c(9L, 2L, 4L, 7L),
  values = c(80, 10, 20, 50), points = 4L), axis_test)
same(fill$columns, 2:9); same(fill$values, c(10, 15, 20, 30, 40, 50, 65, 80))
stopifnot(fill$points == 4L) # No extrapolation into the first cell.
zero <- list(columns = c(2L, 4L, 6L), values = c(10, 0, 30), points = 3L)
same(tasks$investigate_interpolate_empty(zero, list(mz = 0:6)), zero)
short_cell <- tasks$investigate_interpolate_empty(list(columns = c(1L, 3L), values = c(10, 24), points = 2L),
  list(mz = c(0.5, 1.5, 2.25)))
same(short_cell$values, c(10, 18, 24))
# Signed bounds remain valid and a single anchor cannot create invented edges.
signed <- tasks$investigate_interpolate_empty(list(columns = c(1L, 3L), values = c(-4, 4), points = 2L), list(mz = 0:2))
same(signed$values, c(-4, 0, 4))
single <- list(columns = 2L, values = 10, points = 1L)
same(tasks$investigate_interpolate_empty(single, list(mz = 0:4)), single)

# Exact medians, including signed signals, odd/even scan counts, sparse zeros,
# and the disk-spill path. All paths agree with independent stats::median.
for (n in 1:12) for (trial in 1:30) {
  set.seed(n*100 + trial)
  observed <- sample(c(-9, -1, 0, 0, 0, 2, 8), n, replace = TRUE)
  same(tasks$investigate_zero_median(observed[observed != 0], n), median(observed))
}
median_ram <- tasks$investigate_median_spectrum(raw, c(99.9, 101), c(0, 10), width)
median_disk <- tasks$investigate_median_spectrum(raw, c(99.9, 101), c(0, 10), width, cache_bytes = 64)
same(median_ram$data$intensity, apply(median_oracle, 2, median))
same(median_disk$data, median_ram$data)
stopifnot(median_ram$method == "cached sparse values", median_disk$method == "partitioned sparse values")
median_crop <- tasks$investigate_median_spectrum(raw, c(99.9, 101), c(1, 2), width)
same(median_crop$data$intensity, apply(median_oracle[2:3, , drop = FALSE], 2, median))
stopifnot(median_crop$scans == 2L)

# Subsampling is applied within the median's own time window BEFORE I/O.
for (budget in c(64, 1e6)) {
  read_before <- length(signals$read_scans)
  sampled <- tasks$investigate_median_spectrum(raw, c(99.9, 101), c(1, 9), width,
    cache_bytes = budget, scan_step = 3)
  kept_rows <- c(2L, 5L, 8L)
  same(sampled$data$intensity, apply(median_oracle[kept_rows, , drop = FALSE], 2, median))
  same(sampled$selected_scans, raw$scans[kept_rows])
  same(sampled$selected_times, c(1, 4, 7))
  stopifnot(sampled$scans == 3, sampled$available_scans == 9, sampled$scan_step == 3)
  read_scans <- signals$read_scans[seq.int(read_before + 1L, length(signals$read_scans))]
  stopifnot(all(read_scans %in% sampled$selected_scans),
    all(sampled$selected_scans %in% read_scans))
}
small_sample <- tasks$investigate_median_spectrum(raw, c(99.9, 101), c(0, 1), width, scan_step = 10)
same(small_sample$data$intensity, as.numeric(median_oracle[1, ]))
stopifnot(small_sample$scans == 1, small_sample$available_scans == 2)
for (value in list(0, -1, 1.5, NA_real_, Inf)) {
  stopifnot(grepl("median scan factor", error_text(tasks$investigate_median_spectrum(raw,
    c(99.9, 101), c(0, 10), width, scan_step = value))))
}
# Bounds display at most three decimals and retain the original endpoints.
bounds <- tasks$investigate_display_range(c(0.055433333, 49.99313333))
same(bounds, c(0.055, 49.994))
stopifnot(bounds[1] <= 0.055433333, bounds[2] >= 49.99313333)

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
# A retained empty scan contributes zeros, even if all other scans are skipped.
empty_sample <- tasks$investigate_median_spectrum(native, c(99, 104), c(0.5, 5), scan_step = 10)
stopifnot(empty_sample$scans == 1L, empty_sample$available_scans == 3L,
  all(empty_sample$data$intensity == 0))
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

# A large median is displayed as a bounded SVG overview. Keep actual extrema,
# including a one-position positive spike next to a negative spike, unchanged.
large <- list(data = data.frame(mz = seq(50, 1700, length.out = 200000), intensity = numeric(200000)),
  mz_range = c(50, 1700), intensity_range = c(-100, 1000), key = "large-median")
large$data$intensity[c(12345, 12346, 100500, 199990)] <- c(999, -99, 500, 750)
overview <- tasks$investigate_median_view_data(large)
stopifnot(nrow(overview$data) <= 4000L, overview$visible == 200000L,
  overview$shown < overview$visible, max(overview$data$intensity) == 999,
  min(overview$data$intensity) == -99)
rows <- match(overview$data$mz, large$data$mz)
same(overview$data$intensity, large$data$intensity[rows])
svg <- plotly_build(tasks$investigate_plot_median(overview))
stopifnot(svg$x$data[[1]]$type == "scatter", svg$x$layout$dragmode == "zoom",
  svg$x$source == "metabomeister_median", isFALSE(svg$x$data[[1]]$line$simplify),
  "plotly_relayout" %in% svg$x$shinyEvents, !"plotly_hover" %in% svg$x$shinyEvents)
zoom <- list(x = large$data$mz[c(10001, 11000)], y = c(-10, 10))
detail <- tasks$investigate_median_view_data(large, zoom)
stopifnot(detail$visible == 1000, detail$shown == 1000)
same(detail$data$mz, large$data$mz[10000:11001])
same(detail$data$intensity, large$data$intensity[10000:11001])
next_view <- tasks$investigate_median_zoom(list(`xaxis.range[0]` = zoom$x[1],
  `xaxis.range[1]` = zoom$x[2], `yaxis.range` = zoom$y), large)
same(next_view, zoom)
same(tasks$investigate_median_zoom(list(autosize = TRUE), large, zoom), zoom)
same(tasks$investigate_median_zoom(list(`xaxis.range` = c(NA, Inf)), large, zoom), zoom)
same(tasks$investigate_median_zoom(list(`xaxis.autorange` = TRUE, `yaxis.autorange` = TRUE), large, zoom),
  list(x = large$mz_range, y = large$intensity_range))
blank <- tasks$investigate_median_view_data(large, list(x = c(1800, 1900), y = c(0, 1)))
stopifnot(blank$visible == 0, blank$shown == 0)

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
    investigate_median_time_min = 0, investigate_median_time_max = 10, investigate_median_step = 1,
    investigate_slice_time_min = 0, investigate_slice_time_max = 10,
    investigate_mz_group = 2, investigate_time_group = 3, investigate_alignment_width = 0.1, investigate_median_alignment_width = 0.1,
    investigate_slice_axis = "mz",
    investigate_center = 100, investigate_window = 10, investigate_unit = "ppm",
    investigate_auto_z = TRUE, investigate_auto_y = TRUE, investigate_log_z = FALSE,
    investigate_aggregation = "mean")
  session$setInputs(investigate_draw_3d = 1)
  stopifnot(is.null(investigate_grid()), tail(signals$errors, 1)[[1]]$action == "Investigate: draw 3D preview", !session$isClosed())
  session$setInputs(investigate_load = 1)
  stopifnot(investigate_raw()$name == "Fixture.mz5", is.null(investigate_grid()), is.null(investigate_trace()), is.null(investigate_median()))
  same(as.numeric(signals$input_updates$investigate_mz_group$value), 40)
  same(as.numeric(signals$input_updates$investigate_time_group$value), 20)
  same(as.numeric(signals$input_updates$investigate_alignment_width$value), 0.01)
  same(as.numeric(signals$input_updates$investigate_median_alignment_width$value), 0.01)
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
    grepl("Minimum native Δm/z", output$investigate_status),
    grepl("minimum Δt", output$investigate_status))
  stopifnot(!grepl("including empty scans|positive|Measured across", output$investigate_status),
    grepl("Da; time ≈ 180 s", output$investigate_group_widths))
  session$setInputs(investigate_draw_median = 1)
  same(investigate_median()$data$intensity, median_ram$data$intensity)
  median_json <- jsonlite::fromJSON(output$investigate_median_plot, simplifyVector = FALSE)
  same(as.numeric(unlist(median_json$x$data[[1]]$y)), median_ram$data$intensity)
  # Real relayout events zoom the display from the cached full median only.
  reads <- length(signals$batches)
  session$setInputs(`plotly_relayout-metabomeister_median` =
    '{"xaxis.range[0]":99.99,"xaxis.range[1]":100.2,"yaxis.range[0]":0,"yaxis.range[1]":20}')
  same(investigate_median_view(), list(x = c(99.99, 100.2), y = c(0, 20)))
  stopifnot(investigate_median_display()$shown == 2L, length(signals$batches) == reads)
  same(investigate_median()$data$intensity, median_ram$data$intensity)
  plotted <- jsonlite::fromJSON(output$investigate_median_plot, simplifyVector = FALSE)
  stopifnot(plotted$x$data[[1]]$type == "scatter", plotted$x$source == "metabomeister_median")
  # View updates use trace-only proxy messages, never rebuild the widget.
  same(plotted, median_json)
  session$setInputs(investigate_reset_median = 1)
  stopifnot(is.null(investigate_median_view()), investigate_median_display()$shown == 11L,
    length(signals$batches) == reads)
  reads <- length(signals$batches)
  session$setInputs(investigate_mz_group = 3, investigate_alignment_width = 0.2, investigate_draw_median = 2)
  stopifnot(length(signals$batches) == reads)
  same(investigate_median()$alignment_width, 0.1)
  old_grid <- investigate_grid()
  session$setInputs(investigate_median_alignment_width = 0.05, investigate_draw_median = 20)
  same(investigate_median()$alignment_width, 0.05)
  same(investigate_grid(), old_grid)
  session$setInputs(investigate_median_alignment_width = 0.1)
  session$setInputs(investigate_time_max = 2)
  same(as.numeric(signals$input_updates$investigate_median_time_max$value), 2)
  # Simulate receipt of the default update, then customize the median window.
  session$setInputs(investigate_median_time_max = 2, investigate_draw_median = 3)
  stopifnot(investigate_median()$scans == 3L)
  session$setInputs(investigate_median_time_edited = 1, investigate_median_time_min = 1, investigate_median_time_max = 9,
    investigate_median_step = 3)
  stopifnot(grepl("Ranges or scan factor have changed", output$investigate_median_status))
  session$setInputs(investigate_draw_median = 4)
  same(investigate_median()$selected_times, c(1, 4, 7))
  default_update <- signals$input_updates$investigate_median_time_max$value
  reads <- length(signals$batches)
  session$setInputs(investigate_time_max = 10, investigate_draw_median = 5)
  # Changing 3D time cannot overwrite an independent median window or its cache.
  same(signals$input_updates$investigate_median_time_max$value, default_update)
  stopifnot(length(signals$batches) == reads, investigate_median()$scans == 3L)
  session$setInputs(investigate_median_step = 0, investigate_draw_median = 6)
  stopifnot(!session$isClosed(), grepl("median scan factor", tail(signals$errors, 1)[[1]]$message))
  session$setInputs(investigate_median_step = 1, investigate_draw_median = 7)
  stopifnot(investigate_median()$scans == 9L, !session$isClosed())
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
cat("PASS: independent fixed 0.01/four-decimal 3D/median widths; per-scan linear filling against an independent oracle, recorded zeros/edges preserved; exact m/z traces; native sampling; adjacent grouped means; independent-window/subsampled medians with skipped-scan I/O checks, sparse-zero and disk-spill paths; three-decimal bounds and default grouping/physical widths; bounded peak-preserving SVG median with full-detail cached zoom/reset; light palette; combined time-window means and original spectra overlays in nearest-first colour order; bounded batches/grid; actual guarded GUI actions, masked progress functions, cached intensity updates, axis-aware clicks, recovery and handle cleanup.\n")

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
for (pair in list(c("investigate_mz_group", "40"), c("investigate_time_group", "20"), c("investigate_median_step", "10"), c("investigate_alignment_width", "0.01"), c("investigate_median_alignment_width", "0.01"))) {
  control <- Filter(function(tag) identical(tag$attribs$id, pair[1]), walk_tags(actual_ui))[[1]]
  stopifnot(as.character(control$attribs$value) == pair[2], control$name == "input")
}
cat("PASS: actual Investigate UI, sidebar order and unique control/output IDs.\n")

ui_text <- as.character(htmltools::renderTags(actual_ui)$html)
stopifnot(!grepl("Time is shown in minutes|Saved analyses do not contain|investigate_group_status", ui_text),
  !grepl("Time uses the typical scan interval|Last groups may be smaller", ui_text))
cat("PASS: requested redundant UI messages removed.\n")
