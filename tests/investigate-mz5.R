# Rscript --vanilla tests/investigate-mz5.R /path/to/real-fixture.mz5
# Integration check with a genuine mz5/HDF5 file and the installed mzR reader.
# A small reference fixture is available in ProteoWizard's example_data:
# https://github.com/ProteoWizard/pwiz/blob/master/example_data/tiny.pwiz.1.1.mz5
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, file.exists(args[1]))
if (!requireNamespace("mzR", quietly = TRUE)) stop("This integration check requires mzR.")
source("tests/investigate.R")
fixture <- normalizePath(args[1])
hash <- unname(tools::md5sum(fixture))
local({
  raw <- tasks$investigate_open_raw(fixture, basename(fixture))
  on.exit(tasks$investigate_close_raw(raw), add = TRUE)
  header <- mzR::header(raw$handle)
  stopifnot(length(raw$scans) > 0L, all(header$msLevel[raw$scans] == 1L))
  same(raw$times, header$retentionTime[raw$scans] / 60)
  scans <- lapply(raw$scans, function(index) mzR::peaks(raw$handle, scans = index))
  cleaned <- lapply(scans, tasks$investigate_clean_scan)
  nonempty <- which(vapply(cleaned, nrow, integer(1)) > 0L)
  stopifnot(length(nonempty) > 0L)
  first <- cleaned[[nonempty[1]]]
  center <- first[which(first[, 1] > 0)[1], 1]
  width <- 0.01
  limits <- tasks$investigate_window_limits(center, width, "Da")
  values <- lapply(scans, function(scan) {
    scan <- as.matrix(scan)
    good <- is.finite(scan[, 1]) & is.finite(scan[, 2]) &
      scan[, 1] >= limits[1] & scan[, 1] <= limits[2]
    scan[good, 2]
  })
  for (method in c("mean", "sum", "max")) {
    expected <- vapply(values, function(x) if (!length(x)) 0 else switch(method,
      mean = mean(x), sum = sum(x), max = max(x)), numeric(1))
    actual <- tasks$investigate_extract_slice(raw, center, width, "Da", raw$time_range, method)
    same(actual$data$intensity, expected)
  }
  stopifnot(raw$sampling$shared, length(raw$sampling$mz_axis) == 15L,
    raw$sampling$mz_delta == 1, raw$sampling$time_count == 3L)
  full <- matrix(0, nrow = length(scans), ncol = 15L)
  for (row in seq_along(cleaned)) if (nrow(cleaned[[row]])) {
    full[row, match(cleaned[[row]][, 1], raw$sampling$mz_axis)] <- cleaned[[row]][, 2]
  }
  grid <- tasks$investigate_build_grid(raw, raw$mz_range, raw$time_range, 4L, 2L)
  same(dim(grid$z), c(2L, 4L))
  for (row in 1:2) for (column in 1:4) {
    rows <- seq.int((row - 1)*2 + 1, min(row*2, 3))
    columns <- seq.int((column - 1)*4 + 1, min(column*4, 15))
    same(grid$z[row, column], mean(full[rows, columns]))
  }
  for (budget in c(64, 1e6)) {
    med <- tasks$investigate_median_spectrum(raw, raw$mz_range, raw$time_range, cache_bytes = budget)
    same(med$data$intensity, apply(full, 2, median))
  }
  for (step in c(2, 10)) {
    med <- tasks$investigate_median_spectrum(raw, raw$mz_range, raw$time_range, scan_step = step)
    retained <- seq.int(1L, nrow(full), by = step)
    same(med$data$intensity, apply(full[retained, , drop = FALSE], 2, median))
    same(med$selected_scans, raw$scans[retained])
  }
  time_center <- raw$times[nonempty[1]]
  time_window <- max(abs(raw$times - time_center)) * 60 + 1
  combined <- tasks$investigate_extract_time_slice(raw, time_center, raw$mz_range,
    mode = "combine", window_seconds = time_window, mz_bin_width = 1)
  expected <- numeric(nrow(combined$data))
  for (scan in cleaned) {
    if (!nrow(scan)) next
    bins <- pmin(length(expected), floor((scan[, 1] - raw$mz_range[1])/1) + 1L)
    per_scan <- tapply(scan[, 2], bins, mean)
    indices <- as.integer(names(per_scan))
    expected[indices] <- expected[indices] + as.numeric(per_scan)
  }
  same(combined$data$intensity, expected/length(raw$scans))
  overlays <- tasks$investigate_extract_time_slice(raw, time_center, raw$mz_range,
    mode = "overlay", window_seconds = time_window, n_spectra = 5)
  for (spectrum in overlays$spectra) {
    original <- cleaned[[match(spectrum$scan, raw$scans)]]
    original <- original[order(original[, 1]), , drop = FALSE]
    same(spectrum$data$mz, original[, 1])
    same(spectrum$data$intensity, original[, 2])
  }
  distances <- vapply(overlays$spectra, function(spectrum) spectrum$distance_seconds, numeric(1))
  stopifnot(!is.unsorted(distances))
})
same(unname(tools::md5sum(fixture)), hash)
cat("PASS: genuine mz5 MS1 metadata/reading; exact m/z traces; independently checked combined-spectrum averages and original nearest-first overlays; native sampling, adjacent grouped means, exact sparse medians, closed handle and unchanged raw file.\n")
