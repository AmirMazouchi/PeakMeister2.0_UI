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
  grid <- tasks$investigate_build_grid(raw, limits, raw$time_range, 100L, 80L)
  same(max(grid$z), max(c(0, unlist(values))))
  stopifnot(grid$scans == length(raw$scans), all(dim(grid$z) == c(80L, 100L)))
})
same(unname(tools::md5sum(fixture)), hash)
cat("PASS: genuine mz5 file opened by mzR; metadata/time units/MS1 selection; exact mean/sum/max slices independently checked; bounded grid maxima; closed handle; original file unchanged.\n")
