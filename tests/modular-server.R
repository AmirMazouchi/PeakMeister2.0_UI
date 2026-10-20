# Run from the app folder: Rscript --vanilla tests/modular-server.R
# Uses synthetic traces and real report/SQLite/Plotly code. Raw mz5 decoding,
# the MSnbase API and PNG writing are supplied by small test fixtures.
suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(plotly)
  library(DBI)
  library(RSQLite)
  library(writexl)
  library(jsonlite)
})

is_assignment <- function(x, name) {
  is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(name))
}
definition <- Filter(function(x) is_assignment(x, "server"), as.list(parse("server.R")))[[1]]
statements <- as.list(definition[[3]][[3]])[-1]
tasks <- new.env(parent = globalenv())
for (expr in statements) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
    eval(expr, envir = tasks)
  }
}
tasks$incProgress <- function(...) invisible(NULL)
tasks$showNotification <- function(...) invisible(NULL)
# Exact trapezoid integral for these integer-time fixtures, independent of the
# app's peak selection/filtering algorithms. Production continues to use AUC.
tasks$AUC <- function(x, y, from, to, ...) {
  ix <- which(x >= from & x <= to)
  sum(diff(x[ix]) * (head(y[ix], -1) + tail(y[ix], -1)) / 2)
}
same <- function(actual, expected) stopifnot(isTRUE(all.equal(actual, expected, check.attributes = TRUE)))

parameters <- data.frame(number.of.injections = 3, Manual.Indexes = "No",
  apply.smoothing = "No", apply.mass.calibration = "No", plot.format = "Sample")
references <- data.frame(name = c("100_Low", "200_High"), mz = c(100, 200),
  class = "Reference", extraction.window.ppm = 10, smoothing.kernel = "normal",
  smoothing.strength = 2, min.mt.min = 0, max.mt.min = 6,
  peak.fwhm.tolerance.multiplier = 2, peak.space.tolerance.percent = 20,
  injection1 = c(100, 120) / 60, injection2 = c(200, 220) / 60,
  injection3 = c(300, 320) / 60)
masses <- data.frame(name = c("150_Target", "180_Other"), mz = c(150, 180),
  extraction.window.ppm = 10, smoothing.kernel = "normal", smoothing.strength = 2,
  minimim.peak.width.seconds = 2, migration.window.seconds = 4,
  peak.space.tolerance.percent = 20, snr.threshold = 3,
  interference = NA_character_, interference.comigration.threshold.seconds = 2,
  injection1 = c(110, 114) / 60, injection2 = c(210, 214) / 60,
  injection3 = c(310, 314) / 60)
config <- tasks$prepare_analysis_parameters(masses, references, parameters, data.frame())
same(as.numeric(config$mi_df[1, 5:7]), rep(0.5, 3))
manual <- config$mi_df[1, ]; manual[, 5:7] <- 0.6
manual_parameters <- parameters; manual_parameters$Manual.Indexes <- "Yes"
same(as.numeric(tasks$prepare_analysis_parameters(masses, references, manual_parameters, manual)$mi_df[1, 5:7]), rep(0.6, 3))
duplicate <- masses; duplicate$name[1] <- references$name[1]
stopifnot(inherits(try(tasks$prepare_analysis_parameters(duplicate, references, parameters, data.frame()), silent = TRUE), "try-error"))

times <- 0:360
triangle <- function(apices) rowSums(vapply(apices, function(apex) pmax(0, 100 - abs(times - apex) * 10), numeric(length(times)))) + 1
traces <- lapply(list(c(100, 200, 300), c(120, 220, 320), c(110, 210, 310), c(114, 214, 314)), triangle)
eie <- data.frame(mt.seconds = times, setNames(traces, paste(config$name_vec, "intensity")), check.names = FALSE)
same(tasks$smooth_electropherograms(eie, references, masses, parameters, config$name_vec, 7), eie)
smooth_parameters <- parameters; smooth_parameters$apply.smoothing <- "Yes"
smoothed <- tasks$smooth_electropherograms(eie, references, masses, smooth_parameters, config$name_vec, 7)
stopifnot(identical(smoothed$mt.seconds, eie$mt.seconds), !identical(smoothed[, -1], eie[, -1]))
bad_parameters <- parameters; bad_parameters$apply.smoothing <- "invalid"
stopifnot(inherits(try(tasks$smooth_electropherograms(eie, references, masses, bad_parameters, config$name_vec, 7), silent = TRUE), "try-error"))

standard_peaks <- tasks$process_internal_standard_peaks(eie, references, parameters, config$name_vec, 2, 3, 7)
same(unname(standard_peaks[, c(2, 9)]), data.frame(c(100L, 200L, 300L), c(120L, 220L, 320L), check.names = FALSE) |> unname())
metabolite_peaks <- tasks$process_metabolite_peaks(eie, masses, references, standard_peaks[, c(2, 9)], config$mi_df, config$name_vec, 2, 3, 7)
same(as.numeric(metabolite_peaks[, 2]), c(110, 210, 310))
same(as.numeric(metabolite_peaks[, 7]), rep(1000, 3))
comments <- tasks$flag_peaks_below_lod(eie, metabolite_peaks, masses, 2, 2, 3)
stopifnot(all(as.matrix(comments) == ""))
low_peaks <- metabolite_peaks; low_peaks[1, 7] <- 0; low_peaks[2, 7] <- 0.5
low_comments <- tasks$flag_peaks_below_lod(eie, low_peaks, masses, 2, 2, 3)
stopifnot(low_comments[1, 1] == "NPD", low_comments[2, 1] == "<LOD")
interfered <- masses; interfered$interference[1] <- "180_Other"; interfered$interference.comigration.threshold.seconds[1] <- 5
flagged <- tasks$flag_interfered_peaks(comments, standard_peaks, metabolite_peaks, interfered, references, config$name_vec, 2, 3)
stopifnot(all(flagged$comment_df$`150_Target` == "Interfered"), all(flagged$comment_df$`100_Low` == ""))
peaks <- tasks$flag_interfered_peaks(comments, standard_peaks, metabolite_peaks, masses, references, config$name_vec, 2, 3)
plots <- tasks$build_plot_data(eie, peaks$peaks_df, peaks$comment_df, config$name_vec, c(100, 200, 150, 180), 3)
same(names(plots), config$name_vec)
same(plots$`150_Target`$annotation_data$peak.apex.seconds, c(110L, 210L, 310L))

# Minimal spectrum/trace API for the extraction boundary, without MSnbase.
setClass("Spectrum1", slots = c(mz = "numeric", intensity = "numeric", rt = "numeric", tic = "numeric", centroided = "logical", fromFile = "integer"))
setClass("RunFixture", slots = c(assayData = "environment"))
setClass("TraceFixture", slots = c(rtime = "numeric", intensity = "numeric"))
tasks$mz <- function(x) x@mz
tasks$intensity <- function(x) x@intensity
tasks$rtime <- function(x) if (is(x, "RunFixture")) vapply(as.list(x@assayData), function(s) s@rt, numeric(1)) else x@rt
tasks$assayData <- function(x) x@assayData
tasks$centroided <- function(x) x@centroided
tasks$fromFile <- function(x) x@fromFile
`[.trace_list` <- function(x, i) unclass(x)[[i]]
spectra <- new.env(parent = emptyenv())
spectra$s1 <- new("Spectrum1", mz = c(200, -1, 100), intensity = c(5, 9, 10), rt = 1, tic = 24, centroided = FALSE, fromFile = 1L)
raw_run <- new("RunFixture", assayData = spectra)
prepared <- tasks$prepare_extraction_spectra(raw_run, references, masses, 7)
same(prepared$run_data@assayData$s1@mz, c(100, 200))
same(prepared$mzr[, 1], c(100, 200, 150, 180) * (1 - 10 / 1e6))
tasks$chromatogram <- function(run_data, mz, rt, aggregationFun, missing, msLevel) {
  stopifnot(aggregationFun == "mean", missing == 0, msLevel == 1, nrow(mz) == 4)
  structure(lapply(traces, function(x) new("TraceFixture", rtime = as.numeric(times), intensity = x)), class = "trace_list")
}
extracted <- tasks$extract_electropherograms(prepared$run_data, prepared$mzr)
expected_eie <- eie; expected_eie$mt.seconds <- as.numeric(expected_eie$mt.seconds)
same(tasks$build_eie_dataframe(extracted, config$name_vec), expected_eie)

# Exercise every calibration branch with known, shifted lock masses.
same(tasks$calibrate_masses(raw_run, parameters, 7), raw_run)
cal_spectra <- new.env(parent = emptyenv())
for (i in 1:2) cal_spectra[[paste0("s", i)]] <- new("Spectrum1", mz = (1:10) + 0.02,
  intensity = c(1, 2, 10, 2, 1, 1, 2, 10, 2, 1), rt = as.numeric(i), tic = 32, centroided = FALSE, fromFile = 1L)
cal_params <- parameters; cal_params$apply.mass.calibration <- "Yes"
cal_params$ref.mass.window.ppm <- 10000; cal_params$ref.mass.minimum.counts <- 5
cal_params$ref.mass.one <- 3; cal_params$ref.mass.two <- 8
calibrated <- tasks$calibrate_masses(new("RunFixture", assayData = cal_spectra), cal_params, 7)
same(calibrated@assayData$s1@mz, as.numeric(1:10))

work <- tempfile("modular-server-"); dir.create(work)
old_dir <- setwd(work)
png_calls <- list()
tasks$ggarrange <- function(...) list(...)
tasks$ggsave <- function(filename, width, height, plot, path) {
  png_calls[[length(png_calls) + 1L]] <<- file.path(path, filename)
  if (inherits(plot, "ggplot")) invisible(ggplot2::ggplot_build(plot))
}
tasks$read_mz5_data <- function(...) raw_run
tasks$calibrate_masses <- function(run_data, ...) run_data
tasks$plotly_data <- function(value) invisible(value)
writeLines("synthetic", "raw-one.mz5"); writeLines("synthetic", "raw-two.mz5")
files <- data.frame(FileName = c("Run01.mz5", "Run02.mz5"), FilePath = normalizePath(c("raw-one.mz5", "raw-two.mz5")))
testServer(function(input, output, session) {}, {
  tasks$run_new_analysis(files, masses, references, parameters, data.frame())
})
folders <- lapply(c("Run01", "Run02"), function(run) list.dirs(file.path("Results", run), recursive = FALSE))
stopifnot(all(lengths(folders) == 1), length(png_calls) == 8L)
for (i in seq_along(folders)) {
  folder <- folders[[i]]
  tables <- read.csv(file.path(folder, "Metabolite Peak Areas.csv"), check.names = FALSE)
  same(tables$file.name, c(paste0("Run0", i), "", ""))
  same(tables$`150_Target`, rep(1000L, 3))
  saved <- new.env(); load(file.path(folder, "Data", files$FileName[i], "plot_list_data.RData"), envir = saved)
  stopifnot(setequal(ls(saved), c("plot_list", "peaks_df", "mi_df", "peak_mt_df", "peak_area_df", "comment_df")))
  db <- dbConnect(SQLite(), tasks$get_db_path(folder))
  same(dbGetQuery(db, "SELECT DISTINCT file_name FROM annotation_data")$file_name, paste0("Run0", i))
  same(as.numeric(dbGetQuery(db, "SELECT peak_area FROM peak_area_df WHERE metabolite = '150_Target' ORDER BY peak_number")$peak_area), rep(1000, 3))
  dbDisconnect(db)
  timings <- read.csv(file.path(folder, "run_timings.csv"))
  stopifnot(all(timings$status == "completed"), tail(timings$stage, 1) == "TOTAL RUN (includes all stages)",
    sum(timings$stage == "Build eie_df") == 1, sum(timings$stage == "Write SQLite results") == 1)
  plotly_saved <- new.env(); load(file.path(folder, "plotly_objects.RData"), envir = plotly_saved)
  same(names(plotly_saved$plotly_objects), paste0("Run0", i, "_", config$name_vec))
}
stopifnot(!any(grepl("temp.mz5$", list.files(work, recursive = TRUE))))
metabolite_format <- parameters; metabolite_format$plot.format <- "Metabolite"
tasks$save_static_plots(plots, eie, config$mi_df, metabolite_format, config$name_vec, 2, work, "Run01", "raw-one.mz5")
stopifnot(length(png_calls) == 12L, all(basename(tail(unlist(png_calls), 4)) == "Run01.png"))

# Error reports must record the failing stage and close the temporary raw copy.
tasks$build_eie_dataframe <- function(...) stop("fixture extraction failure")
testServer(function(input, output, session) {}, {
  err <- try(tasks$run_new_analysis(files[1, ], masses, references, parameters, data.frame()), silent = TRUE)
  stopifnot(inherits(err, "try-error"))
})
failed_folder <- list.dirs(file.path("Results", "Run01"), recursive = FALSE)
failed_folder <- setdiff(failed_folder, folders[[1]])
failed <- read.csv(file.path(failed_folder, "run_timings.csv"))
stopifnot(tail(failed$status, 1) == "failed", any(failed$stage == "Build eie_df" & failed$status == "failed"),
  !any(grepl("temp.mz5$", list.files(work, recursive = TRUE))))

# Downstream tables: missing values, sample identity and reporting orientation.
area <- data.frame(file.name = c("Run01", "", "", "Run02", "", ""), peak.number = rep(1:3, 2),
  Reference = c(2, 4, 8, 2, 4, 8), Target = c("10", "<LOD", "NPD", "20", "30", "40"), check.names = FALSE)
replaced <- tasks$replace_missing_peak_areas(area, c("Reference", "Target"), "Minimum Values/5")
same(replaced$Target, c(10, 2, 2, 20, 30, 40))
zeroed <- tasks$replace_missing_peak_areas(area, c("Reference", "Target"), "Missing values = 0")
same(zeroed$Target, c(10, 0, 0, 20, 30, 40))
normalized <- tasks$normalize_peak_areas(replaced, "Reference", c("Reference", "Target"))
same(normalized$Reference, rep(1, 6))
metadata <- data.frame(`Data file` = c("Run01", "Run02"), `1` = c("QC-NIST", "103"),
  `2` = c("101", "104"), `3` = c("102", "105"), check.names = FALSE)
merged <- tasks$connect_peak_area_metadata(replaced, metadata)
same(merged$`PBM Sample ID`, c("QC-NIST", "101", "102", "103", "104", "105"))
matrix <- tasks$build_metaboanalyst_matrix(merged, "QC='Control',Numeric='Sample'")
same(rownames(matrix)[1:2], c("PBM Sample ID", "CLASS"))
same(as.character(unlist(matrix["CLASS", ])), c("Control", rep("Sample", 5)))
same(unname(tasks$calculate_cv(data.frame(x = c(1, 2, 3)))), 50)
stopifnot(inherits(tasks$build_control_chart(c(1, 2, 3), "Target", "(Raw Data)"), "plotly"))
setwd(old_dir); unlink(work, recursive = TRUE)
cat("PASS: modular scientific stages with synthetic traces; multi-file reports/SQLite/Plotly/timings; PNG paths in both formats; failure timing/temporary-file cleanup; downstream transformations/metadata/reporting.\n")
