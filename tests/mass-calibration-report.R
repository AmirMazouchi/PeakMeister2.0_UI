# Run from the app folder: Rscript --vanilla tests/mass-calibration-report.R
suppressPackageStartupMessages(library(dplyr))
definition <- Filter(function(expr) is.call(expr) && identical(expr[[1]], as.name("<-")) &&
  identical(expr[[2]], as.name("server")), as.list(parse("server.R")))[[1]]
tasks <- new.env()
for (expr in as.list(body(eval(definition[[3]])))[-1]) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) {
    eval(expr, tasks)
  }
}

setClass("CalibrationSpectrum", slots = c(mz = "numeric", intensity = "numeric"))
setClass("CalibrationRun", slots = c(assayData = "environment"))
spectra <- new.env(parent = emptyenv())
peak_shape <- c(1, 2, 10, 2, 1, 1, 2, 10, 2, 1)
spectra$s1 <- new("CalibrationSpectrum", mz = (1:10) + 0.02, intensity = peak_shape)
spectra$s2 <- new("CalibrationSpectrum", mz = (1:10) * 1.005 + 0.01, intensity = peak_shape)
spectra$s3 <- new("CalibrationSpectrum", mz = (1:10) + 0.02,
                  intensity = replace(peak_shape, 8, 2))
run <- new("CalibrationRun", assayData = spectra)
original <- serialize(run, NULL)
tasks$rtime <- function(run_data) seq_along(ls(run_data@assayData))
tasks$incProgress <- function(...) NULL
tasks$withProgress <- function(message, value, expr) force(expr)
tasks$read_mz5_data <- function(...) stop("The shared Engine data should be reused")

params <- data.frame(apply.mass.calibration = "Yes", ref.mass.window.ppm = 10000,
  ref.mass.minimum.counts = 5, ref.mass.one = 3, ref.mass.two = 8)
mass <- data.frame(name = c("Between", "Below"), mz = c(5, 1))
references <- data.frame(name = "Above", mz = 10)
config <- list(mass_df = mass, is_df = references, parameters_df = params)

same <- function(actual, expected) stopifnot(isTRUE(all.equal(actual, expected, tolerance = 1e-10)))
verify_report <- function() {
  work <- tempfile("mass-calibration-")
  dir.create(work)
  original_dir <- setwd(work)
  on.exit({setwd(original_dir); unlink(work, recursive = TRUE)}, add = TRUE)

  tasks$prepare_extraction_spectra <- function(run_data, ...) {
    if (config$parameters_df$apply.mass.calibration == "Yes") {
      # Translation on the first scan, an affine correction with constant ends on the second.
      same(run_data@assayData$s1@mz, as.numeric(1:10))
      same(run_data@assayData$s2@mz, c(0.99, 1.995, 3:8, 9.005, 10.01))
      same(run_data@assayData$s3@mz, spectra$s3@mz)
    } else {
      same(run_data@assayData$s1@mz, spectra$s1@mz)
    }
    stop("extraction reached")
  }
  raw <- "Sample.mz5"
  writeLines("fixture", raw)
  parent <- file.path("Results", "Sample")
  run_analysis <- function() {
    folder <- tasks$create_results_folder(raw)
    error <- tryCatch(tasks$process_mz5_file(raw, raw, folder, config,
      timing = list(stage = function(...) NULL), shared_data = run), error = conditionMessage)
    stopifnot(identical(error, "extraction reached"), identical(serialize(run, NULL), original))
    path <- file.path(parent, paste0("Mass Calibration - ", basename(folder), ".xlsx"))
    stopifnot(file.exists(path), dir.exists(folder), !length(list.files(folder, pattern = "[.]xlsx$")))
    path
  }
  first <- run_analysis()
  checksum <- tools::md5sum(first)
  stopifnot(identical(readxl::excel_sheets(first),
    c("Summary", "Mass List", "Reference Mass List", "Calibration References")))
  read <- function(path, sheet) readxl::read_excel(path, sheet = sheet)
  summary <- read(first, "Summary")
  stopifnot(identical(summary$Value[3:6], c("Applied", "3", "2", "1")))
  targets <- read(first, "Mass List")
  stopifnot(identical(targets$Metabolite, mass$name))
  same(targets[[2]], mass$mz)
  corrected <- c(5 - 0.02, (5 - 0.01) / 1.005)
  same(as.numeric(targets[1, 3:5]), c(min(corrected), mean(corrected), max(corrected)))
  same(as.numeric(targets[2, 3:5]), c(1 - 0.025, 1 - mean(c(0.02, 0.025)), 1 - 0.02))
  refs <- read(first, "Reference Mass List")
  stopifnot(identical(refs$Metabolite, references$name))
  same(as.numeric(refs[1, 2:5]), c(10, 10 - 0.05, 10 - mean(c(0.02, 0.05)), 10 - 0.02))
  locks <- read(first, "Calibration References")
  stopifnot(identical(names(locks), c("Scan (calibration order)", "Scan status",
    "Reference 1 expected m/z", "Reference 1 measured before m/z",
    "Reference 2 expected m/z", "Reference 2 measured before m/z")))
  same(locks[[1]], c(1, 2, 3))
  same(locks[[3]], rep(3, 3))
  same(locks[[4]], c(3.02, 3.025, 3.02))
  same(locks[[5]], rep(8, 3))
  same(locks[[6]], c(8.02, 8.05, NA_real_))
  stopifnot(grepl("^Skipped", locks[[2]][3]))

  config$parameters_df$apply.mass.calibration <- "No"
  second <- run_analysis()
  stopifnot(first != second, identical(tools::md5sum(first), checksum))
  stopifnot(read(second, "Summary")$Value[3] == "Disabled",
    all(is.na(read(second, "Mass List")[3:5])))
  disabled <- read(second, "Calibration References")
  stopifnot(disabled[[2]] == "Disabled", disabled[[3]] == 3, disabled[[5]] == 8,
    all(is.na(disabled[c(4, 6)])))
}
verify_report()
cat("PASS: variable scan corrections, skipped and disabled calibration, workbook contents and rerun-safe placement.\n")
