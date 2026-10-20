# Rscript --vanilla tests/metadata-location.R
# Real load/connect observers, shared metadata precedence and legacy fallback.
suppressPackageStartupMessages({ library(shiny); library(dplyr); library(tidyr) })
is_assignment <- function(x, name) is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(name))
server <- Filter(function(x) is_assignment(x, "server"), as.list(parse("server.R")))[[1]]
statements <- as.list(server[[3]][[3]])[-1]
assignment <- function(name) Filter(function(x) is_assignment(x, name), statements)[[1]]
event <- function(name) Filter(function(x) is.call(x) && identical(x[[1]], as.name("observe_action")) &&
  identical(x[[2]], as.call(list(as.name("$"), as.name("input"), as.name(name)))), statements)[[1]]
work <- tempfile("metadata-location-"); dir.create(work)
a1 <- file.path(work, "Results", "Run01", "Results 2026-10-20")
a2 <- file.path(work, "Results", "Run01", "Results 2026-10-20 3")
b1 <- file.path(work, "Results", "Run02", "Results 2026-10-20")
legacy <- file.path(work, "Results 2026-10-19")
make_analysis <- function(folder, run) {
  dir.create(file.path(folder, "Data"), recursive = TRUE)
  area <- data.frame(file.name = c(run, ""), peak.number = 1:2, `100_A` = c(10, 20), `200_B` = c(30, 40), check.names = FALSE)
  write.csv(area, file.path(folder, "Metabolite Peak Areas.csv"), row.names = FALSE)
  write.csv(area, file.path(folder, "Metabolite Migration Times.csv"), row.names = FALSE)
}
for (folder in c(a1, a2)) make_analysis(folder, "Run01")
make_analysis(b1, "Run02"); make_analysis(legacy, "OldRun")
write_metadata <- function(path, run, prefix) {
  write.csv(data.frame(`Data file` = run, `1` = paste0(prefix, "1"), `2` = paste0(prefix, "2"), check.names = FALSE), path, row.names = FALSE)
}
shared_a <- file.path(dirname(a1), "metadata.csv")
shared_b <- file.path(dirname(b1), "metadata.csv")
write_metadata(shared_a, "Run01", "SharedA")
write_metadata(shared_b, "Run02", "SharedB")
write_metadata(file.path(a1, "Data", "metadata.csv"), "Run01", "LocalA")
write_metadata(file.path(legacy, "Data", "metadata.csv"), "OldRun", "Legacy")
write_metadata(file.path(work, "metadata.csv"), "Wrong", "MustNotUse")
notices <- new.env(); notices$messages <- character()
testServer(function(input, output, session) {
  names <- c("action_label", "report_action_error", "with_action_error", "observe_action", "observe_updates",
    "metadata_file_paths", "connect_peak_area_metadata", "grab_metabolite_files",
    "peak_areas_data", "migration_times_data", "main_folder2")
  for (name in names) eval(assignment(name), envir = environment())
  showNotification <- function(ui, ...) { notices$messages <- c(notices$messages, as.character(ui)); invisible(NULL) }
  for (name in c("results_folder2", "load_migration_area", "connect_metadata")) eval(event(name), envir = environment())
}, {
  stopifnot(identical(metadata_file_paths(a1), metadata_file_paths(paste0(a1, "/"))))
  stopifnot(identical(metadata_file_paths(a1)[1], shared_a), identical(metadata_file_paths(a2)[1], shared_a))
  stopifnot(identical(metadata_file_paths(legacy), file.path(legacy, "Data", "metadata.csv")))
  session$setInputs(results_folder2 = a1, file_selector2 = "Unrelated.mz5")
  session$setInputs(load_migration_area = 1, connect_metadata = 1)
  stopifnot(identical(peak_areas_data()$`PBM Sample ID`, c("SharedA1", "SharedA2")),
    file.exists(file.path(a1, "merged_peak_area.csv")), !file.exists(file.path(dirname(a1), "merged_peak_area.csv")))
  session$setInputs(results_folder2 = a2)
  session$setInputs(load_migration_area = 2, connect_metadata = 2)
  stopifnot(identical(peak_areas_data()$`PBM Sample ID`, c("SharedA1", "SharedA2")))
  write_metadata(shared_a, "Run01", "EditedA")
  session$setInputs(load_migration_area = 3, connect_metadata = 3)
  stopifnot(identical(peak_areas_data()$`PBM Sample ID`, c("EditedA1", "EditedA2")))
  session$setInputs(results_folder2 = b1)
  session$setInputs(load_migration_area = 4, connect_metadata = 4)
  stopifnot(identical(peak_areas_data()$`PBM Sample ID`, c("SharedB1", "SharedB2")))
  unlink(shared_a)
  session$setInputs(results_folder2 = a1)
  session$setInputs(load_migration_area = 5, connect_metadata = 5)
  stopifnot(identical(peak_areas_data()$`PBM Sample ID`, c("LocalA1", "LocalA2")), !file.exists(shared_a))
  session$setInputs(results_folder2 = a2)
  session$setInputs(load_migration_area = 6)
  previous <- peak_areas_data()
  previous_export <- unname(tools::md5sum(file.path(a2, "merged_peak_area.csv")))
  session$setInputs(connect_metadata = 6)
  stopifnot(identical(peak_areas_data(), previous), grepl(shared_a, tail(notices$messages, 1), fixed = TRUE),
    identical(unname(tools::md5sum(file.path(a2, "merged_peak_area.csv"))), previous_export), !session$isClosed())
  session$setInputs(results_folder2 = legacy)
  session$setInputs(load_migration_area = 7, connect_metadata = 7)
  stopifnot(identical(peak_areas_data()$`PBM Sample ID`, c("Legacy1", "Legacy2")), !session$isClosed())
})
unlink(work, recursive = TRUE)
cat("PASS: shared metadata reused across analysis versions; shared precedence; selected-run isolation; edited shared file reloaded; local/legacy fallback; missing-file guidance preserves data/exports; merged results stay analysis-specific.\n")
