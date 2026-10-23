# Rscript --vanilla tests/saved-analysis-loading.R
# Real SQLite/RData fixtures and the actual saved-analysis load observer.
suppressPackageStartupMessages({
  library(shiny); library(DBI); library(RSQLite); library(dplyr); library(tidyr)
})
is_assignment <- function(x, name) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name(name))
server_code <- Filter(function(x) is_assignment(x, "server"), as.list(parse("server.R")))[[1]]
statements <- as.list(server_code[[3]][[3]])[-1]
assignment <- function(name) Filter(function(x) is_assignment(x, name), statements)[[1]]
tasks <- new.env(parent = globalenv())
for (expr in statements) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.call(expr[[3]]) &&
      identical(expr[[3]][[1]], as.name("function"))) eval(expr, tasks)
}
same <- function(x, y) stopifnot(isTRUE(all.equal(x, y, check.attributes = TRUE)))
error_text <- function(expr) tryCatch({ force(expr); stop("Expected a failure") },
                                     error = conditionMessage)
work <- tempfile("saved-analysis-loading-")
dir.create(work)
make_fixture <- function(run, version, offset = 0) {
  folder <- file.path(work, "Results", run, version)
  path <- file.path(folder, "Data", paste0(run, ".mz5"))
  dir.create(path, recursive = TRUE)
  plot_list <- setNames(lapply(c("100_Reference", "150_Target"), function(metab) list(
    eie_data = data.frame(mt.seconds = 1:5, intensity = (1:5) * 10 + offset),
    label_data = c(metab, sub("_.*", "", metab))
  )), c("100_Reference", "150_Target"))
  save(plot_list, file = file.path(path, "plot_list_data.RData"))
  con <- dbConnect(SQLite(), tasks$get_db_path(folder))
  tasks$initialize_db(con)
  for (metab in names(plot_list)) {
    dbAppendTable(con, "plot_axis_data", data.frame(file_name = run, metabolite = metab,
      x_start = 1, x_end = 5, y_max_annot = 50, y_max_eie = 50,
      label_name = metab, label_mz = sub("_.*", "", metab)))
    dbAppendTable(con, "annotation_data", data.frame(file_name = run, metabolite = metab,
      peak_number = 1:2, comment = c("saved edit", ""), peak_apex_secs = c(2, 4),
      peak_height = c(20, 40)))
    dbAppendTable(con, "peaks_df", data.frame(file_name = run, metabolite = metab,
      peak_number = 1:2, start_seconds = c(1, 3), apex_seconds = c(2, 4),
      end_seconds = c(3, 5), start_intensity = c(10, 30), apex_intensity = c(20, 40),
      end_intensity = c(30, 50), peak_area = c("25", "35")))
    dbAppendTable(con, "peak_area_df", data.frame(file_name = run, metabolite = metab,
      peak_number = 1:2, peak_area = c("25", "35")))
    dbAppendTable(con, "peak_mt_df", data.frame(file_name = run, metabolite = metab,
      peak_number = 1:2, migration_time = c(2, 4)))
  }
  dbDisconnect(con)
  folder
}
first <- make_fixture("Run01", "Results 2026-10-19")
second <- make_fixture("Run01", "Results 2026-10-20", 100)
other <- make_fixture("Run02", "Results 2026-10-20", 200)
con <- tasks$open_db(first)
path <- file.path(first, "Data", "Run01.mz5", "plot_list_data.RData")
original <- tasks$load_saved_plot_list(first, "Run01")
loaded <- tasks$load_plot_list_from_db(con, first, "Run01")
same(loaded$`150_Target`$eie_data[[2]], (1:5) * 10)
same(loaded$`150_Target`$annotation_data$comment, c("saved edit", ""))

# Reproduce the exact reported exception before exercising the repaired lookup.
renamed <- original[2:1]; names(renamed) <- c("old target key", "old reference key")
old_raw <- renamed[["150_Target"]]$eie_data
old_frame <- data.frame(mt.seconds = old_raw[[1]], intensity = old_raw[[2]])
stopifnot(grepl("'names' attribute \\[2\\].*vector \\[0\\]",
                error_text(colnames(old_frame)[2] <- "150_Target intensity")))
plot_list <- renamed
save(plot_list, file = path)
loaded <- tasks$load_plot_list_from_db(con, first, "Run01")
same(loaded$`150_Target`$eie_data[[2]], (1:5) * 10)
same(loaded$`150_Target`$annotation_data$comment, c("saved edit", ""))
# No database or RData repair writes: edits and the original snapshot are retained.
same(tasks$load_saved_plot_list(first, "Run01"), renamed)
same(dbReadTable(con, "annotation_data")$comment, rep(c("saved edit", ""), 2))

matrix_plots <- original
matrix_plots[[1]]$eie_data <- as.matrix(matrix_plots[[1]]$eie_data)
same(tasks$read_saved_electropherogram(matrix_plots, "100_Reference", path)$intensity, (1:5) * 10)
positional <- original; positional[[1]] <- unname(positional[[1]])
same(tasks$read_saved_electropherogram(positional, "100_Reference", path)$intensity, (1:5) * 10)
for (bad in list(NULL, data.frame(), data.frame(time = 1:3),
                data.frame(time = character(), intensity = character()),
                data.frame(time = 1:3, intensity = letters[1:3]))) {
  broken <- original; broken[["150_Target"]]["eie_data"] <- list(bad)
  msg <- error_text(tasks$read_saved_electropherogram(broken, "150_Target", path))
  stopifnot(grepl("150_Target", msg), grepl(path, msg, fixed = TRUE),
            !grepl("'names' attribute", msg, fixed = TRUE))
}
missing <- original[1]
stopifnot(grepl("no matching saved trace", error_text(
  tasks$read_saved_electropherogram(missing, "150_Target", path))))
ambiguous <- renamed; ambiguous[[3]] <- ambiguous[[1]]
stopifnot(grepl("multiple saved traces", error_text(
  tasks$read_saved_electropherogram(ambiguous, "150_Target", path))))
dbDisconnect(con)

load_event <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("observe_action")) &&
  identical(paste(deparse(x[[2]]), collapse = ""), "list(db_con(), input$file_selector)"), statements)
stopifnot(length(load_event) == 1L)
guard_names <- c("action_label", "report_action_error", "with_action_error", "observe_action", "observe_updates")
state_names <- c("plot_list_data_values", "run_metadata", "previous_metadata", "modified_peak_plots",
                 "plot_render_trigger", "editing_context", "edit_history_trigger", "main_folder", "db_con")
counts <- new.env(); counts$loads <- 0L; counts$errors <- character()
harness <- function(input, output, session) {
  results_folders_reactive <- reactive(c(first, second, other))
  for (name in c(guard_names, state_names)) eval(assignment(name), environment())
  for (name in c("load_saved_plot_list", "read_saved_electropherogram", "read_plot_data_from_db",
                 "load_wide_data_from_db", "pending_peak_plots", "clear_loaded_analysis")) eval(assignment(name), environment())
  original_loader <- tasks$load_plot_list_from_db
  load_plot_list_from_db <- function(...) {
    counts$loads <- counts$loads + 1L
    original_loader(...)
  }
  report_action_error <- function(action, error) {
    counts$errors <- c(counts$errors, conditionMessage(error))
  }
  eval(load_event[[1]], environment())
}
testServer(harness, {
  connection <- tasks$open_db(first)
  db_con(connection)
  session$setInputs(results_folder = first, file_selector = "Run01.mz5")
  stopifnot(counts$loads == 1L, !length(counts$errors), length(plot_list_data_values$plot_list) == 2L)
  same(plot_list_data_values$comment_df$`150_Target`, c("saved edit", ""))
  # Same raw file, different dated analysis: the connection alone triggers loading.
  dbDisconnect(connection); connection <- tasks$open_db(second)
  db_con(connection); session$setInputs(results_folder = second)
  stopifnot(counts$loads == 2L)
  same(plot_list_data_values$plot_list$`150_Target`$eie_data[[2]], (1:5) * 10 + 100)
  # Switching run and connection together executes one load.
  dbDisconnect(connection); connection <- tasks$open_db(other)
  db_con(connection); session$setInputs(results_folder = other, file_selector = "Run02.mz5")
  stopifnot(counts$loads == 3L, run_metadata$file_name == "Run02.mz5")
  same(plot_list_data_values$plot_list$`150_Target`$eie_data[[2]], (1:5) * 10 + 200)
  # An unrecoverable mismatch shows one useful error and clears stale edit state.
  plot_list <- missing; save(plot_list, file = path)
  previous_metadata$comment_df <- data.frame(old = 1)
  dbDisconnect(connection); connection <- tasks$open_db(first)
  db_con(connection); session$setInputs(results_folder = first, file_selector = "Run01.mz5")
  stopifnot(counts$loads == 4L, length(counts$errors) == 1L,
    grepl("150_Target", counts$errors), is.null(plot_list_data_values$plot_list),
    is.null(run_metadata$file_name), is.null(previous_metadata$comment_df))
  # A valid selection still opens in the same session after that failure.
  dbDisconnect(connection); connection <- tasks$open_db(second)
  db_con(connection); session$setInputs(results_folder = second)
  stopifnot(counts$loads == 5L, run_metadata$file_name == "Run01.mz5")
  dbDisconnect(connection)
})
# The subsequently uploaded file contained exactly these six NULL objects (with
# mi_df saved twice by the old Undo handler). Reject it before per-trace lookup.
empty_snapshot <- list2env(setNames(rep(list(NULL), 6),
  c("plot_list", "comment_df", "mi_df", "peaks_df", "peak_mt_df", "peak_area_df")))
save(list = c("plot_list", "comment_df", "mi_df", "peaks_df", "peak_mt_df", "peak_area_df", "mi_df"),
     envir = empty_snapshot, file = path)
stopifnot(grepl("No saved plot_list", error_text(tasks$load_saved_plot_list(first, "Run01"))))
unlink(work, recursive = TRUE)
cat("PASS: saved SQLite/RData loads; exact reported failure; unique-label recovery; malformed/missing/ambiguous traces; single load/error; edits preserved; analysis/run switching and recovery.\n")
