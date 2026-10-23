# Rscript --vanilla tests/undo-safety.R
# Exercise the actual Undo observer and prevent empty state from overwriting data.
suppressPackageStartupMessages({ library(shiny); library(DBI); library(RSQLite) })
is_assignment <- function(x, name) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(x[[2]], as.name(name))
server_code <- Filter(function(x) is_assignment(x, "server"), as.list(parse("server.R")))[[1]]
statements <- as.list(server_code[[3]][[3]])[-1]
assignment <- function(name) Filter(function(x) is_assignment(x, name), statements)[[1]]
undo_event <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("observe_action")) &&
                       identical(x[[2]], quote(input$undo)), statements)[[1]]
same <- function(x, y) stopifnot(isTRUE(all.equal(x, y, check.attributes = TRUE)))
work <- tempfile("undo-safety-"); path <- file.path(work, "Data", "Run01.mz5")
dir.create(path, recursive = TRUE)
snapshot_path <- file.path(path, "plot_list_data.RData")
make_snapshot <- function(value) list(
  plot_list = list(Target = list(eie_data = data.frame(mt.seconds = 1:3, intensity = value + 1:3),
                                label_data = c("Target", "100"))),
  peaks_df = data.frame(Target = 1:3), mi_df = data.frame(name = "Target"),
  peak_mt_df = data.frame(Target = 1:3), peak_area_df = data.frame(Target = 1:3),
  comment_df = data.frame(Target = rep("", 3)))
snapshot <- make_snapshot(10)
save(list = names(snapshot), file = snapshot_path, envir = list2env(snapshot))
before <- readBin(snapshot_path, "raw", n = file.info(snapshot_path)$size)
counts <- new.env(); counts$messages <- character(); counts$errors <- character()
showNotification <- function(ui, ...) { counts$messages <- c(counts$messages, as.character(ui)); invisible(NULL) }
names_to_load <- c("action_label", "report_action_error", "with_action_error", "observe_action",
  "read_saved_electropherogram", "validate_analysis_snapshot", "save_analysis_snapshot",
  "save_plot_data", "undo_peak_deletion", "step_peak_history", "peak_history_position", "clear_loaded_analysis", "initialize_db",
  "editing_context", "edit_history_trigger", "peak_edit_in_progress", "plot_list_data_values", "run_metadata", "previous_metadata", "modified_peak_plots", "plot_render_trigger", "db_con")
harness <- function(input, output, session) {
  for (name in names_to_load) eval(assignment(name), environment())
  main_folder <- reactive(input$results_folder)
  report_action_error <- function(action, error) { counts$errors <- c(counts$errors, conditionMessage(error)) }
  eval(undo_event, environment())
}
testServer(harness, {
  connection <- dbConnect(SQLite(), file.path(work, "results.sqlite"))
  initialize_db(connection); db_con(connection)
  dbAppendTable(connection, "annotation_data", data.frame(file_name = "Run01", metabolite = "Target",
    peak_number = 1, comment = "retain database row", peak_apex_secs = 2, peak_height = 20))
  database_before <- lapply(dbListTables(connection), function(table) dbReadTable(connection, table))
  put_state <- function(state, value) for (key in names(value)) state[[key]] <- value[[key]]
  put_state(run_metadata, snapshot); run_metadata$file_name <- "Run01.mz5"
  put_state(plot_list_data_values, snapshot); plot_list_data_values$file_name <- "Run01"
  # This is the reported loss path: a loaded file and no previous deletion state.
  session$setInputs(results_folder = work, file_selector = "Run01.mz5", undo = 1)
  stopifnot(length(counts$messages) == 1L, !length(counts$errors))
  same(readBin(snapshot_path, "raw", n = file.info(snapshot_path)$size), before)
  same(run_metadata$plot_list, snapshot$plot_list)
  same(plot_list_data_values$plot_list, snapshot$plot_list)
  # A save with no snapshot is refused before opening/truncating the existing file.
  failed <- try(save_analysis_snapshot(path, NULL, NULL, NULL, NULL, NULL, NULL), silent = TRUE)
  stopifnot(inherits(failed, "try-error"))
  same(readBin(snapshot_path, "raw", n = file.info(snapshot_path)$size), before)
  # Nor may an invalid edit erase SQLite rows before snapshot validation fails.
  plot_list_data_values$plot_list <- NULL
  stopifnot(inherits(try(save_plot_data("Target"), silent = TRUE), "try-error"))
  same(lapply(dbListTables(connection), function(table) dbReadTable(connection, table)), database_before)
  same(readBin(snapshot_path, "raw", n = file.info(snapshot_path)$size), before)
  put_state(plot_list_data_values, snapshot)
  # Old in-memory history cannot be mistaken for persistent committed actions.
  put_state(previous_metadata, make_snapshot(20)); previous_metadata$file_name <- "Run01.mz5"
  session$setInputs(undo = 2)
  same(readBin(snapshot_path, "raw", n = file.info(snapshot_path)$size), before)
  same(run_metadata$plot_list, snapshot$plot_list)
  # A valid snapshot remains saveable without NULLs or duplicate names.
  restored <- make_snapshot(30)
  do.call(save_analysis_snapshot, c(list(subfolder_path = path), restored))
  loaded <- new.env(); saved_names <- load(snapshot_path, loaded)
  stopifnot(length(saved_names) == 6L, !anyDuplicated(saved_names))
  same(as.list(loaded)[names(restored)], restored)
  # Clearing/switching analyses invalidates Undo; a later click leaves the file intact.
  good_bytes <- readBin(snapshot_path, "raw", n = file.info(snapshot_path)$size)
  clear_loaded_analysis(); session$setInputs(undo = 6)
  same(readBin(snapshot_path, "raw", n = file.info(snapshot_path)$size), good_bytes)
  stopifnot(length(counts$messages) == 3L, is.null(run_metadata$plot_list))
  dbDisconnect(connection)
})
unlink(work, recursive = TRUE)
cat("PASS: empty Undo and obsolete one-slot state cannot overwrite saved files; invalid direct edits preserve SQLite; snapshot validation and cleared-session safety remain intact.\n")
