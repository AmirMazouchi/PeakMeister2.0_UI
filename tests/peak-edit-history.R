# Rscript --vanilla tests/peak-edit-history.R
# Real SQLite/RData/CSV files and actual GUI handlers; no raw mz5 decoding.
suppressPackageStartupMessages({ library(shiny); library(DBI); library(RSQLite); library(dplyr); library(tidyr); library(plotly); library(jsonlite) })
is_assignment <- function(x, name) is.call(x) && identical(x[[1]], as.name("<-")) && identical(x[[2]], as.name(name))
definition <- Filter(function(x) is_assignment(x, "server"), as.list(parse("server.R")))[[1]]
statements <- as.list(definition[[3]][[3]])[-1]
assignment <- function(name) Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) &&
  identical(paste(deparse(x[[2]]), collapse = ""), name), statements)[[1]]
event <- function(name) Filter(function(x) is.call(x) && identical(x[[1]], as.name("observe_action")) && identical(x[[2]], as.call(list(as.name("$"), as.name("input"), as.name(name)))), statements)[[1]]
functions <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("<-")) && is.call(x[[3]]) && identical(x[[3]][[1]], as.name("function")), statements)
helpers <- new.env(parent = globalenv()); for (expr in functions) eval(expr, helpers)
same <- function(x, y) stopifnot(isTRUE(all.equal(x, y, check.attributes = TRUE)))
work <- tempfile("peak-history-"); dir.create(work); old <- setwd(work)
make_fixture <- function(run) {
  folder <- file.path(".", "Results", run, "Results 2026-10-19")
  data <- file.path(folder, "Data", paste0(run, ".mz5")); dir.create(data, recursive = TRUE)
  names <- c("100_Target", "200_Other"); apices <- c(5, 15, 25)
  trace <- data.frame(mt.seconds = 0:30, intensity = 1 + vapply(0:30, function(t) sum(pmax(0, 30 - abs(t - apices) * 10)), numeric(1)))
  plot_list <- setNames(lapply(names, function(name) {
    eie <- trace; colnames(eie)[2] <- paste(name, "intensity")
    list(eie_data = eie, annotation_data = data.frame(peak.number = 1:3, comment = "", peak.apex.seconds = apices, peak.height.counts = 31),
      integration_data = do.call(rbind, lapply(1:3, function(i) data.frame(peak.number = i, mt.seconds = (apices[i]-2):(apices[i]+2), intensity = c(11,21,31,21,11), baseline = 1))),
      label_data = c(name, sub("_.*", "", name)), x_axis_data = c(0,30), y_axis_data = c(31,31))
  }), names)
  peaks_df <- as.data.frame(do.call(cbind, lapply(names, function(name) {
    out <- data.frame(start.seconds = apices-2, apex.seconds = apices, end.seconds = apices+2,
      start_intensity = 11, apex_intensity = 31, end_intensity = 11, peak.area = 80)
    colnames(out) <- paste0(name, ".", colnames(out)); out
  })))
  peak_area_df <- data.frame(file.name = c(run,"",""), peak.number = 1:3, check.names = FALSE)
  peak_mt_df <- peak_area_df
  for (name in names) { peak_area_df[[name]] <- 80; peak_mt_df[[name]] <- apices/60 }
  comment_df <- setNames(data.frame(a = rep("",3), b = rep("",3)), names)
  mi_df <- data.frame(name = names, left_is = "Reference", right_is = "none", description = "Reference", peak1 = 1, peak2 = 1, peak3 = 1)
  helpers$save_analysis_snapshot(data, plot_list, peaks_df, mi_df, peak_mt_df, peak_area_df, comment_df)
  helpers$save_analysis_database(folder, run, plot_list, peaks_df, peak_area_df, peak_mt_df, comment_df, mi_df)
  helpers$export_result_csvs(folder, peak_area_df, peak_mt_df)
  writeLines("original PNG", file.path(folder, "original.png"))
  folder
}
original <- make_fixture("Run01"); other <- make_fixture("Run02")
signals <- new.env(); signals$messages <- character(); signals$errors <- character(); signals$updates <- list()
showNotification <- function(ui, ...) { signals$messages <- c(signals$messages, as.character(ui)); invisible(NULL) }
updateSelectInput <- function(session, inputId, ...) { signals$updates[[inputId]] <- list(...); invisible(NULL) }
state_names <- c("plot_list_data_values", "run_metadata", "previous_metadata", "modified_peak_plots", "plot_render_trigger", "db_con",
  "editing_context", "edit_history_trigger", "peak_edit_in_progress", "main_folder", "filtered_plot_names", "line_positions", "zoom_state")
load_event <- Filter(function(x) is.call(x) && identical(x[[1]], as.name("observe_action")) && identical(paste(deparse(x[[2]]), collapse=""), "list(db_con(), input$file_selector)"), statements)[[1]]
harness <- function(input, output, session) {
  for (expr in functions) eval(expr, environment())
  selected_saved_run <- function() sub("\\.mz5$", "", input$file_selector)
  selected_results_filename <- function() input$file_selector
  uploadedmz5 <- function() data.frame(FileName = input$file_selector)
  results_folders_reactive <- reactive(list_results_folders(input$file_selector, selected_saved_run()))
  for (name in state_names) eval(assignment(name), environment())
  report_action_error <- function(action, error) { signals$errors <- c(signals$errors, conditionMessage(error)) }
  AUC <- function(x,y,...) sum(diff(x)*(head(y,-1)+tail(y,-1))/2)
  for (name in c("delete_peak","deletemiddlepeaks","deleteallpeaks","manual_integrate","confirm_peak_number","undo","redo","interferences","baseline_peak_number_adjustment")) eval(event(name), environment())
  eval(event("results_folder"), environment())
  eval(load_event, environment())
  eval(assignment("output$edit_history_status"), environment())
  eval(assignment("output$selected_plot"), environment())
  eval(assignment("output$plot_table"), environment())
  eval(assignment("output$peak_info_table"), environment())
  eval(assignment("output$modified_peak_plots_table"), environment())
  session$onSessionEnded(function() { con <- isolate(db_con()); if (!is.null(con) && dbIsValid(con)) dbDisconnect(con) })
}
manifest <- function(folder) {
  files <- list.files(folder, recursive = TRUE, all.files = TRUE, full.names = TRUE)
  files <- files[!grepl("results.sqlite-(wal|shm)$", files)]
  setNames(unname(tools::md5sum(files)), substring(files, nchar(folder)+2))
}
edited <- NULL
testServer(harness, {
  db_con(open_db(original)); session$setInputs(results_folder = original, file_selector = "Run01.mz5", plot_table_rows_selected = 1L)
  baseline <- capture_peak_state(); original_manifest <- manifest(original)
  session$setInputs(undo = 1, redo = 1, delete_peak = 1)
  stopifnot(!length(signals$errors), nrow(dbGetQuery(db_con(),"SELECT * FROM annotation_data WHERE comment = 'NPD'")) == 0,
    length(list.dirs(dirname(original), recursive = FALSE)) == 1L)
  # Copy failure cannot touch the original or leave a partially copied analysis.
  assign("file.copy", function(...) stop("copy fixture failure"), envir = environment(ensure_editable_analysis))
  session$setInputs(peak_info_table_rows_selected = c(1L,2L), delete_peak = 2)
  rm("file.copy", envir = environment(ensure_editable_analysis))
  stopifnot(any(grepl("copy fixture failure", signals$errors)), length(list.dirs(dirname(original), recursive = FALSE)) == 1L)
  signals$errors <- character()
  session$setInputs(delete_peak = 3)
  edited <<- main_folder()
  stopifnot(edited != original, file.exists(file.path(edited,".peakmeister-edit-copy.rds")), peak_history_position(db_con(),"Run01")$undo == 1L,
    identical(signals$updates$results_folder$selected, edited), !length(signals$errors))
  same(manifest(original), original_manifest)
  stopifnot(all(plot_list_data_values$comment_df$`100_Target`[1:2] == "NPD"))
  packed <- unserialize(dbGetQuery(db_con(),"SELECT before_state FROM peak_edit_history")$before_state[[1]])
  stopifnot(is.null(packed$plot_list$`100_Target`$eie_data))
  # The browser catches up to the new folder; later edits stay in the same copy.
  session$setInputs(results_folder = edited)
  stopifnot(identical(modified_peak_plots$names, "Run01_100_Target"))
  line_positions$red <- 6/60; line_positions$blue <- 8/60
  session$setInputs(peak_info_table_rows_selected = 1L, manual_integrate = 1)
  stopifnot(!length(signals$errors), peak_history_position(db_con(),"Run01")$undo == 2L,
    plot_list_data_values$plot_list$`100_Target`$annotation_data$peak.apex.seconds[1] == 6,
    plot_list_data_values$peak_mt_df$`100_Target`[1] == 6/60)
  session$userData$clicked_position <- list(x = 25/60, y = 3)
  session$setInputs(baseline_peak_number = "3", baseline_peak_number_adjustment = 1)
  stopifnot(!length(signals$errors), peak_history_position(db_con(),"Run01")$undo == 3L)
  state_after_three <- capture_peak_state()
  # Repeated Undo traverses mixed edits and restores actual SQL, snapshots and CSVs.
  session$setInputs(undo = 2)
  stopifnot(peak_history_position(db_con(),"Run01")$undo == 2L)
  session$setInputs(undo = 3)
  stopifnot(peak_history_position(db_con(),"Run01")$undo == 1L, plot_list_data_values$comment_df$`100_Target`[1] == "NPD")
  session$setInputs(undo = 4)
  stopifnot(peak_history_position(db_con(),"Run01")$undo == 0L, peak_history_position(db_con(),"Run01")$redo == 3L)
  same(plot_list_data_values$comment_df, baseline$comment_df)
  stopifnot(all(dbGetQuery(db_con(),"SELECT peak_area FROM peak_area_df")$peak_area == "80"))
  saved <- new.env(); load(file.path(edited,"Data","Run01.mz5","plot_list_data.RData"), saved)
  same(saved$comment_df, baseline$comment_df)
  csv <- read.csv(file.path(edited,"Metabolite Peak Areas.csv"), check.names = FALSE)
  stopifnot(all(csv$`100_Target` == 80))
  mt_csv <- read.csv(file.path(edited,"Metabolite Migration Times.csv"), check.names = FALSE)
  same(mt_csv$`100_Target`, c(5,15,25)/60)
  for (click in 2:4) session$setInputs(redo = click)
  same(plot_list_data_values$peaks_df, state_after_three$peaks_df)
  same(plot_list_data_values$plot_list$`100_Target`$integration_data, state_after_three$plot_list$`100_Target`$integration_data)
  # All baselines in a plot count as one step, as do middle/all deletions and flags.
  adjust_plot_baseline(2)
  stopifnot(peak_history_position(db_con(),"Run01")$undo == 4L)
  session$setInputs(undo = 5)
  session$setInputs(deletemiddlepeaks = 1)
  # That middle peak was already deleted: a no-op must retain the Redo branch.
  stopifnot(peak_history_position(db_con(),"Run01")$undo == 3L, peak_history_position(db_con(),"Run01")$redo == 1L)
  session$setInputs(plot_table_rows_selected = 2L, deletemiddlepeaks = 2)
  stopifnot(peak_history_position(db_con(),"Run01")$undo == 4L, peak_history_position(db_con(),"Run01")$redo == 0L)
  session$setInputs(deleteallpeaks = 1)
  stopifnot(peak_history_position(db_con(),"Run01")$undo == 5L)
  session$setInputs(plot_table_rows_selected = 1L, peak_info_table_rows_selected = 2L, interferences = 1)
  stopifnot(peak_history_position(db_con(),"Run01")$undo == 6L)
  # A real SQL write failure rolls back partial table writes and records no action.
  before_failed <- capture_peak_state(); position <- peak_history_position(db_con(),"Run01")
  data_before <- dbReadTable(db_con(), "annotation_data")
  assign("dbAppendTable", function(conn, name, value, ...) {
    if (name == "peak_area_df") stop("database fixture failure")
    DBI::dbAppendTable(conn, name, value, ...)
  }, envir = environment(save_plot_data))
  session$setInputs(peak_info_table_rows_selected = 1L, delete_peak = 4)
  rm("dbAppendTable", envir = environment(save_plot_data))
  same(dbReadTable(db_con(),"annotation_data"), data_before)
  same(peak_history_position(db_con(),"Run01"), position)
  same(plot_list_data_values$comment_df, before_failed$comment_df)
  stopifnot(any(grepl("database fixture failure", signals$errors)))
  signals$errors <- character()
  # Failed Undo snapshot saving must also leave the SQL cursor/data unchanged.
  normal_save <- save_current_peak_snapshot
  failed_once <- FALSE
  assign("save_current_peak_snapshot", function() {
    normal_save()
    if (!failed_once) { failed_once <<- TRUE; stop("snapshot fixture failure after replacement") }
  }, envir = environment(step_peak_history))
  session$setInputs(undo = 6)
  assign("save_current_peak_snapshot", normal_save, envir = environment(step_peak_history))
  same(peak_history_position(db_con(),"Run01"), position)
  same(dbReadTable(db_con(),"annotation_data"), data_before)
  restored_file <- new.env(); load(file.path(edited,"Data","Run01.mz5","plot_list_data.RData"), restored_file)
  same(restored_file$comment_df, before_failed$comment_df)
  stopifnot(any(grepl("snapshot fixture failure", signals$errors)))
  signals$errors <- character()
  same(manifest(original), original_manifest)
  # A successful action after failures still works in this session.
  session$setInputs(undo = 7)
  stopifnot(peak_history_position(db_con(),"Run01")$undo == 5L, !length(signals$errors))
})
# Reopening an edited folder restores its history, including the Redo branch.
testServer(harness, {
  db_con(open_db(edited)); session$setInputs(results_folder = edited, file_selector = "Run01.mz5", plot_table_rows_selected = 1L)
  stopifnot(peak_history_position(db_con(),"Run01")$undo == 5L, peak_history_position(db_con(),"Run01")$redo == 1L)
  stopifnot(length(modified_peak_plots$names) > 0L)
  session$setInputs(redo = 1)
  stopifnot(peak_history_position(db_con(),"Run01")$undo == 6L)
  # A different run has no access to this run's history.
  dbDisconnect(db_con()); db_con(open_db(other)); clear_loaded_analysis()
  session$setInputs(results_folder = other, file_selector = "Run02.mz5", undo = 1)
  stopifnot(peak_history_position(db_con(),"Run02")$undo == 0L, !file.exists(file.path(other,".peakmeister-edit-copy.rds")))
})
setwd(old); unlink(work, recursive = TRUE)
cat("PASS: original analysis preserved; one copy per editing session; grouped/mixed edits; persistent multi-step Undo/Redo; redo branching; database/snapshot/CSV restoration; copy/SQL/file failures and retries; isolated run histories.\n")
