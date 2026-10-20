# Run from the application folder: Rscript --vanilla tests/gui-startup.R
# Exercises the actual UI, README renderer and plot interaction handlers with
# small in-memory plot data; no mz5 processing or result files are required.
library(shiny)
library(shinydashboard)
library(DT)
library(shinyjs)
library(waiter)
library(plotly)

is_assignment <- function(x, name) {
  is.call(x) && identical(x[[1]], as.name("<-")) &&
    identical(paste(deparse(x[[2]]), collapse = ""), name)
}
ui_code <- as.list(parse("ui.R"))
server_code <- as.list(parse("server.R"))
server_definition <- Filter(function(x) is_assignment(x, "server"), server_code)[[1]]
statements <- as.list(server_definition[[3]][[3]])[-1]
assignment <- function(name) Filter(function(x) is_assignment(x, name), statements)[[1]]
handler <- function(event) {
  Filter(function(x) is.call(x) && identical(x[[1]], as.name("observeEvent")) &&
    identical(paste(deparse(x[[2]]), collapse = ""), event), statements)
}
warnings <- character()
record_warning <- function(w) {
  warnings <<- c(warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
}
withCallingHandlers({
  title <- "PeakMeister v2.0"
  ui <- eval(Filter(function(x) is_assignment(x, "ui"), ui_code)[[1]])
  stopifnot(inherits(ui, "shiny.tag"))
  # Shiny must send the empty value when these dynamic selectors are cleared.
  # A data-nonempty flag suppresses that event and breaks reselecting the same run.
  for (id in c("existing_run", "results_folder", "results_folder2")) {
    config <- Filter(function(tag) identical(tag$attribs[["data-for"]], id),
      htmltools::tagQuery(ui)$find("script")$selectedTags())
    stopifnot(length(config) == 1L, is.null(config[[1]]$attribs[["data-nonempty"]]))
    select <- htmltools::tagQuery(ui)$find(paste0("#", id))$selectedTags()
    stopifnot(length(select) == 1L, grepl('value=""', as.character(select[[1]]), fixed = TRUE))
  }
}, warning = record_warning)
stopifnot(!any(grepl("Navigation containers", warnings, fixed = TRUE)))

# Run the exact About setup, including its session cleanup callback.
first <- which(vapply(statements, is_assignment, logical(1), name = "readme_dir"))
last <- which(vapply(statements, is_assignment, logical(1), name = "output$readmeContent"))
about <- as.call(c(list(as.name("{")), statements[first:last]))
resource <- NULL
directory <- NULL
withCallingHandlers(testServer(function(input, output, session) {
  eval(about, envir = environment())
}, {
  rendered <- output$readmeContent
  stopifnot(grepl("<iframe", rendered$html, fixed = TRUE))
  resource <<- readme_resource
  directory <<- readme_dir
  stopifnot(identical(shiny::resourcePaths()[[resource]], normalizePath(directory)))
  html <- paste(readLines(file.path(directory, "README.html"), warn = FALSE), collapse = "\n")
  stopifnot(grepl("Software for Pre-Processing MSI-CE-MS Datasets", html, fixed = TRUE),
    grepl("data:(text|application)/javascript", html),
    !grepl("README_files/", html, fixed = TRUE),
    grepl("MathJax", html, fixed = TRUE), grepl("id=\"refs\"", html, fixed = TRUE))
}), warning = record_warning)
stopifnot(!dir.exists(directory), !resource %in% names(shiny::resourcePaths()))

fixture <- list(
  eie_data = data.frame(mt.seconds = c(0, 30, 60, 90), intensity = c(1, 5, 3, 1)),
  annotation_data = data.frame(peak.apex.seconds = 30, peak.height.counts = 5,
    peak.number = 1, comment = ""),
  integration_data = data.frame(mt.seconds = c(0, 30, 60), baseline = 1,
    intensity = c(1, 5, 3)),
  label_data = c("Target", "100"), x_axis_data = c(0, 90), y_axis_data = c(6, 0))
event_block <- as.call(c(list(as.name("{")),
  lapply(c("plot_function_plotly", "plot_event_data", "zoom_state",
    "baseline_adjustment_mode", "plotbaseline_adjustment_mode",
    "output$selected_plot", "output$adjust_baseline_plot"), assignment),
  handler('plot_event_data("plotly_relayout", "peakmeister_selected")')[1],
  handler('plot_event_data("plotly_click", "peakmeister_selected")'),
  handler('plot_event_data("plotly_click", "peakmeister_baseline")')))
withCallingHandlers(testServer(function(input, output, session) {
  plot_render_trigger <- reactiveVal(0)
  db_con <- function() TRUE
  main_folder <- function() tempdir()
  filtered_plot_names <- function() "Fixture_Target"
  read_plot_data_from_db <- function(...) fixture
  adjusted <- reactiveValues(individual = 0L, all = 0L, baseline = NULL)
  adjust_individual_baseline <- function(...) adjusted$individual <- adjusted$individual + 1L
  adjust_plot_baseline <- function(value) {
    adjusted$all <- adjusted$all + 1L
    adjusted$baseline <- value
  }
  eval(event_block, envir = environment())
}, {
  session$flushReact()
  stopifnot(length(session$userData$plotlyShinyEventIDs) == 0L)
  session$setInputs(file_selector = "Fixture.mz5", plot_table_rows_selected = 1L,
    peak_info_table_rows_selected = 1L)
  main <- jsonlite::fromJSON(output$selected_plot)
  modal <- jsonlite::fromJSON(output$adjust_baseline_plot)
  stopifnot(main$x$source == "peakmeister_selected", modal$x$source == "peakmeister_baseline")
  session$flushReact()
  stopifnot("plotly_click-peakmeister_selected" %in% session$userData$plotlyShinyEventIDs,
    "plotly_relayout-peakmeister_selected" %in% session$userData$plotlyShinyEventIDs,
    "plotly_click-peakmeister_baseline" %in% session$userData$plotlyShinyEventIDs)
  session$setInputs(`plotly_relayout-peakmeister_selected` =
    '{"xaxis.range[0]":0.2,"xaxis.range[1]":1,"yaxis.range[0]":0,"yaxis.range[1]":6}')
  stopifnot(zoom_state$xmin == 0.2, zoom_state$xmax == 1)
  session$setInputs(`plotly_relayout-peakmeister_selected` = '{"xaxis.autorange":true}')
  stopifnot(is.null(zoom_state$xmin), is.null(zoom_state$xmax))
  baseline_adjustment_mode(TRUE)
  plotbaseline_adjustment_mode(TRUE)
  session$setInputs(`plotly_click-peakmeister_selected` = '{"x":0.5,"y":5}')
  stopifnot(adjusted$individual == 1L, adjusted$all == 0L,
    session$userData$clicked_position$y == 5)
  # Modal clicks must not trigger the individual-baseline handler on the main plot.
  baseline_adjustment_mode(TRUE)
  session$setInputs(`plotly_click-peakmeister_baseline` = '{"x":0.5,"y":2}')
  stopifnot(adjusted$individual == 1L, adjusted$all == 1L, adjusted$baseline == 2)
}), warning = record_warning)
stopifnot(length(warnings) == 0L)
cat("PASS: full UI without navigation warnings; clearable run/analysis dropdowns; standalone README with embedded assets, equations and references; session resource cleanup; no pre-plot event warnings; zoom/reset and isolated main/dialog baseline clicks.\n")
