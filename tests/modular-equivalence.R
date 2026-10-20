# Compare with a saved copy of server.R from commit 3604a34.
# Rscript --vanilla tests/modular-equivalence.R /path/to/previous/server.R
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, file.exists(args[1]))
baseline_path <- normalizePath(args[1])
source('tests/modular-server.R')

before <- Filter(function(x) is_assignment(x, 'server'), as.list(parse(baseline_path)))[[1]]
before_statements <- as.list(before[[3]][[3]])[-1]
before_init <- Filter(function(x) is.call(x) && identical(x[[1]], as.name('observeEvent')) && identical(x[[2]], quote(input$initialize)), before_statements)[[1]][[3]]
normalize <- function(x) {
  if (is.character(x) && any(grepl('\n', x, fixed = TRUE))) {
    return(vapply(x, function(y) paste(trimws(strsplit(y, '\n', fixed = TRUE)[[1]]), collapse = '\n'), character(1), USE.NAMES = FALSE))
  }
  if (is.call(x)) return(as.call(lapply(as.list(x), normalize)))
  x
}
for (pair in list(c('delete_peak','delete_selected_peaks'), c('deletemiddlepeaks','delete_middle_peaks'),
  c('deleteallpeaks','delete_all_peaks'), c('undo','undo_peak_deletion'),
  c('regenerateplots','regenerate_changed_plots'), c('interferences','mark_peak_interferences'),
  c('load_batchcorr_packages','load_batch_correction_packages'), c('batch_corr','apply_batch_correction'),
  c('drug_databases','load_drug_databases'), c('confirm_parse_drugbank','parse_and_save_drugbank'))) {
  original <- Filter(function(x) is.call(x) && identical(x[[1]], as.name('observeEvent')) && identical(x[[2]], as.call(list(as.name('$'), as.name('input'), as.name(pair[1])))), before_statements)[[1]][[3]]
  stopifnot(identical(normalize(original), normalize(body(tasks[[pair[2]]]))))
}

# Restore the extracted builder after the intentional failing-stage test.
eval(Filter(function(x) is_assignment(x, 'build_eie_dataframe'), statements)[[1]], envir = tasks)
comparison <- tempfile('modular-equivalence-'); dir.create(comparison)
original_dir <- getwd(); setwd(comparison)
writeLines('synthetic raw fixture', 'one.mz5'); writeLines('synthetic raw fixture', 'two.mz5')
input_files <- data.frame(FileName = c('Run01.mz5','Run02.mz5'),
  FilePath = normalizePath(c('one.mz5','two.mz5')), Source = 'upload')
baseline <- new.env(parent = tasks)
baseline$analysis_mode <- function() 'new'
baseline$uploadedmz5 <- function() input_files
baseline$massData <- function() masses
baseline$refMassListData <- function() references
baseline$parametersData <- function() parameters
baseline$userMTI <- function() data.frame()
baseline$readMSData <- function(...) raw_run
baseline$isMSnbaseVerbose <- function() FALSE
baseline$env_binding_unlock <- function(...) invisible(NULL)
baseline_runner <- eval(call('function', as.pairlist(alist()), before_init), envir = baseline)
testServer(function(input, output, session) {}, {
  baseline_runner()
  tasks$run_new_analysis(input_files, masses, references, parameters, data.frame())
})
for (run in c('Run01','Run02')) {
  folders <- list.dirs(file.path('Results',run), recursive = FALSE)
  stopifnot(length(folders) == 2)
  first <- folders[!grepl(' 2$', folders)]; second <- setdiff(folders, first)
  for (filename in c('Migration Index Summary.csv','Metabolite Peak Areas.csv','Metabolite Migration Times.csv')) {
    same(readLines(file.path(first,filename)), readLines(file.path(second,filename)))
  }
  a <- new.env(); b <- new.env()
  load(file.path(first,'Data',paste0(run,'.mz5'),'plot_list_data.RData'), a)
  load(file.path(second,'Data',paste0(run,'.mz5'),'plot_list_data.RData'), b)
  same(ls(a),ls(b))
  for (name in ls(a)) same(a[[name]], b[[name]])
  ca <- dbConnect(SQLite(), tasks$get_db_path(first)); cb <- dbConnect(SQLite(), tasks$get_db_path(second))
  same(dbListTables(ca), dbListTables(cb))
  for (table in dbListTables(ca)) {
    arows <- dbReadTable(ca, table); brows <- dbReadTable(cb, table)
    if (table == 'edit_log') { arows$timestamp <- NULL; brows$timestamp <- NULL }
    same(arows,brows)
  }
  dbDisconnect(ca); dbDisconnect(cb)
  a <- new.env(); b <- new.env()
  load(file.path(first,'plotly_objects.RData'), a); load(file.path(second,'plotly_objects.RData'), b)
  same(names(a$plotly_objects), names(b$plotly_objects))
  for (name in names(a$plotly_objects)) {
    ap <- plotly_build(a$plotly_objects[[name]]); bp <- plotly_build(b$plotly_objects[[name]])
    same(ap$x$data, bp$x$data); same(ap$x$layout, bp$x$layout)
  }
  at <- read.csv(file.path(first,'run_timings.csv')); bt <- read.csv(file.path(second,'run_timings.csv'))
  same(at[,c('file','stage','status')],bt[,c('file','stage','status')])
}
setwd(original_dir); unlink(comparison, recursive = TRUE)
cat('PASS: old/new two-file pipeline produces identical CSVs, saved analysis objects, all SQLite tables, Plotly traces/layouts and stage/status sequences; ten moved GUI action bodies match.\n')
