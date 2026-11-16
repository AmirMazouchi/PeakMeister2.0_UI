# Rscript --vanilla tests/shared-analysis-equivalence.R /path/to/pre-change/server.R
# Compare the existing modular pipeline with shared Engine-data processing.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L,file.exists(args[1]))
baseline_path <- normalizePath(args[1])
source('tests/modular-server.R')
before <- Filter(function(x)is_assignment(x,'server'),as.list(parse(baseline_path)))[[1]]
before_statements <- as.list(before[[3]][[3]])[-1]
for(name in c('calibrate_masses','prepare_extraction_spectra','extract_electropherograms',
  'build_eie_dataframe','smooth_electropherograms','process_internal_standard_peaks',
  'process_metabolite_peaks','flag_peaks_below_lod','flag_interfered_peaks')) {
  old <- Filter(function(x)is_assignment(x,name),before_statements)[[1]][[3]]
  current <- Filter(function(x)is_assignment(x,name),statements)[[1]][[3]]
  stopifnot(identical(old,current))
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
baseline$process_mz5_file <- eval(Filter(function(x)is_assignment(x,'process_mz5_file'),before_statements)[[1]][[3]],envir=baseline)
baseline$isMSnbaseVerbose <- function() FALSE
baseline$env_binding_unlock <- function(...) invisible(NULL)
baseline_old_runner <- eval(Filter(function(x)is_assignment(x,'run_new_analysis'),before_statements)[[1]][[3]],envir=baseline)
baseline_runner <- function()baseline_old_runner(input_files,masses,references,parameters,data.frame())
testServer(function(input, output, session) {}, {
  baseline_runner()
  before_shared <- serialize(raw_run,NULL)
  tasks$run_new_analysis(input_files, masses, references, parameters, data.frame(),
    loaded=setNames(rep(list(list(data=raw_run)),nrow(input_files)),input_files$FilePath))
  stopifnot(identical(serialize(raw_run,NULL),before_shared))
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
cat('PASS: old/new two-file pipeline produces identical CSVs, saved analysis objects, all SQLite tables, Plotly traces/layouts and stage/status sequences; shared Engine objects remain untouched. Editing/history safety is tested separately.\n')
