# Rscript --vanilla tests/shared-engine-mz5.R /path/to/real-fixture.mz5
# Exercise the production snapshot and analysis-copy boundary using real MSnbase data.
args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L,file.exists(args[1]))
raw_path <- normalizePath(args[1])
suppressPackageStartupMessages(library(MSnbase))
is_assignment <- function(x,name)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(name))
definition <- Filter(function(x)is_assignment(x,'server'),as.list(parse('server.R')))[[1]]
statements <- as.list(definition[[3]][[3]])[-1]
tasks <- new.env(parent=globalenv())
for(expression in statements)if(is.call(expression)&&identical(expression[[1]],as.name('<-'))&&
  is.call(expression[[3]])&&identical(expression[[3]][[1]],as.name('function')))eval(expression,tasks)
run <- readMSData(raw_path,msLevel=1,mode='inMemory',verbose=FALSE)
rlang::env_binding_unlock(run@assayData) # read_mz5_data() unlocks Engine's assayData.
original <- tasks$investigate_from_engine(run,basename(raw_path))
before <- serialize(run,NULL)
raw_before <- serialize(original$handle,NULL)
tasks$incProgress <- function(...)NULL
tasks$withProgress <- function(message,value,expr)force(expr)
tasks$read_mz5_data <- function(...)stop('Unexpected second file read')
tasks$calibrate_masses <- function(run_data,...) {
  # Deliberately alter the mutable assayData exactly where calibration owns it.
  name <- ls(run_data@assayData)[1]
  run_data@assayData[[name]]@mz <- run_data@assayData[[name]]@mz + 100
  run_data
}
tasks$prepare_extraction_spectra <- function(run_data,...) {
  name <- ls(run_data@assayData)[1]
  stopifnot(all.equal(run_data@assayData[[name]]@mz,run@assayData[[name]]@mz+100)==TRUE)
  stop('CALIBRATION_COPY_CHECKED')
}
folder <- tempfile('engine-copy-');dir.create(folder)
config <- list(mass_df=data.frame(),is_df=data.frame(),parameters_df=data.frame())
error <- tryCatch(tasks$process_mz5_file(raw_path,basename(raw_path),folder,config,
  timing=list(stage=function(...)NULL),shared_data=run),error=function(e)conditionMessage(e))
if (!identical(error,'CALIBRATION_COPY_CHECKED')) stop(error)
stopifnot(identical(error,'CALIBRATION_COPY_CHECKED'),identical(serialize(run,NULL),before),
  identical(serialize(original$handle,NULL),raw_before))
stopifnot(length(original$times)==length(spectra(run)))
ordered <- spectra(run)[order(rtime(run))]
for(i in seq_along(original$times)) {
  stopifnot(identical(unname(original$handle[[i]][,1]),unname(ordered[[i]]@mz)),
    identical(unname(original$handle[[i]][,2]),unname(ordered[[i]]@intensity)))
}
unlink(folder,recursive=TRUE)
cat('PASS: real Engine MSnbase snapshot; production analysis reuses loaded data without reading; calibration mutates only its own copy; original uncalibrated Investigate spectra remain identical.\n')
