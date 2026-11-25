# Rscript --vanilla tests/investigate-cache.R
# Independent extrema oracle, persistence and actual shared-data/live GUI handlers.
suppressPackageStartupMessages({library(shiny); library(plotly); library(DBI); library(RSQLite)})
same <- function(a,b) stopifnot(isTRUE(all.equal(a,b,check.attributes=TRUE,tolerance=1e-10)))
is_assignment <- function(x,name) is.call(x) && identical(x[[1]],as.name('<-')) &&
  identical(paste(deparse(x[[2]]),collapse=''),name)
definition <- Filter(function(x)is_assignment(x,'server'),as.list(parse('server.R')))[[1]]
statements <- as.list(definition[[3]][[3]])[-1]
assignment <- function(name) Filter(function(x)is_assignment(x,name),statements)[[1]]
tasks <- new.env(parent=globalenv())
for (expression in statements) if (is.call(expression) && identical(expression[[1]],as.name('<-')) &&
  is.call(expression[[3]]) && identical(expression[[3]][[1]],as.name('function'))) eval(expression,tasks)
workspace <- tempfile('fine-spectrum-');dir.create(workspace)
times <- c(0,1.2,1.5,2,2,2.9,4,4.4,5,5.3)
scans <- lapply(seq_along(times),function(i) {
  x <- c(162.0002,162.0022,162.0032,162.0052,162.0092)
  y <- c(i,-i,0,2*i,3*i)
  if (i == 7L) { x <- x[-2];y <- y[-2] }
  if (i == 8L) { x <- x[-5];y <- y[-5] }
  cbind(x,y)
})
raw <- tasks$investigate_memory_raw(scans,times,'Fixture.mz5')
cache <- tasks$investigate_cache_build(raw,file.path(workspace,'Results','Fixture'))
stopifnot(!is.null(cache),all(file.exists(cache$paths)),
  identical(basename(cache$paths),c('max_mz_spec_dmz0.001.sqlite','min_mz_spec_dmz0.001.sqlite')))
same(cache$meta$settings$breaks,c((0:31)/6,5.3))
stopifnot(sum(cache$meta$counts)==length(times))
# Oracle bins readings globally, then fills only empty cells between measured
# nonzero neighbours. Zeros measured in a cell block interpolation across it.
settings <- cache$meta$settings
centers <- (settings$first+seq_len(settings$bins)-0.5)/1000
oracle <- t(vapply(scans,function(scan){
  bins <- pmin(settings$bins,pmax(1L,floor(scan[,1]*1000+1e-8)-settings$first+1L))
  means <- tapply(scan[,2],bins,mean)
  values <- numeric(settings$bins); occupied <- as.integer(names(means));values[occupied] <- means
  for (j in seq_len(length(occupied)-1L)) {
    a <- occupied[j]; b <- occupied[j+1L]
    if(b>a+1L && values[a]!=0 && values[b]!=0) {
      inside <- seq.int(a+1L,b-1L)
      values[inside] <- values[a]+(values[b]-values[a])*(centers[inside]-centers[a])/(centers[b]-centers[a])
    }
  }
  values
},numeric(settings$bins)))
for (statistic in c('max','min')) for (interval in list(c(0,5.3),c(1,2),c(2,3),c(3,5),c(5,5.3))) {
  selected <- times>=interval[1] & (times<interval[2] | (interval[2]==max(times) & times==max(times)))
  result <- tasks$investigate_cache_spectrum(cache,statistic,interval)
  same(result$data$intensity,apply(oracle[selected,,drop=FALSE],2,get(statistic)))
  same(result$scans,sum(selected))
}
empty <- tasks$investigate_cache_spectrum(cache,'min',c(3,4))
stopifnot(empty$scans==0L,all(is.na(empty$data$intensity)))
same(tasks$investigate_window_choices(cache,'full'),matrix(c(0,5.3),nrow=1))
for (width in c(20/60, 30/60, 1)) {
  windows <- tasks$investigate_window_choices(cache,width,c(1,5))
  stopifnot(all(abs(windows[,2]-windows[,1]-width)<1e-9),
            all(windows[,1]>=1),all(windows[,2]<=5))
}
interval <- c(1.19,4.41)
for (statistic in c('max','min')) {
  result <- tasks$investigate_cache_spectrum(cache,statistic,interval,raw=raw)
  selected <- times>=interval[1] & times<interval[2]
  same(result$data$intensity,apply(oracle[selected,,drop=FALSE],2,get(statistic)))
  same(result$scans,sum(selected))
}
stopifnot(inherits(try(tasks$investigate_cache_spectrum(cache,'max',c(1.19,1.21)),silent=TRUE),'try-error'))
same(diff(centers),rep(0.001,length(centers)-1L))
# Shared native axes must STILL use the fine grid, not native point spacing.
stopifnot(raw$sampling$shared == FALSE,settings$dmz==0.001)
native <- tasks$investigate_memory_raw(rep(scans[1],2),c(1.2,2),'Native.mz5')
native_cache <- tasks$investigate_cache_build(native,file.path(workspace,'Native'))
stopifnot(native$sampling$shared,native_cache$meta$settings$bins>nrow(scans[[1]]))
# Reopen without any spectra reader. A matching upload must not be read again.
no_read <- raw;no_read$read_peaks <- function(...)stop('Unexpected raw read')
reopened <- tasks$investigate_cache_build(no_read,dirname(cache$paths[1]))
same(reopened,cache)
stopifnot(is.null(tasks$investigate_cache_open(dirname(cache$paths[1]),native)))
changed <- raw;changed$fingerprint <- 'different-content-same-name'
stopifnot(is.null(tasks$investigate_cache_open(dirname(cache$paths[1]),changed)))
saved <- tasks$investigate_cache_open(dirname(cache$paths[1]))
stopifnot(!is.null(saved),is.null(saved$meta$raw$read_peaks),is.null(saved$meta$raw$handle))
same(tasks$investigate_cache_spectrum(saved,'min',c(0,5.3))$data$intensity,
  apply(oracle,2,min))
# Generation mismatch and interrupted/incomplete output must never be reused.
con <- dbConnect(SQLite(),cache$paths[2]);meta <- cache$meta;meta$complete <- FALSE
dbExecute(con,'UPDATE metadata SET payload=? WHERE id=1',params=list(list(serialize(meta,NULL))))
dbDisconnect(con)
stopifnot(is.null(tasks$investigate_cache_open(dirname(cache$paths[1]))))
cache <- tasks$investigate_cache_build(raw,dirname(cache$paths[1]))
stopifnot(!is.null(cache))
# A narrow viewport fetches fewer mass blocks than the full grid.
wide <- raw;wide$mz_range <- c(50,1000);wide$fingerprint <- paste0(raw$fingerprint,'-wide')
wide_cache <- tasks$investigate_cache_build(wide,file.path(workspace,'wide'))
full <- tasks$investigate_cache_spectrum(wide_cache,'max',c(0,5.3))
narrow <- tasks$investigate_cache_spectrum(wide_cache,'max',c(0,5.3),c(162,165))
stopifnot(narrow$blocks_read < full$blocks_read,nrow(narrow$data)<4000L)
same(narrow$data$intensity,full$data$intensity[match(narrow$data$mz,full$data$mz)])
cat('PASS: shared 0.001 grid; paired extrema; all scans; exact boundary ownership; partial/empty slices; zeros/internal interpolation; indexed region queries; source/settings validation; completed-cache reopen and interrupted-cache rejection.\n')

# Evaluate the actual Matabo section and Engine upload observer with an Engine
# reader fixture. This checks that Investigate never opens or reads a file.
code <- readLines('server.R',warn=FALSE)
first <- which(grepl('^  #####Matabo3.6 Investigate raw',code))
last <- which(grepl('^  #####Matabo3.6 End Investigate',code))
block <- parse(text=code[(first+1):(last-1)])
guard_names <- c('action_label','report_action_error','with_action_error','observe_action','observe_updates','plot_event_data')
signals <- new.env();signals$reads <- 0L;signals$autoruns <- 0L;signals$errors <- list();signals$updates <- list();signals$messages <- list()
upload_event <- Filter(function(x)is.call(x) && identical(x[[1]],as.name('observe_action')) &&
  identical(paste(deparse(x[[2]]),collapse=''),'input$mz5Files'),statements)[[1]]
harness <- function(input,output,session){
  for (name in guard_names) eval(assignment(name),environment())
  report_action_error <- function(action,error) {signals$errors[[length(signals$errors)+1L]] <- conditionMessage(error)}
  engine_loaded <- reactiveVal(list())
  uploadedmz5 <- reactiveVal(data.frame(FileName=character(),FilePath=character(),Source=character()))
  massData <- reactiveVal(data.frame(mz=162));refMassListData <- reactiveVal(data.frame(mz=162))
  analysis_mode <- function()if (identical(input$analysis_mode,'existing')) 'existing' else 'new'
  request_new_analysis <- function(){signals$autoruns <- signals$autoruns+1L}
  read_mz5_data <- function(...) {signals$reads <- signals$reads+1L;list(original=TRUE)}
  mz5_result_parent <- function(name)file.path(workspace,'Results',sub('\\.mz5$','',name))
  eval(block,environment())
  investigate_from_engine <- function(...)raw
  investigate_open_raw <- function(...)stop('Investigate opened a file')
  old_input <- session$sendInputMessage
  session$sendInputMessage <- function(id,message){signals$updates[[id]] <- message;old_input(id,message)}
  old_custom <- session$sendCustomMessage
  session$sendCustomMessage <- function(type,message){signals$messages[[length(signals$messages)+1L]] <- list(type=type,message=message);old_custom(type,message)}
  eval(upload_event,environment())
}
original_dir <- getwd();setwd(workspace)
testServer(harness,{
  session$setInputs(investigate_median_statistic='max',investigate_window_width='full',
    investigate_mz_min=162,investigate_mz_max=162.01,investigate_median_time_min=1.2,
    investigate_median_time_max=5.3,investigate_median_step=1,investigate_median_alignment_width=0.001,
    autorun=TRUE,mz5Files=data.frame(name='Fixture.mz5',datapath='one-engine-upload'))
  # The filename/status must be flushed before even the Engine reader starts.
  # This catches the reported "File Reading Complete" then blank picker stall.
  stopifnot(signals$reads==0L,signals$autoruns==0L,
    grepl('Fixture.mz5',signals$updates$investigate_engine_file$options,fixed=TRUE),
    grepl('Preparing Engine upload: Fixture.mz5',output$investigate_status,fixed=TRUE))
  for (attempt in 1:20) {
    later::run_now(0.1);session$flushReact()
    if (signals$autoruns==1L) break
  }
  stopifnot(signals$reads==1L,signals$autoruns==1L,identical(investigate_raw()$fingerprint,raw$fingerprint),
    !is.null(investigate_cache()),isTRUE(investigate_median()$cached))
  same(as.numeric(signals$updates$investigate_time_min$value),1.2)
  same(as.numeric(signals$updates$investigate_mz_group$value),40)
  before_plot <- output$investigate_median_plot
  generation <- investigate_plot_generation()
  session$setInputs(`plotly_relayout-metabomeister_median`='{"xaxis.range":[162.002,162.008],"yaxis.range":[-3,12]}')
  expected_view <- list(x=c(162.002,162.008),y=c(-3,12))
  same(investigate_median_view(),expected_view)
  session$setInputs(investigate_window_width='1')
  for(position in c(1L,3L,1L)) {
    session$setInputs(investigate_window_live=position)
    same(investigate_median_view(),expected_view)
    same(investigate_plot_generation(),generation)
  }
  same(output$investigate_median_plot,before_plot)
  same(investigate_median()$time_range,c(0,1))
  session$setInputs(investigate_extrema_start=1.19,investigate_extrema_end=4.41,
    investigate_window_width='full')
  session$setInputs(investigate_extrema_apply=1)
  same(investigate_median()$time_range,c(1.19,4.41))
  same(investigate_median()$scans,sum(times>=1.19 & times<4.41))
  same(investigate_median_view(),expected_view)
  stopifnot(signals$reads==1L)
  # The cache is already persisted: a second use of the same Engine upload
  # reuses both its raw object and completed summaries.
  engine_load_upload('one-engine-upload','Fixture.mz5')
  stopifnot(signals$reads==1L)
  session$setInputs(investigate_window_width='full',investigate_median_statistic='min')
  same(investigate_median()$statistic,'min')
  # Median keeps its own independent window/scan stride and native algorithm.
  session$setInputs(investigate_median_statistic='median',investigate_median_time_min=2,
    investigate_median_time_max=4.4,investigate_median_step=2,investigate_draw_median=1)
  same(investigate_median()$selected_times,c(2,2.9,4.4))
  stopifnot(!isTRUE(investigate_median()$cached),signals$reads==1L)
  session$setInputs(analysis_mode='existing',existing_run='Fixture',investigate_median_statistic='max')
  stopifnot(is.null(investigate_raw()$read_peaks),isTRUE(investigate_median()$cached),signals$reads==1L)
  same(investigate_median()$data$intensity,apply(oracle,2,max))
  # Reproduce files-on-disk/empty-GUI state and recover without an mz5 read.
  investigate_activate(NULL)
  session$setInputs(investigate_cache_run='Fixture',investigate_median_statistic='max',
    investigate_load_saved_spectra=1)
  stopifnot(isTRUE(investigate_median()$cached),signals$reads==1L,
    is.null(investigate_raw()$read_peaks))
  same(investigate_median()$data$intensity,apply(oracle,2,max))
  session$setInputs(investigate_median_statistic='min',investigate_load_saved_spectra=2)
  same(investigate_median()$data$intensity,apply(oracle,2,min))
  stopifnot(signals$reads==1L)
  # A new cache must not gate the filename/raw-data handshake. Run its deferred
  # callback with an intentional failure and verify that loaded data survives.
  assign('investigate_cache_open',function(...)NULL,envir=environment(engine_load_upload))
  assign('investigate_cache_build',function(...)stop('cache-build regression fixture'),envir=environment(engine_load_upload))
  engine_load_upload('one-engine-upload','Fixture.mz5')
  same(signals$updates$investigate_engine_file$selectize,NULL)
  stopifnot(grepl('Fixture.mz5',signals$updates$investigate_engine_file$options,fixed=TRUE),
    !is.null(investigate_raw()),is.null(investigate_cache()),signals$reads==1L)
  session$flushReact()
  for (attempt in 1:10) {
    later::run_now(0.1)
    session$flushReact()
    if (length(signals$errors)) break
  }
  stopifnot(!is.null(investigate_raw()),signals$reads==1L,
    identical(tail(signals$errors,1)[[1]],'cache-build regression fixture'))
  signals$errors <- head(signals$errors,-1L)
  # Reader errors still publish the uploaded name and an actionable status.
  investigate_activate(NULL)
  assign('read_mz5_data',function(...)stop('reader regression fixture'),envir=environment(engine_load_upload))
  session$setInputs(analysis_mode='new',mz5Files=data.frame(name='Broken.mz5',datapath='broken-upload'))
  stopifnot(grepl('Broken.mz5',signals$updates$investigate_engine_file$options,fixed=TRUE))
  for (attempt in 1:20) {
    later::run_now(0.1);session$flushReact()
    if (length(signals$errors)) break
  }
  stopifnot(grepl('Broken.mz5',signals$updates$investigate_engine_file$options,fixed=TRUE),
    grepl('Could not load Broken.mz5',output$investigate_status,fixed=TRUE),
    identical(tail(signals$errors,1)[[1]],'reader regression fixture'))
  signals$errors <- head(signals$errors,-1L)
})
stopifnot(!length(signals$errors))
setwd(original_dir)
proxy <- Filter(function(x)identical(x$type,'plotly-calls'),signals$messages)
stopifnot(length(proxy)>0)
# The curve-only live path emits restyle; only explicit reset may relayout.
stopifnot(all(vapply(proxy,function(x)x$message$method=='restyle',logical(1))))
unlink(workspace,recursive=TRUE)
cat('PASS: actual Engine upload reuse; automatic controls/status; live cached curves without widget rebuild or x/y zoom changes; independent Median; no extra raw-file read.\n')
