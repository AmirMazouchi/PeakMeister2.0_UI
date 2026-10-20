# Run from the application folder: Rscript --vanilla tests/saved-run-selection.R
# Tests actual Engine observers and saved-result selection with directory fixtures.
# Raw mz5 reading is forbidden and the scientific processing tail is excluded.
library(shiny)
invisible(parse('ui.R'))
expr<-parse('server.R')
is_assign<-function(x,name) is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(name))
is_event<-function(x,name) is.call(x)&&identical(x[[1]],as.name('observeEvent'))&&identical(paste(deparse(x[[2]]),collapse=''),paste0('input$',name))
text<-function(x) paste(deparse(x),collapse=' ')
server_expr<-Filter(function(x)is_assign(x,'server'),as.list(expr))[[1]]
statements<-as.list(server_expr[[3]][[3]])[-1]
get_assign<-function(name) Filter(function(x)is_assign(x,name),statements)[[1]]
get_event<-function(name) Filter(function(x)is_event(x,name),statements)[[1]]
first<-which(vapply(statements,is_assign,logical(1),name='empty_mz5_files'))
last<-which(vapply(statements,is_event,logical(1),name='initialize'))
engine<-statements[first:(last-1)]
engine<-Filter(function(x)!grepl('^output\\$fileTable',text(x)),engine)
extra_names<-c('plot_list_data_values','run_metadata','previous_metadata','modified_peak_plots','plot_render_trigger','get_db_path','db_exists','list_results_folders','selected_results_filename','results_folders_reactive','last_results_filename','latest_results_folder','last_visualization_filename','db_con','clear_loaded_analysis','peak_areas_data','migration_times_data')
extras<-lapply(extra_names,get_assign)
selector<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name('observeEvent'))&&grepl('^list\\(selected_results_filename\\(\\), results_folders_reactive\\(\\)\\)',text(x[[2]])),statements)
stopifnot(length(selector)==1)
clear_file<-Filter(function(x)is.call(x)&&identical(x[[1]],as.name('observeEvent'))&&identical(text(x[[2]]),'selected_results_filename()'),statements)
init<-get_event('initialize');parts<-as.list(init[[3]])[-1]
init[[3]]<-as.call(c(list(as.name('{')),parts[1:3],list(quote(events$analyses<-events$analyses+1L))))
block<-as.call(c(list(as.name('{')),engine,extras,selector,clear_file,list(get_event('results_folder2'),init)))
events<-new.env();events$updates<-list();events$notices<-character();events$nav<-NULL;events$analyses<-0L
showNotification<-function(ui,...) {events$notices<-c(events$notices,as.character(ui));invisible(NULL)}
updateSelectInput<-function(session,inputId,...) {events$updates[[inputId]]<-list(...);invisible(NULL)}
updateTabItems<-function(session,inputId,selected){events$nav<-selected;invisible(NULL)}
mode_server<-function(input,output,session){
 eval(block,envir=environment())
 reset_raw_selection<-function()NULL
 request_new_analysis<-function()NULL
}
readMSData<-function(...)stop('Existing analysis must never read raw data')
run_saved_selection_checks<-function(){
 work<-tempfile('saved-run-selection-');dir.create(work);old<-setwd(work)
 on.exit({setwd(old);unlink(work,recursive=TRUE)},add=TRUE)
 make_result<-function(folder,filename){
  path<-file.path(folder,'Data',filename);dir.create(path,recursive=TRUE,showWarnings=FALSE)
  placeholder<-list();save(placeholder,file=file.path(path,'plot_list_data.RData'))
 }
 make_result('./Results/Run01/Results 2026-10-19','Run01.mz5')
 make_result('./Results/Run01/Results 2026-10-20 2','Run01.mz5')
 make_result('./Results/Run01/Results 2026-10-20 10','Run01.mz5')
 make_result('./Results/Run02/Results only','Run02.mz5')
 make_result('./Run01/Results prior-layout','Run01.mz5')
 make_result('./Results legacy','Run01.mz5')
 make_result('./Results/Run01/Results wrong-file','Run02.mz5')
 make_result('./Results/Renamed/Results mismatch','Original.mz5')
 make_result('./Results/Space run/Results first','Space run.MZ5')
 dir.create('./Results/Run01/Results incomplete')
 dir.create('./Results/Empty')
 testServer(mode_server,{
  session$setInputs(analysis_mode='new')
  stopifnot(setequal(saved_runs(),c('Run01','Run02','Space run')))
  stopifnot(nrow(uploadedmz5())==0L,length(results_folders_reactive())==0L)
  for(bad in c('..','../Run01','Run01/../../outside'))stopifnot(length(saved_run_files(bad))==0L)
  stopifnot(identical(saved_run_files('Space run'),'Space run.MZ5'))
  session$setInputs(analysis_mode='existing',existing_run='Run01')
  stopifnot(identical(selected_saved_run(),'Run01'))
  stopifnot(nrow(uploadedmz5())==1L,uploadedmz5()$FileName=='Run01.mz5',
   uploadedmz5()$Source=='saved_run',is.na(uploadedmz5()$FilePath))
  stopifnot(setequal(results_folders_reactive(),c('./Results/Run01/Results 2026-10-19','./Results/Run01/Results 2026-10-20 2','./Results/Run01/Results 2026-10-20 10')))
  stopifnot(events$updates$results_folder$selected=='./Results/Run01/Results 2026-10-20 10')
  stopifnot(latest_results_folder(character())=='')
  session$setInputs(results_folder='./Results/Run01/Results 2026-10-19')
  make_result('./Results/Run01/Results 2026-10-21','Run01.mz5')
  session$elapse(5001)
  stopifnot(events$updates$results_folder$selected=='./Results/Run01/Results 2026-10-19')
  stopifnot(isTRUE(check_saved_analysis('./Results/Run01/Results 2026-10-19','Run01.mz5')))
  for(folder in c('./Results legacy','./Run01/Results prior-layout','./Results/Run02/Results only'))
   stopifnot(inherits(try(check_saved_analysis(folder,'Run01.mz5'),silent=TRUE),'try-error'))
  stopifnot(identical(output$show_visualization_file_selector,"false"))
  session$setInputs(initialize=1)
  stopifnot(events$analyses==0L)
  plot_list_data_values$file_name<-'old';previous_metadata$comment_df<-data.frame(a=1)
  session$setInputs(existing_run='Run02',file_selector='Run02.mz5')
  stopifnot(uploadedmz5()$FileName=='Run02.mz5',
   identical(results_folders_reactive(),'./Results/Run02/Results only'),
   events$updates$results_folder$selected=='./Results/Run02/Results only',
   is.null(plot_list_data_values$file_name),is.null(previous_metadata$comment_df))
  session$setInputs(clearFiles=1)
  stopifnot(nrow(uploadedmz5())==0L,length(selected_saved_run())==0L,length(results_folders_reactive())==0L)
  session$setInputs(existing_run='Missing')
  stopifnot(nrow(uploadedmz5())==0L,length(results_folders_reactive())==0L)
  session$setInputs(existing_run='Run01',file_selector='Run01.mz5')
  session$setInputs(analysis_mode='new')
  stopifnot(nrow(uploadedmz5())==0L,length(selected_saved_run())==0L,
   identical(output$show_visualization_file_selector,"true"))
  raw<-tempfile(fileext='.mz5');writeLines('test upload',raw)
  session$setInputs(mz5Files=data.frame(name='New.mz5',datapath=raw),autorun=FALSE)
  stopifnot(uploadedmz5()$Source=='upload',uploadedmz5()$FileName=='New.mz5')
  session$setInputs(initialize=2)
  stopifnot(events$analyses==1L)
  unlink(raw)
  session$setInputs(analysis_mode='existing',existing_run='Run02')
  unlink('./Results/Run02',recursive=TRUE)
  stopifnot(length(saved_run_files('Run02'))==0L,
   length(list_results_folders('Run02.mz5',selected_saved_run()))==0L)
 })
 # Missing Results root is a normal empty state.
 unlink('Results',recursive=TRUE)
 testServer(mode_server,{
  session$setInputs(analysis_mode='existing')
  stopifnot(length(saved_runs())==0L,length(results_folders_reactive())==0L)
 })
}
run_saved_selection_checks()
cat('PASS: run discovery; invalid/incomplete/mismatched exclusion; saved filenames without raw reads; per-run version filtering; switching/clearing; latest-by-date/numeric-version selection; Initialize guard; new-upload workflow; missing/deleted runs.\n')
