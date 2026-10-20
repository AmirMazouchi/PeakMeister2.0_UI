# Results are stored under Results/<mz5 filename without extension>/Results <date> [number].
# Visualization lists saved analyses matching its selected mz5 filename, including
# compatible legacy Results folders. Engine selects filenames only when reopening.
# Multi-file new analyses write a separate Results folder for each input file.

# Analysis modes added: New analysis reads uploaded mz5 data; Open existing
# analysis accepts filenames only and opens the corresponding saved results.
# The legacy pseudo-file workflow is retained. Source checks prevent name-only
# selections and pseudo files from reaching Initialize Run. Use with matching ui.R.

# Timing instrumentation added for Initialize Run (2026-10-16).
# Records elapsed seconds per input file for reading, calibration, extraction,
# smoothing, peak processing, plotting, RData saving, SQLite writes and exports.
# Shows stage messages in the progress panel and R console, and saves
# run_timings.csv in the Results folder after each stage, including total time
# and partial timings on errors or interrupts. Optional stages are labelled
# "if enabled"; their timings include checking whether that stage is enabled.
# The timing helper is embedded below: no additional R files are required.
# This adds measurement only; the scientific processing steps are unchanged.

#PeakMeister2.0_UI was written by Erick Helmeczi and Liam Surry 
#Initialized 2024-11-12 @McMaster University
#Department of Chemistry and Chemical Biology

#Set maximum file upload size. 
#Set maximum file upload size to 150 GB
options(shiny.maxRequestSize = 150 * 1024^3)
options(shiny.timeout = 600)  #Set timeout to 600 seconds (10 minutes)

#Initialize server
server <- function(input, output, session){ 
  
  #Initialize plotly_data as an empty reactive variable
  plotly_data <- reactiveVal(list())
  
  #Reactive values to store marker positions
  line_positions <- reactiveValues(red = 2, blue = 4)
  
  ####1. About tab####
  #Functions for "About" tab page -> reading information from markdown files
  output$disclaimerContent <- renderUI({
    includeMarkdown("docs/DISCLAIMER.md")
  })
  output$readmeContent <- renderUI({
    includeMarkdown(rmarkdown::render("docs/README.md"))
  })
  output$updatesContent <- renderUI({
    includeMarkdown("docs/UPDATES.md")
  })
  output$licenseContent <- renderUI({
    includeMarkdown("docs/LICENSE.md")
  })
####2. User supplied inputs tab####
 #Load data from user supplied Excel file using sheet names
  
  #####2.1 Initialize environment#####
  #Initialize reactive values as empty
  #Reactive values for mass list
  massData <- reactiveVal(data.frame(
    name = character(),
    mz = numeric(),
    extraction.window.ppm = numeric(),
    interference = character(),
    interference.comigration.threshold.seconds = numeric(),
    minimim.peak.width.seconds = numeric(),
    migration.window.seconds = numeric(),
    peak.space.tolerance.percent = numeric(),
    snr.threshold = numeric(),
    smoothing.kernel = character(),
    smoothing.strength = numeric(),
    stringsAsFactors = FALSE
  ))
  
  #Reactive values for reference mass list
  refMassListData <- reactiveVal(data.frame(
    name = character(),
    mz = numeric(),
    extraction.window.ppm = numeric(),
    class = character(),
    min.mt.min = numeric(),
    max.mt.min = numeric(),
    peak.fwhm.tolerance.multiplier = numeric(),
    peak.space.tolerance.percent = numeric(),
    smoothing.kernel = character(),
    smoothing.strength = numeric(),
    stringsAsFactors = FALSE
  ))
  
  #reactive values for paramaters
  parametersData <- reactiveVal(data.frame(
    number.of.injections = numeric(),
    ref.mass.one = numeric(),
    ref.mass.two = numeric(),
    ref.mass.window.ppm = numeric(),
    ref.mass.minimum.counts = numeric(),
    apply.mass.calibration = character(),
    apply.smoothing = character(),
    plot.format = character(),
    Manual.Indexes = character(),
    stringsAsFactors = FALSE
  ))
  #May need to include peak migfration times into these? Not sure???
  
  #Reactive values for manual MTI list
  userMTI <- reactiveVal(data.frame(
    name = character(),
    left_is = character(),
    right_is = character(),
    description = character(),
    stringsAsFactors = FALSE
  ))
  
  #Render editable tables
  output$massListData <- DT::renderDataTable({ massData() }, editable = TRUE)
  output$refMassListData <- DT::renderDataTable({ refMassListData() }, editable = TRUE)
  output$parametersData <- DT::renderDataTable({ parametersData() }, editable = TRUE)
  output$userMTIdata <- DT::renderDataTable({userMTIdata()}, editable = TRUE)
  
  #####2.2 Data Table Input#####
  #Look for an excel file/read parameters file.
  observeEvent(input$massList,
               {
                 req(input$massList)
                 #Load data from the uploaded file -> excel sheet must have sheet names "Mass List, Reference Mass List, and Paramaters"
                 #Mass List
                 mass_df <- readxl::read_excel(input$massList$datapath, sheet = "Mass List") %>%
                   as.data.frame()
                 #Reference Mass List
                 is_df <- readxl::read_excel(input$massList$datapath, sheet = "Reference Mass List") %>%
                   as.data.frame()
                 #Paramaters
                 parameters_df <- readxl::read_excel(input$massList$datapath, sheet = "Parameters") %>%
                   as.data.frame()
                 #Store the data in reactive values
                 massData(mass_df)
                 refMassListData(is_df)
                 parametersData(parameters_df)
                 })
  
  #Handle manual input when no file is uploaded -> allows users to input their own data manually 
  observeEvent(input$massList, {
    if (is.null(input$massList)) {
      massData(data.frame())  
      refMassListData(data.frame())  
      parametersData(data.frame())
    }
  })
  
  #####2.3 Manual Mass List Input#####
  #Add a row in the table 
  observeEvent(input$addMassRow, {
    newRow <- data.frame(
      #Input values for table 
      name = input$name,
      mz = input$mz,
      extraction.window.ppm = input$extraction.window.ppm,
      interference = input$interference,
      interference.comigration.threshold.seconds = input$interference.comigration.threshold.seconds,
      minimim.peak.width.seconds = input$minimim.peak.width.seconds,
      migration.window.seconds = input$migration.window.seconds,
      peak.space.tolerance.percent = input$peak.space.tolerance.percent,
      snr.threshold = input$snr.threshold,
      smoothing.kernel = input$smoothing.kernel,
      smoothing.strength = input$smoothing.strength,
      peak.1.mt = input$peak.1.mt,
      peak.2.mt = input$peak.2.mt,
      peak.3.mt = input$peak.3.mt,
      peak.4.mt = input$peak.4.mt,
      peak.5.mt = input$peak.5.mt,
      peak.6.mt = input$peak.6.mt,
      peak.7.mt = input$peak.7.mt,
      peak.8.mt = input$peak.8.mt,
      peak.9.mt = input$peak.9.mt,
      peak.10.mt = input$peak.10.mt,
      peak.11.mt = input$peak.11.mt,
      peak.12.mt = input$peak.12.mt,
      peak.13.mt = input$peak.13.mt,
      stringsAsFactors = FALSE
    )
    #Update dataset 
    updatedData <- rbind(massData(), newRow)
    massData(updatedData)
    
    #Clear inputs after adding allowing users to enter new ones
    updateTextInput(session, "name", value = "")
    updateNumericInput(session, "mz", value = 0)
    updateNumericInput(session, "extraction.window.ppm", value = 0)
    updateTextInput(session, "interference", value = "")
    updateNumericInput(session, "interference.comigration.threshold.seconds", value = 0)
    updateNumericInput(session, "minimim.peak.width.seconds", value = 0)
    updateNumericInput(session, "migration.window.seconds", value = 0)
    updateNumericInput(session, "peak.space.tolerance.percent", value = 0)
    updateNumericInput(session, "snr.threshold", value = 0)
    updateTextInput(session, "smoothing.kernel", value = "")
    updateNumericInput(session, "smoothing.strength", value = 0)
    
    #Clear peaks inputs
    lapply(1:13, function(i) {
      updateNumericInput(session, paste0("peak.", i, ".mt"), value = 0)
    })
  })
  
  #Delete selected row in Mass List server function 
  observeEvent(input$deletemassrow, {
    selected_row <- input$massListData_rows_selected
    if (length(selected_row) > 0) {
      updatedData <- massData()[-selected_row, ]
      massData(updatedData)} else {
      showNotification("Please select a row to delete.", type = "error")}})
  #Server function for clearing data table
  observeEvent(input$clearMassList, {
    massData(data.frame())
    showNotification("Table has been reset.", type = "message")})
  #Data table output
  output$massListData <- DT::renderDataTable({
    DT::datatable(massData(), selection = 'single', editable = TRUE)})
  
  #####2.4 Manual Reference Mass Input#####
  observeEvent(input$addRefMassRow, {
    newRow <- data.frame(
      #Input values for table
      name = input$name,
      mz = input$mz,
      extraction.window.ppm = input$extraction.window.ppm,
      class = input$class,
      min.mt.min = input$min.mt.min,
      max.mt.min = input$max.mt.min,
      peak.fwhm.tolerance.multiplier = input$peak.fwhm.tolerance.multiplier,
      peak.space.tolerance.percent = input$peak.space.tolerance.percent,
      smoothing.kernel = input$smoothing.kernel,
      smoothing.strength = input$smoothing.strength,
      peak.1.mt = input$peak.1.mt,
      peak.2.mt = input$peak.2.mt,
      peak.3.mt = input$peak.3.mt,
      peak.4.mt = input$peak.4.mt,
      peak.5.mt = input$peak.5.mt,
      peak.6.mt = input$peak.6.mt,
      peak.7.mt = input$peak.7.mt,
      peak.8.mt = input$peak.8.mt,
      peak.9.mt = input$peak.9.mt,
      peak.10.mt = input$peak.10.mt,
      peak.11.mt = input$peak.11.mt,
      peak.12.mt = input$peak.12.mt,
      peak.13.mt = input$peak.13.mt,
      stringsAsFactors = FALSE
    )
    #Update dataset
    updatedData <- rbind(refMassListData(), newRow)
    refMassListData(updatedData)
    
    #Clear inputs after adding
    updateTextInput(session, "name", value = "")
    updateNumericInput(session, "mz", value = 0)
    updateNumericInput(session, "extraction.window.ppm", value = 0)
    updateTextInput(session, "class", value = "")
    updateNumericInput(session, "min.mt.min", value = 0)
    updateNumericInput(session, "max.mt.min", value = 0)
    updateNumericInput(session, "peak.fwhm.tolerance.multiplier", value = 0)
    updateNumericInput(session, "peak.space.tolerance.percent", value = 0)
    updateTextInput(session, "smoothing.kernel", value = "")
    updateNumericInput(session, "smoothing.strength", value = 0)
    
    #Clear peaks inputs
    lapply(1:13, function(i) {
      updateNumericInput(session, paste0("ref.peak.", i, ".mt"), value = 0)
    })
  })
  
  #Deleting a selected row in reference mass list
  observeEvent(input$deleterefmassrow, {
    selected_row <- input$refMassListData_rows_selected
    if (length(selected_row) > 0) {
      updatedData <- refMassListData()[-selected_row, ]
      refMassListData(updatedData)
    } else {
      showNotification("Please select a row to delete.", type = "error")}})
  #Clearing the reference mass list
  observeEvent(input$clearRefMassList, {
    refMassListData(data.frame())
    showNotification("Reference Mass List has been cleared.", type = "message")})
  #Data table output
  output$refMassListData <- DT::renderDataTable({
    DT::datatable(refMassListData(), selection = 'single', editable = TRUE)})
  
  #####2.5 Manual Parameters Input#####
  observeEvent(input$addParamRow, {
    newRow <- data.frame(
      #Input values for table
      number.of.injections = input$number.of.injections,
      ref.mass.one = input$ref.mass.one,
      ref.mass.two = input$ref.mass.two,
      ref.mass.window.ppm = input$ref.mass.window.ppm,
      ref.mass.minimum.counts = input$ref.mass.minimum.counts,
      apply.mass.calibration = input$apply.mass.calibration,
      apply.smoothing = input$apply.smoothing,
      plot.format = input$plot.format,
      Manual.Indexes = input$Manual.Indexes,
      stringsAsFactors = FALSE
    )
    #Update dataset
    updatedData <- rbind(parametersData(), newRow)
    parametersData(updatedData)
    
    #Clear input values when user adds more paramaters 
    updateNumericInput(session, "number.of.injections", value = 0)
    updateNumericInput(session, "ref.mass.one", value = 0)
    updateNumericInput(session, "ref.mass.two", value = 0)
    updateNumericInput(session, "ref.mass.window.ppm", value = 0)
    updateNumericInput(session, "ref.mass.minimum.counts", value = 0)
    updateCheckboxInput(session, "apply.mass.calibration", value = FALSE)
    updateCheckboxInput(session, "apply.smoothing", value = FALSE)
    updateTextInput(session, "plot.format", value = "")
    updateCheckboxInput(session, "Manual.Indexes", value = FALSE)
  })
  
  #Deleting selected parameter row
  observeEvent(input$deleteparamrow, {
    selected_row <- input$parametersData_rows_selected
    if (length(selected_row) > 0) {
      updatedData <- parametersData()[-selected_row, ]
      parametersData(updatedData)
    } else {
      showNotification("Please select a row to delete.", type = "error")}})
  #Clearing the parameters list
  observeEvent(input$clearparam, {
    parametersData(data.frame())
    showNotification("Parameters List has been cleared.", type = "message")})
  #Data table output
  output$parametersData <- DT::renderDataTable({
    DT::datatable(parametersData(), selection = 'single', editable = TRUE)})
  
  #####2.6 Project information input#####
  #The idea behind this tab is to automatically generate a metadata file on the project information that is included with the analysis. Makes environment more reproducible and traceable. Additionally, project information will be saved into the generated results folder so results are easier to identify/organize on your computer. 
  
  #Generate reactive values for information
  projectName <- reactiveVal("")
  projectDescription <- reactiveVal("")
  projectSupervisor <- reactiveVal("")
  projectContact <- reactiveVal("")
  #Update reactive values when someone clicks Save project info)
  observeEvent(input$saveProjectInfo, {
    projectName(input$projectName)  
    projectDescription(input$projectDescription)  
    projectSupervisor(input$projectSupervisor)  
    projectContact(input$projectContact)  
  })
  #Display output information from reactive variables
  output$currentProjectInfo <- renderText({
    paste("Project Name:", projectName(), 
          "\nProject Description:", projectDescription(), 
          "\nProject Supervisor:", projectSupervisor(), 
          "\nProject Contact:", projectContact())
  })
  
  #Server function for clearing project information when button is clicked
  observeEvent(input$clearproject, {
    #Clear the input fields
    updateTextInput(session, "projectName", value = "")
    updateTextInput(session, "projectDescription", value = "")
    updateTextInput(session, "projectSupervisor", value = "")
    updateTextInput(session, "projectContact", value = "")
    #Reset reactive values
    projectName("")
    projectDescription("")
    projectSupervisor("")
    projectContact("")
    showNotification("Project information has been cleared.", type = "message")
  })
  
  #####2.7 Manual MTI input#####
  
  #Upload excel sheet titled user supplied migration indexes
  observeEvent(input$userMTI, {
    req(input$userMTI)
    userMTI_df <- readxl::read_excel(input$userMTI$datapath, sheet = "User Supplied Migration Indexes") %>% as.data.frame()
    
    #Store file in reactive value
    userMTI(userMTI_df)
  })
  
  #Manual Data Input
  observeEvent(input$addMTIrow, {
    newRow <- data.frame(
      name = input$name,
      left_is = input$left_is,
      right_is = input$right_is,
      description = input$description,
      `1` = input$peak.1.mt,
      `2` = input$peak.2.mt,
      `3` = input$peak.3.mt,
      `4` = input$peak.4.mt,
      `5` = input$peak.5.mt,
      `6` = input$peak.6.mt,
      `7` = input$peak.7.mt,
      `8` = input$peak.8.mt,
      `9` = input$peak.9.mt,
      `10` = input$peak.10.mt,
      `11` = input$peak.11.mt,
      `12` = input$peak.12.mt,
      `13` = input$peak.13.mt,
      stringsAsFactors = FALSE
    )
    #Update dataset
    updatedData <- rbind(userMTI(), newRow)
    userMTI(updatedData)
    
    #Clear inputs values when user adds more things
    updateTextInput(session, "name", value = "")
    updateTextInput(session, "left_is", value = "")
    updateTextInput(session, "right_is", value = "")
    updateTextInput(session, "description", value = "")
    lapply(1:13, function(i) {
      updateNumericInput(session, paste0("ref.peak.", i, ".mt"), value = 0)
    })
  })
  
  #Delete selected MTI row
  observeEvent(input$deleteMTIrow, {
    selected_row <- input$userMTIdata_rows_selected
    if (length(selected_row) > 0) {
      updatedData <- userMTI()[-selected_row, ]
      userMTI(updatedData)
    } else {
      showNotification("Please select a row to delete.", type = "error")
    }
  })
  
  #Clearing the MTI list
  observeEvent(input$clearMTIlist, {
    userMTI(data.frame())
    showNotification("User Supplied Migration Indexes list has been cleared.", type = "message")
  })
  
  #Data table output for userMTI information
  output$userMTIdata <- DT::renderDataTable({
    DT::datatable(userMTI(), selection = 'single', editable = TRUE)})
  
  
  #####2.8 Save input information#####
  #Function for saving datatable once user has uploaded their values
  output$downloadData <- downloadHandler(
    #Save data datble as excel sheet with name "Mass List and Parameters"
    #Will add function to save it with a project identifier as well for organizational purposes
    filename = function() {
      paste("Mass List and Parameters", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      #Create a list of data frames for each sheet
      write_xlsx(list(
        Parameters = parametersData(),
        'Mass List' = massData(),
        'Reference Mass List' = refMassListData()
      ), path = file)
    }
  )
  
  ####3. Engine tab####
  #####3.1 mz5 file upload#####
  ######3.1.1 Upload real .mz5 files######
  empty_mz5_files <- function() {
    data.frame(FileName = character(), FilePath = character(), Source = character(),
               stringsAsFactors = FALSE)
  }
  uploadedmz5 <- reactiveVal(empty_mz5_files())
  pseudomz5names <- reactiveVal(data.frame(FileName = "", stringsAsFactors = FALSE))
  existing_mz5_names <- reactiveVal(character())
  existing_analysis_status <- reactiveVal("")
  analysis_mode <- reactive({
    if (identical(input$analysis_mode, "existing")) "existing" else "new"
  })

  reset_raw_selection <- function() shinyjs::reset("mz5Files")
  request_new_analysis <- function() shinyjs::click("initialize")

  # Reject ambiguous names before constructing paths or matching stored records.
  valid_mz5_names <- function(names) {
    is.character(names) && length(names) > 0L && !anyNA(names) &&
      all(nzchar(names)) && !anyDuplicated(names) &&
      all(grepl("\\.mz5$", names, ignore.case = TRUE)) &&
      !any(grepl("[/\\\\]", names))
  }

  mz5_result_parent <- function(filename) {
    if (length(filename) != 1L || !valid_mz5_names(filename)) stop("Select a valid mz5 filename.")
    stem <- sub("\\.mz5$", "", filename, ignore.case = TRUE)
    if (stem %in% c("", ".", "..")) stop("The mz5 filename must have a nonempty name.")
    file.path(".", "Results", stem)
  }

  create_results_folder <- function(filename) {
    parent <- mz5_result_parent(filename)
    if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE)) {
      stop("Could not create the results parent folder: ", parent)
    }
    count <- 1L
    repeat {
      leaf <- paste("Results", Sys.Date(), if (count == 1L) "" else count)
      folder <- file.path(parent, trimws(leaf))
      if (!file.exists(folder) && dir.create(folder, showWarnings = FALSE)) return(folder)
      if (!file.exists(folder)) stop("Could not create the Results folder: ", folder)
      count <- count + 1L
    }
  }

  prepare_results_folder <- function(folder, masses, references, parameters, indexes = NULL) {
    for (name in c("Plots", "Data")) dir.create(file.path(folder, name), showWarnings = FALSE)
    if (parameters$plot.format == "Metabolite") {
      for (name in c(references$name, masses$name)) {
        dir.create(file.path(folder, "Plots", name), showWarnings = FALSE)
      }
    }
    write_xlsx(list(Parameters = parameters, 'Mass List' = masses,
                    'Reference Mass List' = references),
               path = file.path(folder, paste0("Mass List and Parameters", Sys.Date(), ".xlsx")))
    if (!is.null(indexes)) {
      write.csv(indexes, file.path(folder, "Migration Index Summary.csv"), row.names = FALSE)
    }
  }

  # Only inspect saved result files/database; never open the selected raw files.
  check_saved_analysis <- function(folder, names) {
    if (!valid_mz5_names(names)) stop("Select real .mz5 files with distinct filenames.")
    if (length(folder) != 1L || is.na(folder) || !folder %in% list_results_folders(names[1])) {
      stop("Select an available saved Results folder.")
    }
    paths <- file.path(folder, "Data", names, "plot_list_data.RData")
    missing <- names[!file.exists(paths)]
    if (length(missing)) {
      stop("No saved results found for: ", paste(missing, collapse = ", "),
           ". Select the matching Results folder and original filenames.")
    }
    # If no database exists, the existing migration dialog will offer to build it.
    if (db_exists(folder)) {
      con <- DBI::dbConnect(RSQLite::SQLite(), get_db_path(folder), flags = RSQLite::SQLITE_RO)
      on.exit(DBI::dbDisconnect(con), add = TRUE)
      stored <- DBI::dbGetQuery(con, "SELECT DISTINCT file_name FROM plot_axis_data")$file_name
      missing <- names[!sub("\\.mz5$", "", names) %in% stored]
      if (length(missing)) stop("Files are missing from the saved database: ", paste(missing, collapse = ", "))
    }
    invisible(TRUE)
  }

  output$fileTable <- DT::renderDataTable({
    DT::datatable(uploadedmz5()[, "FileName", drop = FALSE], options = list(pageLength = 5))
  })
  output$existing_selection_summary <- renderText({
    names <- existing_mz5_names()
    if (!length(names)) "No files selected." else paste(names, collapse = ", ")
  })
  output$existing_analysis_status <- renderText(existing_analysis_status())

  observeEvent(input$analysis_mode, {
    uploadedmz5(empty_mz5_files())
    existing_mz5_names(character())
    existing_analysis_status("")
    reset_raw_selection()
  }, ignoreInit = TRUE, priority = 100)

  observeEvent(input$mz5Files, {
    if (analysis_mode() != "new") return()
    req(input$mz5Files)
    newFiles <- data.frame(FileName = input$mz5Files$name,
                           FilePath = input$mz5Files$datapath, Source = "upload",
                           stringsAsFactors = FALSE)
    # Preserve uploads from this mode only; placeholders cannot become raw inputs.
    current <- uploadedmz5()
    current <- current[current$Source == "upload", , drop = FALSE]
    uploadedmz5(rbind(current, newFiles))
    if (isTRUE(input$autorun)) request_new_analysis()
  })

  observeEvent(input$existing_mz5_selection, {
    if (analysis_mode() != "existing") return()
    names <- unlist(input$existing_mz5_selection$names, use.names = FALSE)
    uploadedmz5(empty_mz5_files())
    existing_mz5_names(character())
    existing_analysis_status("")
    if (!valid_mz5_names(names)) {
      showNotification("Select real .mz5 files with distinct filenames.", type = "error")
      return()
    }
    existing_mz5_names(names)
    uploadedmz5(data.frame(FileName = names, FilePath = NA_character_,
                           Source = "name_only", stringsAsFactors = FALSE))
    existing_analysis_status("Filenames selected. In Visualization, select a file and one of its saved analyses.")
  })

  observeEvent(input$open_existing_analysis, {
    if (analysis_mode() != "existing") return()
    if (!nrow(uploadedmz5())) {
      showNotification("Select real .mz5 files first (names only).", type = "error")
      return()
    }
    updateTabItems(session, "main_navigation", selected = "processing")
  })

  observeEvent(input$clearFiles, {
    uploadedmz5(empty_mz5_files())
    pseudomz5names(data.frame(FileName = "", stringsAsFactors = FALSE))
    existing_mz5_names(character())
    existing_analysis_status("")
    reset_raw_selection()
  })

  ######3.1.2 Upload file names and generate pseudo .mz5 files######
  #Render the editable data table
  output$mz5fileTable <- DT::renderDataTable({
    datatable(pseudomz5names(), editable = TRUE, options = list(pageLength = 5))
  })
  
  #Observe changes in the table and update the reactive variable with the input information
  observeEvent(input$mz5fileTable_cell_edit, {
    info <- input$mz5fileTable_cell_edit
    newFileNames <- pseudomz5names()
    newFileNames[info$row, info$col] <- info$value
    pseudomz5names(newFileNames)})
  
  #Update the reactive value with filenames from the uploaded Excel file
  observeEvent(input$excelFile, {
    req(input$excelFile)
    excel_data <- read_excel(input$excelFile$datapath)
    newFiles <- data.frame(
      FileName = excel_data$FileName,
      stringsAsFactors = FALSE)
    pseudomz5names(newFiles)
  })
  
  #Generate pseudo .mz5 files based on filenames in the data table. 
  #These are empty .mz5 files soley used for accessing the metadata specific to each file without having to upload these large files every single time to process the data.
  observeEvent(input$generateFiles, {
    req(pseudomz5names())
    
    #Generate empty .mz5 files
    newFiles <- data.frame(
      FileName = paste0(pseudomz5names()$FileName, ".mz5"),
      FilePath = sapply(pseudomz5names()$FileName, function(name) {
        file_path <- tempfile(fileext = ".mz5")
        file.create(file_path)
        file_path
      }), Source = "pseudo", stringsAsFactors = FALSE)
    
    #Store empty .mz5 files in the reactive variable for use elsewhere in the app.
    uploadedmz5(newFiles)
    
    showNotification("Pseudo .mz5 files generated", type = "message")
  })
  
  #Add empty rows to the table
  observeEvent(input$addRow, {
    currentData <- pseudomz5names()
    newRow <- data.frame(FileName = "", stringsAsFactors = FALSE)
    updatedData <- rbind(currentData, newRow)
    pseudomz5names(updatedData)
  })
  
  
  #####3.2 Initialze Run Button - Main Content######
  #Initialize run button - main engine content 
  observeEvent(input$initialize, {
    if (analysis_mode() != "new") {
      shiny::showNotification("Choose New analysis to process raw files.", type = "error")
      return(NULL)
    }
    files <- uploadedmz5()
    if (!nrow(files) || !valid_mz5_names(files$FileName) || any(files$Source != "upload") ||
        anyNA(files$FilePath) || !all(file.exists(files$FilePath))) {
      shiny::showNotification("Upload real .mz5 files for a new analysis. Name-only selections and pseudo files cannot be processed.", type = "error")
      return(NULL)
    }
    
    # Prepare the first file's output; subsequent files get their own folders.
    file_name <- create_results_folder(files$FileName[1])

    # Stage-level timings use elapsed time, including file I/O and database waits.
    # The CSV is refreshed at every stage boundary so completed stages survive errors.
    new_run_timer <- function(path, notify = function(stage, file) NULL) {
      rows <- data.frame(file = character(), stage = character(), status = character(),
                         started_at = character(), elapsed_seconds = numeric(),
                         detail = character(), stringsAsFactors = FALSE)
      run_start <- proc.time()[["elapsed"]]
      run_wall <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
      active <- NULL
      finished <- FALSE
      persist <- function() write.csv(rows, path, row.names = FALSE)
      record <- function(file, stage, status, wall, elapsed, detail = "") {
        rows[nrow(rows) + 1L, ] <<- list(file, stage, status, wall, elapsed, detail)
        persist()
        message(sprintf("[Timing] %s | %s | %s | %.3f seconds", file, stage, status, elapsed))
      }
      close_stage <- function(status = "completed", detail = "") {
        if (!is.null(active)) {
          record(active$file, active$stage, status, active$wall,
                 proc.time()[["elapsed"]] - active$start, detail)
          active <<- NULL
        }
      }
      persist()
      list(
        stage = function(stage, file = "") {
          stopifnot(!finished)
          close_stage()
          active <<- list(file = file, stage = stage,
                          wall = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                          start = proc.time()[["elapsed"]])
          message(sprintf("[Timing] Starting %s%s", stage,
                          if (nzchar(file)) paste0(" | ", file) else ""))
          notify(stage, file)
        },
        end_stage = function() close_stage(),
        finish = function(status = "completed", detail = "") {
          if (finished) return(invisible(NULL))
          close_stage(status, detail)
          record("", "TOTAL RUN (includes all stages)", status, run_wall,
                 proc.time()[["elapsed"]] - run_start, detail)
          finished <<- TRUE
          invisible(rows)
        }
      )
    }
    timing <- new_run_timer(file.path(file_name, "run_timings.csv"),
      notify = function(stage, file) {
        if (nzchar(file)) shiny::setProgress(detail = stage)
      })
    tryCatch({
    timing$stage("Validate parameters and calculate migration indexes")

    #Create a vector of metabolite and internal standard names
    name_vec <- c(refMassListData()$name, massData()$name)
    
    #Check for duplicate metabolite and internal standard names
    if (sum(duplicated(name_vec)) > 0) {
      stop("Mass List and Parameters contain duplicated metabolite or internal standard names.")
    }
    
    #Determine the number of metabolites and internal standards
    num_of_metabolites <- nrow(massData())
    num_of_is <- nrow(refMassListData())
    
    #Get the number of injections
    num_of_injections <- parametersData()$number.of.injections[1]
    
    prepare_results_folder(file_name, massData(), refMassListData(), parametersData())

    #####3.3 Migration Index Calculation#####
    
    #Summarize user supplied migration time data
    metabolites_mt_df <- data.frame(name = massData()$name, massData()[, c((ncol(massData()) - num_of_injections + 1):ncol(massData()))])
    is_mt_df <- subset(refMassListData(), refMassListData()$class == "Reference")
    is_mt_df <- data.frame(name = is_mt_df$name, is_mt_df[, c((ncol(is_mt_df) - num_of_injections + 1):ncol(is_mt_df))])
    
    #Determine IS on the left
    is_left_vec <- character(num_of_metabolites)
    
    for (m in 1:nrow(metabolites_mt_df)) {
      is_temp <- (metabolites_mt_df[m, 2] - (is_mt_df[, 2]))
      is_temp <- set_names(is_temp, is_mt_df$name)
      is_temp <- is_temp[which(sign(is_temp) == 1 | sign(is_temp) == 0)] %>%
        which.min() %>%
        names()
      
      if (length(is_temp) == 0) {
        is_temp <- "none"
      }
      
      is_left_vec[m] <- is_temp
    }
    
    #Determine IS on the right
    is_right_vec <- character(num_of_metabolites)
    
    for (m in 1:nrow(metabolites_mt_df)) {
      is_temp <- (metabolites_mt_df[m, 2] - (is_mt_df[, 2]))
      is_temp <- set_names(is_temp, is_mt_df$name)
      is_temp <- is_temp[which(sign(is_temp) == -1 | sign(is_temp) == 0)] %>%
        which.max() %>%
        names()
      
      if (length(is_temp) == 0) {
        is_temp <- "none"
      }
      
      is_right_vec[m] <- is_temp
    }
    
    #Compute migration index
    mi_df <- matrix(NA, ncol = (num_of_injections + 1), nrow = nrow(metabolites_mt_df)) %>%
      as.data.frame()
    mi_df[, 1] <- metabolites_mt_df$name
    
    summary_vec <- character(nrow(metabolites_mt_df))
    
    for (m in 1:nrow(metabolites_mt_df)) {
      for (i in 2:(num_of_injections + 1)) {
        left <- is_left_vec[m]
        right <- is_right_vec[m]
        
        if (left == "none") {
          mi_df[m, i] <- metabolites_mt_df[m, i] / is_mt_df[which(is_mt_df$name == right), i]
          summary_vec[m] <- right
          next
        }
        
        if (right == "none") {
          mi_df[m, i] <- metabolites_mt_df[m, i] / is_mt_df[which(is_mt_df$name == left), i]
          summary_vec[m] <- left
          next
        }
        
        if (metabolites_mt_df[m, 2] / is_mt_df[which(is_mt_df$name == left), 2] < 1.01) {
          mi_df[m, i] <- metabolites_mt_df[m, i] / is_mt_df[which(is_mt_df$name == left), i]
          summary_vec[m] <- left
          is_right_vec[m] <- "none"
          next
        }
        
        if (metabolites_mt_df[m, 2] / is_mt_df[which(is_mt_df$name == right), 2] > 0.99) {
          mi_df[m, i] <- metabolites_mt_df[m, i] / is_mt_df[which(is_mt_df$name == right), i]
          summary_vec[m] <- right
          is_left_vec[m] <- "none"
          next
        }
        
        mi_df[m, i] <- (metabolites_mt_df[m, i] - is_mt_df[which(is_mt_df$name == left), i]) / 
          (is_mt_df[which(is_mt_df$name == right), i] - is_mt_df[which(is_mt_df$name == left), i])
        summary_vec[m] <- "mi"
      }
    }
    
    #Store results in one data frame
    colnames(mi_df)[2:ncol(mi_df)] <- c(1:num_of_injections)
    mi_df <- cbind(name = massData()$name,
                   left_is = is_left_vec,
                   right_is = is_right_vec,
                   description = summary_vec,
                   mi_df[2:ncol(mi_df)])
    #####3.4 User Supplied MTIs#####
    
    user_mti_df <- userMTI()
    
    if (parametersData()$Manual.Indexes == "Yes") {
        
        names <- unique(c(user_mti_df$name, user_mti_df$left_is, user_mti_df$right_is, user_mti_df$description))
        
        if (all(names %in% c(name_vec, "mi", "none")) == FALSE) {
          stop("Names detected in User Supplied Migration Indexes.csv not found in Mass List and Parameters.")
        }
        
        if (ncol(user_mti_df) != (4 + num_of_injections)) {
          stop(paste("Make sure User Supplied Migration Indexes.csv contains migration time indexes for all", num_of_injections, "injections."))
        }
        
        for (i in 1:nrow(user_mti_df)) {
          row <- which(mi_df$name == user_mti_df$name[i])
          
          if (length(row) == 0) {
            stop(paste("Metabolite", user_mti_df$name[i], "from User Supplied Migration Indexes.csv is not in Mass List and Parameters."))
          }
          
          mi_df[row, 2:ncol(mi_df)] <- user_mti_df[i, 2:ncol(user_mti_df)]
        }
      }
    
    #Write the migration index summary to a CSV file
    write.csv(mi_df, 
              file.path(file_name, "Migration Index Summary.csv"), 
              row.names = FALSE)
    
    #####3.5 Data File Analysis#####
    #Create vectors for data file directories and name
    data_files <- uploadedmz5()$FilePath
    data_file_names <- uploadedmz5()$FileName  
    plotly_objects <- list()
    
    
    #Define the function to save plotly_objects to .RDA file
    save_plotly_objects <- function(file_name) {
      save(plotly_objects, file = file.path(file_name, "plotly_objects.RData"))
    }
    
    
    
    for (d in seq_along(data_files)){
      if (d > 1L) {
        file_name <- create_results_folder(data_file_names[d])
        timing <- new_run_timer(file.path(file_name, "run_timings.csv"),
          notify = function(stage, file) {
            if (nzchar(file)) shiny::setProgress(detail = stage)
          })
        timing$stage("Prepare output folders")
        prepare_results_folder(file_name, massData(), refMassListData(), parametersData(), mi_df)
      }
      # Each folder contains only this file's reports and plot objects.
      plotly_objects <- list()
      #Initialize progress bar
      withProgress(message = paste("Processing data file:", data_file_names[d]), value = 0,{ total_steps <- 7
 
      
      print(paste(d, ". ", "Analyzing Data File: ", data_file_names[d], sep = ""))
      
      #Create a subfolder for the current file inside the Data folder
      subfolder_path <- file.path(file_name, "Data", data_file_names[d])
      dir.create(path = subfolder_path, showWarnings = FALSE, recursive = TRUE)
      
      
      timing$stage("Copy mz5 file", data_file_names[d])
      #Make a copy of the data file as data will be written directly to this file during mass calibration
      file <- gsub(".mz5", "", data_files[d])
      file.copy(data_files[d], to = paste(file, "temp.mz5", sep = "_"))
      
      timing$stage("Read mz5 into memory", data_file_names[d])
      #Read copied data file and update progress bar
      incProgress(1/total_steps, detail = paste("Reading Data File"))
      print("Reading Data File")
      
      run_data <- readMSData(
        file = paste(file, "temp.mz5", sep = "_"),
        pdata = NULL,
        msLevel = 1,
        verbose = isMSnbaseVerbose(),
        centroided. = FALSE,
        smoothed. = FALSE,
        cache. = 0,
        mode = "inMemory"
      )
      
      print("File Reading Complete")
      
      #Unlock "assayData" environment 
      env_binding_unlock(run_data@assayData)
      
      timing$stage("Mass calibration (if enabled)", data_file_names[d])
      #####3.6 Perform Mass Calibration#####
      
      #This is a stupid line but not sure way around it currently
      mass_df <- massData()
      is_df <- refMassListData()
      parameters_df <- parametersData()
      

      calibration_response <- parameters_df$apply.mass.calibration
      
      #Confirm response is "Yes" or "No". Otherwise, produce an error.
      
      if(calibration_response != "Yes" & calibration_response != "No"){
        stop(paste("apply.mass.calibration parameter should be 'Yes' or 'No'. ", "'", calibration_response, "' is not an acceptable input.", sep = ""))
      }
      
      #Check if mass calibration should be applied
      
      if (calibration_response == "No"){
        
        print("Skipping Mass Calibration")
        
      }
      
      if (calibration_response == "Yes"){
        
        #Update progress bar
        incProgress(1/total_steps, detail = paste("Performing Mass Calibration"))
        
        print("Performing Mass Calibration")
        
        run_data <- local({
          
          #Define mass window and minimum lock mass counts
          
          mass_window <- parameters_df$ref.mass.window.ppm[1]
          
          minimum_counts <- parameters_df$ref.mass.minimum.counts[1]
          
          #Define a function to calculate 1. experimental m/z of lock masses, 2. the corresponding mass error, and 3. the index of the lock mass
          #The lock mass value is currently the most intense point in the mass window plus an adjustment to better predict the lock mass value.
          
          calibration_parameters <- function(spectrum, lock_mass, minimum_counts, mass_window) {
            
            lock_mass_range <- c((lock_mass - lock_mass * mass_window / 1000000), 
                                 (lock_mass + lock_mass * mass_window / 1000000))
            
            #Find the index of the most intense point in the lock mass window
            
            max_intensity_index <- spectrum %>%
              filter(mz >= lock_mass_range[1] & mz <= lock_mass_range[2]) %>%
              slice_max(intensity) %>%
              pull(index)
            
            if(length(max_intensity_index) == 1){
              if(spectrum$intensity[max_intensity_index] > minimum_counts){
                
                #Create a data frame with all four points
                line_points <- data.frame("mz" = c(spectrum$mz[(max_intensity_index - 2):(max_intensity_index - 1)], 
                                                   spectrum$mz[(max_intensity_index + 1):(max_intensity_index + 2)]),
                                          "intensity" = c(spectrum$intensity[(max_intensity_index - 2):(max_intensity_index - 1)], 
                                                          spectrum$intensity[(max_intensity_index + 1):(max_intensity_index + 2)]))
                
                #Create a model for all four points
                model_left <- lm(formula = intensity ~ mz, data = line_points[c(1,2),])
                model_right <- lm(formula = intensity ~ mz, data = line_points[c(3,4),])
                
                #Get the coefficients for both models
                coefficients_left <- c(model_left$coefficients["mz"], model_left$coefficients["(Intercept)"])
                coefficients_right <- c(model_right$coefficients["mz"], model_right$coefficients["(Intercept)"])
                
                #Calculate the slope and intercept
                slope <- coefficients_left[1] - coefficients_right[1]
                intercept <- coefficients_right[2] - coefficients_left[2]
                
                #Solve for the experimental_mz
                experimental_mz <- solve(slope, intercept)
                
                experimental_mass_diff <- lock_mass - experimental_mz
                
                return(c(experimental_mz, experimental_mass_diff, max_intensity_index))
              }
            }
          }
          
          #Loop through each spectrum to build a model and perform the mass calibration
          
          for (s in 1:end(rtime(run_data))[1]){
            
            spectrum_name <- ls(run_data@assayData)[s]
            
            spectrum <- data.frame(index = 1:length(run_data@assayData[[spectrum_name]]@mz),
                                   mz = run_data@assayData[[spectrum_name]]@mz,
                                   intensity = run_data@assayData[[spectrum_name]]@intensity)
            
            #Lower lock mass
            
            lock_mass <- parameters_df$ref.mass.one[1]
            
            cal_para_1 <- calibration_parameters(spectrum = spectrum,
                                                 lock_mass = lock_mass,
                                                 minimum_counts = minimum_counts,
                                                 mass_window = mass_window)
            
            #Upper lock mass
            
            lock_mass <- parameters_df$ref.mass.two[1]
            
            cal_para_2 <- calibration_parameters(spectrum = spectrum,
                                                 lock_mass = lock_mass,
                                                 minimum_counts = minimum_counts,
                                                 mass_window = mass_window)
            
            if (is.null(cal_para_1) | is.null(cal_para_2)){
              next
            }
            
            ######3.6.1 Develop correction model######
            
            model_data <- data.frame("x" = c(cal_para_1[1], cal_para_2[1]),
                                     "y" = c(cal_para_1[2], cal_para_2[2]))
            
            model <- lm(y ~ x, model_data)
            
            ######3.6.2 Apply correction model######
            
            correction_vector <- c(model[["coefficients"]][["x"]] * spectrum$mz) + model[["coefficients"]][["(Intercept)"]]
            
            correction_vector[1:cal_para_1[3]] <- cal_para_1[2]
            
            correction_vector[cal_para_2[3]:length(correction_vector)] <- cal_para_2[2]
            
            run_data@assayData[[spectrum_name]]@mz <- run_data@assayData[[spectrum_name]]@mz + correction_vector
            
          }
          
          run_data
          
        })
      }
      
      
      print("Mass Calibration Complete")
      
      timing$stage("Prepare extraction windows and spectra", data_file_names[d])
     #####3.7 Extract Electropherograms#####
      
      #Update progress bar
      incProgress(1/total_steps, detail = paste("Extracting Electropherograms"))
      
      print("Extracting Electropherograms")
      
      #Define mass error in ppm
      
      mass_error_vec <- c(is_df$extraction.window.ppm, mass_df$extraction.window.ppm)
      
      #Create a matrix of minimum and maximum m/z values for each internal standard and metabolite
      mz_vec <- c(is_df$mz, mass_df$mz)
      min <- mz_vec - mz_vec * mass_error_vec / 1000000
      max <- mz_vec + mz_vec * mass_error_vec / 1000000
      mzr <- matrix(c(min, max), ncol = 2)
      
      #Process/check for correcting any unsorted MZ values in extracted spectra in the MSnExp obj. No clue why this occurs but it can be corrected.
      #Extract the spectra from the assayData environment stored in the run_data MSnExp
      assay_data_env <- assayData(run_data)
      
      #List all the objects in the assayData environment
      spectra_names <- ls(assay_data_env)
      
      #Create synthetic environment to merge with the real assayData environment component in run_data after spectra have been sorted
      sorted_assay_data_env <- new.env()
      
      #Loop through each spectrum in run_data and check if m/z values are unsorted AND if m/z values are negative. 
      for (i in spectra_names) {
        spectrum <- get(i, assay_data_env) 
        
        mz_values <- mz(spectrum)
        intensity_values <- intensity(spectrum)
        
        #Find any negative m/z values and replace them.
        valid_idx <- mz_values >= 0
        mz_values <- mz_values[valid_idx]
        intensity_values <- intensity_values[valid_idx]
        
        if (length(mz_values) == 0) {
          cat("Spectrum", i, "has no valid MZ values after removing negatives\n")
          next
        }
        
        
        #Check if m/z values are unsorted
        if (is.unsorted(mz_values)) {
          cat("Spectrum", i, "has unsorted MZ values\n")
          
          #Sort the problematic spectrum
          
          order_idx <- order(mz_values)
          
          #Create a sorted Spectrum1 object which will then be merged into the original run_data assayData. This allows for correction of the unsorted spectra.
          sorted_spectrum <- new("Spectrum1",
                                 mz = mz_values[order_idx],
                                 intensity = intensity_values[order_idx],
                                 rt = rtime(spectrum),
                                 tic = sum(intensity_values[order_idx]), 
                                 centroided = centroided(spectrum),
                                 fromFile = fromFile(spectrum))
          
          #Confirm spectra is sorted
          cat("Sorted Spectrum", i, "\n")  
        } else {
          #If already sorted, keep the original spectrum
          sorted_spectrum <- spectrum  
        }
        
        #Store the sorted spectrum in the new environment
        assign(i, sorted_spectrum, sorted_assay_data_env)
      }
     
     #Directly update the assayData slot of run_data
     slot(run_data, "assayData") <- sorted_assay_data_env
      
      timing$stage("Extract electropherograms", data_file_names[d])
      #Extract electropherograms
      electropherograms <- chromatogram(run_data,
                                        mz = mzr,
                                        rt = c(0,end(rtime(run_data))),
                                        aggregationFun = "mean",
                                        missing = 0,
                                        msLevel = 1)
      
      timing$stage("Build eie_df", data_file_names[d])
      eie_df <- data.frame("mt.seconds" = electropherograms[1]@rtime)
      
      for (n in 1:length(name_vec)){
        temp_df <- data.frame(electropherograms[n]@intensity)
        colnames(temp_df) <- paste(name_vec[n], "intensity", sep = " ")
        eie_df <- cbind(eie_df,temp_df)
      }
      
      print("Extraction Complete")
      
      timing$stage("Smooth electropherograms (if enabled)", data_file_names[d])
      #####3.8 Smooth Intensity Vectors#####
      
      #Check if mass calibration should be applied
      
      smoothing_response <- parameters_df$apply.smoothing
      
      #Confirm response is "Yes" or "No". Otherwise, produce an error.
      
      if(smoothing_response != "Yes" & smoothing_response != "No"){
        stop(paste("apply.smoothing parameter should be 'Yes' or 'No'. ", "'", smoothing_response, "' is not an acceptable input.", sep = ""))
      }
      
      if (smoothing_response == "No"){
        
        print("Skipping Mass Calibration")
        
      }
      
      if (smoothing_response == "Yes"){
        
        #Update progress bar
        incProgress(1/total_steps, detail = paste("Smoothing Electropherograms"))
        
        print("Smoothing Electropherograms")
        
        
        smoothing_kernel_vec <- c(is_df$smoothing.kernel, mass_df$smoothing.kernel)
        smoothing_strength_vec <- c(is_df$smoothing.strength, mass_df$smoothing.strength)
        
        for (n in 1:length(name_vec)){ 
          Smooth <- with(eie_df, 
                         ksmooth(x = mt.seconds, 
                                 y = eie_df[,n + 1], 
                                 kernel = smoothing_kernel_vec[n], 
                                 bandwidth = smoothing_strength_vec[n]))
          eie_df[,n + 1] <- Smooth[["y"]]
        }
        
        #Clean-up global environment
        
        rm(list = c("electropherograms", "mzr", "Smooth", "temp_df", "max", "min", 
                    "n", "mass_error_vec", "run_data", "smoothing_strength_vec","smoothing_kernel_vec","assay_data_env","spectra_names","sorted_assay_data_env"))
      }
      
      renderText("Electropherograms Smoothing Complete")
      print("Electropherograms Smoothing Complete")
      
      timing$stage("Internal standard peak detection, integration and filtering", data_file_names[d])
      #####3.9 Internal Standard Peak Detection, Integration, and Filtering#####
      
      #Update progress bar
      incProgress(1/total_steps, detail = paste("Performing Peak Picking and Filtering for Internal Standards"))
      
      print("Performing Peak Picking and Filtering for Internal Standards")
      
      ######3.9.1 Peak detection######
      
      is_peaks_df <- local({
        
        n <- 1
        
        for (s in 1:num_of_is){
          
          rle_output <- eie_df[,s+1] %>%
            diff() %>%
            sign() %>%
            rle()
          
          consecutive_runs <- which(rle_output$lengths > n & rle_output$values == 1)
          consecutive_runs <- subset(consecutive_runs, (consecutive_runs + 1) %in% (which(rle_output$lengths > n)) == TRUE)
          
          run_lengths <- cumsum(rle_output$lengths) + 1
          
          start <- eie_df$mt.seconds[run_lengths[consecutive_runs - 1]]
          apex <- eie_df$mt.seconds[run_lengths[consecutive_runs]]
          end <- eie_df$mt.seconds[run_lengths[consecutive_runs + 1]]
          
          #For FWHM calculations I will also add intensity values here as well
          
          start_intensity <- eie_df[run_lengths[consecutive_runs - 1], (s+1)]
          apex_intensity <- eie_df[run_lengths[consecutive_runs], (s+1)]
          end_intensity <- eie_df[run_lengths[consecutive_runs + 1], (s+1)]
          
          #Account for peaks that start immediately during the analysis
          
          if(length(start) != length(apex)){
            start <- append(start, 0, 0)
            start_intensity <- append(start_intensity, 0, 0)
          }
          
          #Create a data frame containing the start, apex, and end migration times of each peak
          
          peak_df <- data.frame(start,
                                apex,
                                end,
                                start_intensity,
                                apex_intensity,
                                end_intensity)
          
          ######3.9.2 Migration time filtering######
          
          #Filter peaks that are outside migration time limits
          
          peak_df <- subset(peak_df, peak_df$apex >= is_df$min.mt.min[s] * 60 & peak_df$apex <= is_df$max.mt.min[s] * 60)
          
          ######3.9.3 Integrate peaks######
          
          peak_area_vector = c(1:nrow(peak_df))
          
          #If the length of peak_area_vector is less than the number of injections, 
          #print an error and suggest a solution
          
          if (length(peak_area_vector) < num_of_injections){
            cat(paste("Warning: ", name_vec[s], " EIE contains insufficient data. If an error occurs try:
        1. Increasing the extraction.window.ppm parameter 
        2. Adjusting min.mt.min and max.mt.min parameters", sep =""))
          }
          
          for (p in 1:nrow(peak_df)){
            
            peak_area_vector[p] <- AUC(eie_df$mt.seconds,
                                       eie_df[,s+1],
                                       method = "trapezoid",
                                       from = peak_df[p,1],
                                       to = peak_df[p,3],
                                       absolutearea = FALSE,
                                       na.rm = FALSE)
            
            #Perform baseline correction
            
            peak_area_vector[p] <- peak_area_vector[p] - (peak_df[p,3] - peak_df[p,1]) * min(peak_df[p,4], peak_df[p,6])
          }
          
          peak_df <- cbind(peak_df, peak_area_vector)
          
          #rename peak_df columns
          
          colnames(peak_df) <- c(paste(name_vec[s], "start.seconds", sep = "."),
                                 paste(name_vec[s], "apex.seconds", sep = "."),
                                 paste(name_vec[s], "end.seconds", sep = "."),
                                 paste(name_vec[s], "start_intensity", sep = "."),
                                 paste(name_vec[s], "apex_intensity", sep = "."),
                                 paste(name_vec[s], "end_intensity", sep = "."),
                                 paste(name_vec[s], "peak.area", sep = "."))
          
          #Retain peak_df for future filtering steps
          
          peak_df_fill <- peak_df
          
          ######3.9.4 FWHM filtering######
          
          #Do not not apply the FWHM filter if the number of injections is equal to 1
          
          if (parameters_df$number.of.injections != 1) {
            
            #Find the peak intensity at half the peak height
            
            intensity_fwhm <- peak_df[,4] + (peak_df[,5] - peak_df[,4])/2 
            
            #Find the migration times closest to these intensities within each peak
            
            fwhm_vec <- vector()
            single_eie <- eie_df[,c(1,s+1)]
            
            for(p in 1:nrow(peak_df)){
              
              single_eie_temp <- subset(single_eie, single_eie$mt.seconds >= peak_df[p,1] & single_eie$mt.seconds <= peak_df[p,2])
              fwhm_mt_left <- single_eie_temp$mt.seconds[which.min(abs(single_eie_temp[,2] - intensity_fwhm[p]))]
              
              single_eie_temp <- subset(single_eie, single_eie$mt.seconds <= peak_df[p,3] & single_eie$mt.seconds >= peak_df[p,2])
              fwhm_mt_right <- single_eie_temp$mt.seconds[which.min(abs(single_eie_temp[,2] - intensity_fwhm[p]))]
              
              fwhm_vec <- append(fwhm_vec, fwhm_mt_right - fwhm_mt_left)
              
            }
            
            peak_df$fwhm <- fwhm_vec
            
            #Determine the fwhm of the peaks (n = num_of_injections) with the greatest area
            
            df_temp <- peak_df[order(-peak_df[,8]),]
            df_temp <- df_temp[1:num_of_injections,]
            
            fwhm_cutoff <- median(df_temp$fwhm)
            
            peak_df <- subset(peak_df, peak_df$fwhm <= (fwhm_cutoff * is_df$peak.fwhm.tolerance.multiplier[s]))
            peak_df <- peak_df[,c(1:7)]
            
          }
          
          #subset peak_df so that only the peaks (n = number.of.injections) with the greatest area are kept
          
          cut_off <- sort(peak_df[,7], decreasing = TRUE)[num_of_injections]
          
          peak_df <- subset(peak_df, peak_df[,7] >= cut_off)
          
          ######3.9.5 Peak space filtering######
          
          #Do not apply peak space filtering if the number of injections is equal to 1
          
          if(num_of_injections != 1){
            
            #Determine the upper and lower migration time limits for space between peaks
            
            median_space <- peak_df[,2] %>%
              diff() %>%
              median()
            
            median_space_tol <- is_df$peak.space.tolerance.percent[s] / 100
            
            median_space_lower_lim <- median_space - median_space * median_space_tol
            median_space_upper_lim <- median_space + median_space * median_space_tol
            
            #Check if peaks migrate within the tolerance limits
            
            peak_space_tol_check <- between(diff(peak_df[,2]), median_space_lower_lim, median_space_upper_lim)
            
            if(all(peak_space_tol_check) != TRUE){
              bad_space <- which(peak_space_tol_check == FALSE)
            }else{
              bad_space <- NA
            }
            
            ######3.9.6 Scenario 1######
            #Only one bad space is detected
            
            if(length(bad_space) == 1 & is.na(bad_space[1]) == FALSE){
              
              false_peak_diff <- peak_df[num_of_injections, 2] - peak_df[(num_of_injections - 1), 2]
              
              #This algorithm always assumes the final peak is false - likely due to carryover
              #It is possible the first peak is false but this seems less likely
              #Final peak only removed if case 1 does not produce a duplicate
              
              #Case 1- An interior peak is missing (usually a blank)
              
              if(bad_space != (num_of_injections - 1)){
                
                #Since the we know the bad space is not at the end, use the space after the bad space to find the expected apex
                
                expected_peak_apex <- peak_df[bad_space[1],2] + peak_df[(bad_space[1] + 2),2] - peak_df[(bad_space[1] + 1),2]
                
              }
              
              #Case 2 -  Final space is false and less than median - suspect that true final peak was missed
              
              if(bad_space == (num_of_injections - 1) & false_peak_diff < median_space_upper_lim){
                expected_peak_apex <- peak_df[(num_of_injections - 1 ),2] + median_space
              }
              
              #Case 3 - Final space is false and greater than median - suspect that peak 1 was missed
              
              if(bad_space == (num_of_injections - 1) & false_peak_diff > median_space_lower_lim){
                expected_peak_apex <- peak_df[1,2] - median_space
              }
              
              peak <- which.min(abs(peak_df_fill[,2] - expected_peak_apex))
              
              #If the nearest peak is too far from the expected migration time use a place holder
              
              if (abs(peak_df_fill[peak,2] - expected_peak_apex) > (median_space / 2)){
                
                nearest_mt <- (eie_df$mt.seconds - expected_peak_apex) %>%
                  abs() %>%
                  which.min(.)
                
                peaks <- data.frame(eie_df[nearest_mt,1],
                                    eie_df[nearest_mt,1],
                                    eie_df[nearest_mt,1],
                                    eie_df[nearest_mt,s + 1],
                                    eie_df[nearest_mt,s + 1],
                                    eie_df[nearest_mt,s + 1],
                                    0)
                
                colnames(peaks) <- colnames(peak_df)
                
                peak_df <- rbind(peak_df, peaks)
                
                peak_df <- peak_df[order(peak_df[,2]),]
                
              }else{
                
                peak_df <- rbind(peak_df, peak_df_fill[peak,])
                
                peak_df <- peak_df[order(peak_df[,2]),]
                
              }
              
              #Remove bad peak
              
              if (any(duplicated(peak_df[,2]))){
                
                peak_df <- peak_df[-c(which(duplicated(peak_df[,2]))),]
                
              }else{
                
                peak_df <- peak_df[1:(num_of_injections),]
                
              }
            }
            
            ######3.9.7 Scenario 2######
            #Two or three bad spaces are detected 
            
            if(length(bad_space) == 2 | length(bad_space) == 3){
              
              peaks_to_check <- c(bad_space, (bad_space + 1)) %>%
                unique() %>%
                sort()
              
              #Outside peaks must be checked last
              
              peaks_to_check <- c(peaks_to_check[-1], peaks_to_check[1])
              
              #Loop through each suspect bad peak and remove them iteratively until the number of bad spaces reaches 1
              
              for (p in 1:length(peaks_to_check)){
                
                num_bad_space <- peak_df[,2] %>%
                  .[-c(peaks_to_check[p])] %>%
                  diff(.) %>%
                  between (median_space_lower_lim, median_space_upper_lim) 
                num_bad_space <- length(which(num_bad_space == FALSE))
                
                if (num_bad_space == 1){
                  peak_df <- peak_df[-c(peaks_to_check[p]),]
                  break
                }
              }
              
              #To avoid errors where removing a peak results in 0 bad spaces only fill
              #in gap if nrow(peak_df) == number of injections - 1
              
              if(nrow(peak_df) == (num_of_injections - 1)){
                
                #Find the peak gap and calculate the expected migration time for the missing peak
                
                gap <- which(between(diff(peak_df[,2]), median_space_lower_lim, median_space_upper_lim) == FALSE)
                
                expected_peak_apex <- (peak_df[(gap[1] + 1),2] - peak_df[gap[1],2])/2 + peak_df[gap[1],2]
                
                #Find the nearest peak in the peak_df_fill data frame to the expected migration time
                #Avoid duplicate peaks by not using exisitng peaks in peak_df
                
                peak_df_fill <- subset(peak_df_fill, !(peak_df_fill[,2] %in% peak_df[,2]))
                
                peak <- which.min(abs(peak_df_fill[,2] - expected_peak_apex))
                peak_df <- rbind(peak_df, peak_df_fill[peak,])
                
                peak_df <- peak_df[order(peak_df[,2]),]
                
              }
            }
          }
          
          #Summarize peak_df data in is.peak_df
          
          if(s == 1){
            is_peaks_df = peak_df
          }else{
            is_peaks_df = cbind(is_peaks_df, peak_df)  
          }
        }
        
        is_peaks_df
        
      })
      
      #Make a data frame containing the apex migration times of the internal standards
      #to be used to filter metabolite peaks
      
      is_mt_df <- is_peaks_df[,seq(from = 2,
                                   to = ncol(is_peaks_df),
                                   by = 7)]
      
      print("Peak Picking and Filtering for Internal Standards Complete")
      
      timing$stage("Metabolite peak detection, integration and filtering", data_file_names[d])
      #####3.10 Metabolite  Peak Detection, Integration, and Filtering#####
      
      #Update progress bar
      incProgress(1/total_steps, detail = paste("Performing Peak Picking and Filtering for Analytes"))
      
      print("Performing Peak Picking and Filtering for Analytes")
      
      ######3.10.1 Peak detection######
      
      metabolite_peaks_df <- local ({
        
        for (m in (num_of_is + 1):length(name_vec)){
          
          peak_df <- data.frame()
          
          #Determine the start, apex, and end of peaks. Use the user defined value "n" to detect peaks.
          #If n results in fewer peaks then injection, decrease n by 1 and repeat
          
          n <- 1
          
          rle_output <- eie_df[,m + 1] %>%
            diff() %>%
            sign() %>%
            rle()
          
          consecutive_runs <- which(rle_output$lengths > n & rle_output$values == 1)
          consecutive_runs <- subset(consecutive_runs, (consecutive_runs + 1) %in% (which(rle_output$lengths > n)) == TRUE)
          
          run_lengths <- cumsum(rle_output$lengths) + 1
          
          start <- eie_df$mt.seconds[run_lengths[consecutive_runs - 1]]
          apex <- eie_df$mt.seconds[run_lengths[consecutive_runs]]
          end <- eie_df$mt.seconds[run_lengths[consecutive_runs + 1]]
          
          #I will also add intensity values here as well
          
          start_intensity <- eie_df[run_lengths[consecutive_runs - 1], (m+1)]
          apex_intensity <- eie_df[run_lengths[consecutive_runs], (m+1)]
          end_intensity <- eie_df[run_lengths[consecutive_runs + 1], (m+1)]
          
          #Account for peaks that start immediately during the analysis
          
          if(length(start) != length(apex)){
            start <- append(start, 0, 0)
            start_intensity <- append(start_intensity, 0, 0)
          }
          
          #Create a data frame containing the start, apex, and end migration times of each 
          #peak in addition to required intensities for FWHM calculations
          
          peak_df <- data.frame(start,
                                apex,
                                end,
                                start_intensity,
                                apex_intensity,
                                end_intensity)
          
          ######3.10.2 Filter peaks by peak width######
          
          #Define a minimum peak width cut off in seconds. Remove peaks with a width <= cutoff
          
          min_width_cut_off <- mass_df$minimim.peak.width.seconds[m - num_of_is]
          
          peak_df <- subset(peak_df, (peak_df$end - peak_df$start) >= min_width_cut_off)
          
          ######3.10.3 Integrate peaks######
          
          peak_area_vector = c(1:nrow(peak_df))
          
          #If the length of peak_area_vector is less than the number of injections, 
          #print an error and suggest a solution
          
          if (length(peak_area_vector) < num_of_injections){
            cat(paste("Warning: ", name_vec[m], " EIE contains insufficient data. If an error occurs try:
        1. Increasing the extraction.window.ppm parameter 
        2. Decreasing the minimim.peak.width.seconds parameter", sep =""))
          }
          
          for (p in 1:nrow(peak_df)){
            
            peak_area_vector[p] <- AUC(eie_df$mt.seconds,
                                       eie_df[,m+1],
                                       method = "trapezoid",
                                       from = peak_df[p,1],
                                       to = peak_df[p,3],
                                       absolutearea = FALSE,
                                       na.rm = FALSE)
            
            peak_area_vector[p] <- peak_area_vector[p] - (peak_df[p,3] - peak_df[p,1]) * min(peak_df[p,4], peak_df[p,6])
          }
          
          peak_df <- cbind(peak_df, peak_area_vector)
          
          #rename peak_df columns
          
          colnames(peak_df) <- c(paste(name_vec[m], "start.seconds", sep = "."),
                                 paste(name_vec[m], "apex.seconds", sep = "."),
                                 paste(name_vec[m], "end.seconds", sep = "."),
                                 paste(name_vec[m], "start_intensity", sep = "."),
                                 paste(name_vec[m], "apex_intensity", sep = "."),
                                 paste(name_vec[m], "end_intensity", sep = "."),
                                 paste(name_vec[m], "peak.area", sep = "."))
          
          ######3.10.4 Filter peaks######
          
          #Filter peaks based on smallest mt difference 
          
          #Determine the expected migration times of the metabolites
          
          n <- m - num_of_is
          left_is <- paste(mi_df$left_is[n], ".apex.seconds", sep ="")
          right_is <- paste(mi_df$right_is[n], ".apex.seconds", sep ="")
          rmt_is <- paste(mi_df$description[n], ".apex.seconds", sep ="")
          
          if(mi_df$description[n] == "mi"){
            
            mf_vec <- mi_df[n,5:ncol(mi_df)] %>%
              unlist() %>%
              unname()
            
            expected_mt <- mf_vec * (is_mt_df[,which(colnames(is_mt_df) == right_is)] - is_mt_df[,which(colnames(is_mt_df) == left_is)]) + is_mt_df[,which(colnames(is_mt_df) == left_is)]
            
          }else{
            
            mf_vec <- mi_df[n,5:ncol(mi_df)] %>%
              unlist() %>%
              unname()
            
            expected_mt <- mf_vec * is_mt_df[,which(colnames(is_mt_df) == rmt_is)] 
            
          }
          
          #Filter peak_df for peaks within migration time tolerance
          
          migration_window <- mass_df$migration.window.seconds[m - num_of_is]
          
          for (i in 1:num_of_injections){
            
            peaks <- peak_df %>%
              filter(., peak_df[,2] <= expected_mt[i] + migration_window & 
                       peak_df[,2] >= expected_mt[i] - migration_window)
            
            #If more than one peak is found choose the nearest one
            
            if(nrow(peaks) > 1){
              peaks <- (peak_df[,2] - expected_mt[i]) %>%
                abs() %>%
                which.min(.)
              peaks <- peak_df[peaks,]
            }
            
            #If no peak is found, generate a place holder
            
            if(nrow(peaks) == 0){
              
              nearest_mt <- (eie_df$mt.seconds - expected_mt[i]) %>%
                abs() %>%
                which.min(.)
              
              peaks <- data.frame(eie_df[nearest_mt,1],
                                  eie_df[nearest_mt,1],
                                  eie_df[nearest_mt,1],
                                  eie_df[nearest_mt,m + 1],
                                  eie_df[nearest_mt,m + 1],
                                  eie_df[nearest_mt,m + 1],
                                  0)
              
              colnames(peaks) <- colnames(peak_df)
            }
            
            if(i == 1){
              filtered_peaks_df <- peaks
            }else{
              filtered_peaks_df <- rbind(filtered_peaks_df, peaks)
            }
          }
          
          ######3.10.5 Filter peaks outside of run time limits######
          
          #Set a place holder for peaks where the expected migration time > total run time
          
          total_run_time <- eie_df$mt.seconds[nrow(eie_df)]
          
          late_peaks <- (expected_mt > total_run_time) %>%
            which()
          
          #Find migration times to use as placeholders that do not belong to other identified peaks
          
          mt <- tail(eie_df$mt.seconds, n = 15)
          
          mt <- mt[!(mt %in% filtered_peaks_df[,2])]
          
          for (i in late_peaks){
            mt_temp <- mt[i]
            filtered_peaks_df[i,] <- c(mt_temp,
                                       mt_temp,
                                       mt_temp,
                                       eie_df[which(eie_df$mt.seconds == mt_temp) ,m + 1],
                                       eie_df[which(eie_df$mt.seconds == mt_temp) ,m + 1],
                                       eie_df[which(eie_df$mt.seconds == mt_temp) ,m + 1],
                                       0)
          }
          
          
          ######3.10.6 Filter duplicated peaks######
          
          #If the same peak is assigned to multiple injection numbers, reapply rmt filter with more austere rmt tolerances
          #New rmt tolerance will be the original / count, which starts at 2 and increases by 1 each iteration
          
          count = 2
          
          while (any(duplicated(filtered_peaks_df[,2])) & count < 100){
            
            strict_migration_window <- migration_window/count
            
            #find rows with duplicated values 
            
            duplicate_location <- filtered_peaks_df[,2] %>%
              duplicated() %>%
              which()
            
            duplicate_rows <- which(filtered_peaks_df[,2] %in% filtered_peaks_df[duplicate_location,2])
            
            #reapply filtering for these peaks with the more strict rmt tolerance
            
            for (r in duplicate_rows){
              
              peaks <- peak_df %>%
                filter(., peak_df[,2] <= expected_mt[r] + strict_migration_window & 
                         peak_df[,2] >= expected_mt[r] - strict_migration_window)
              
              #If more than one peak is found choose the nearest one
              
              if(nrow(peaks) > 1){
                peaks <- (peak_df[,2] - expected_mt[r]) %>%
                  abs() %>%
                  which.min(.)
                peaks <- peak_df[peaks,]
              }
              
              #If no peak is found, generate a place holder
              
              if(nrow(peaks) == 0){
                
                nearest_mt <- (eie_df$mt.seconds - expected_mt[r]) %>%
                  abs() %>%
                  which.min(.)
                
                peaks <- data.frame(eie_df[nearest_mt,1],
                                    eie_df[nearest_mt,1],
                                    eie_df[nearest_mt,1],
                                    eie_df[nearest_mt,m + 1],
                                    eie_df[nearest_mt,m + 1],
                                    eie_df[nearest_mt,m + 1],
                                    0)
                
                colnames(peaks) <- colnames(peak_df)
              }
              
              filtered_peaks_df[r,] <- peaks
              
            }
            
            count = count + 1
          }
          
          #If the duplicate peak filter fails (count = 100) then generate a place holder peak_df
          #and generate a warning that the filter failed. This filter fails when two or more expected migration 
          #times are too close to each other.
          
          if (count == 100){
            
            mt.temp <- eie_df$mt.seconds[seq(1, num_of_injections * 10, 10)]
            
            filtered_peaks_df[,c(1:3)] <- mt.temp
            filtered_peaks_df[,c(4:7)] <- 0
            
          }
          
          ######3.10.7 Filter using peak spaces######
          
          #Do not apply peak space filtering if the number of injections is equal to 1
          
          if(num_of_injections != 1){
            
            #Get peak space tolerance
            
            space_tol <- mass_df$peak.space.tolerance.percent[m - num_of_is] / 100
            
            #Make a vector containing all the expected space lengths
            
            space_vec <- expected_mt %>%
              diff()
            
            #Define upper and lower peak space limits
            
            space_lower_lim <- space_vec - space_vec * space_tol
            space_upper_lim <- space_vec + space_vec * space_tol
            
            #Check if peaks migrate within the tolerance limits
            
            peak_space_tol_check <- between(diff(filtered_peaks_df[,2]), space_lower_lim, space_upper_lim)
            
            if(all(peak_space_tol_check) != TRUE){
              bad_space <- which(peak_space_tol_check == FALSE)
            }else{
              bad_space <- NA
            }
            
            #identify which peaks are potentially incorrectly assigned (bad peaks)
            #these are peaks before and after each bad space
            
            bad_peaks <- c(bad_space, bad_space + 1) %>%
              unique() %>%
              sort()
            
            #Define a count that will be used to modify the peak space tolerance
            
            count = 4
            
            #Keep unaltered filtered peaks data frame for the event that the peak space algorithm fails
            
            filtered_peaks_df_retain <- filtered_peaks_df
            
            #set count limit to determine when the algorithm fails
            
            count_limit = 100
            
            #Identify bad peaks, and replace them with peaks meeting peak space criteria
            #If the number of bad peaks is equal to the number of injections, do not apply this filter
            
            while(length(bad_peaks) > 0 & length(bad_peaks) < (num_of_injections - 1) & count < count_limit){
              
              #define remaining peaks which are correctly assigned (good peaks)
              
              good_peaks <- c(1:num_of_injections) %>%
                setdiff(., c(bad_peaks))
              
              #Use the median of the expected peak space times to find peaks
              
              peak_tolerance <- expected_mt %>%
                diff() %>%
                median () / count
              
              for (b in 1:length(bad_peaks)){
                
                #find the nearest good peak neighbor for each bad peak
                
                nearest_good_peak <- (good_peaks - bad_peaks[b]) %>%
                  abs() %>%
                  which.min()
                
                #calculate the expected migration time 
                
                expected_mt <- filtered_peaks_df[good_peaks[nearest_good_peak], 2] - 
                  (good_peaks[nearest_good_peak] - bad_peaks[b]) * median(space_vec)
                
                #find peaks nearest to the expected migration time within the tolerance
                
                peaks <- peak_df %>%
                  filter(., peak_df[,2] <= expected_mt + peak_tolerance & peak_df[,2] >= expected_mt - peak_tolerance)
                
                #if more than one peak is found, select the closest one
                
                if(nrow(peaks) > 1){
                  peaks <- (peak_df[,2] - expected_mt) %>%
                    abs() %>%
                    which.min(.)
                  peaks <- peak_df[peaks,]
                }
                
                #if no peaks are found, define a place holder
                
                if(nrow(peaks) == 0){
                  
                  nearest_mt <- (eie_df$mt.seconds - expected_mt) %>%
                    abs() %>%
                    which.min(.)
                  
                  peaks <- data.frame(eie_df[nearest_mt,1],
                                      eie_df[nearest_mt,1],
                                      eie_df[nearest_mt,1],
                                      eie_df[nearest_mt,m + 1],
                                      eie_df[nearest_mt,m + 1],
                                      eie_df[nearest_mt,m + 1],
                                      0)
                  
                  colnames(peaks) <- colnames(peak_df)
                  
                  filtered_peaks_df[bad_peaks[b],] <- peaks
                }
                
                #if only one peak is found
                
                filtered_peaks_df[bad_peaks[b],] <- peaks
                
              }
              
              count = count + 1
              
              duplicate_location <- filtered_peaks_df[,2] %>%
                duplicated() %>%
                which()
              
              #bad peaks correspond to any rows that are not unique
              
              bad_peaks <- which(filtered_peaks_df[,2] %in% filtered_peaks_df[duplicate_location,2])
              
            }
            
            #if algorithm failed, revert back to filtered_peaks_df
            
            if (count == count_limit){
              
              filtered_peaks_df <- filtered_peaks_df_retain
              
            }
          }
          
          #Summarize filtered.peak_df data in metabolite_peak_df
          
          if(m == (num_of_is + 1)){
            metabolite_peaks_df = filtered_peaks_df
          }else{
            metabolite_peaks_df = cbind(metabolite_peaks_df, filtered_peaks_df)  
          }
          
        }
        
        metabolite_peaks_df
        
      })
      
      timing$stage("Filter peaks below LOD", data_file_names[d])
      ######3.10.8 Filter peaks below LOD######
      
      #Build a data frame to store comments for each metabolite peak
      
      comment_df <- matrix(nrow = num_of_injections, ncol = num_of_metabolites, "") %>%
        as.data.frame
      
      colnames(comment_df) <- mass_df$name
      
      #Loop through each metabolite and see if its area is below the LOD threshold
      
      for (m in 1:num_of_metabolites){
        
        #Determine the noise of the electropherogram
        #Fine the noise levels in 60 seconds intervals
        
        region_start <- seq(1, nrow(eie_df), 60)
        region_end <- seq(60, nrow(eie_df), 60)
        length(region_start) <- length(region_end)
        
        #Generate a vector to store noise data
        
        noise_vec <- rep(NA, length(region_start))
        
        #Define a function to calculate noise
        
        noise_calculation <- function(temp_noise) {
          mean(temp_noise) + mass_df$snr.threshold[m] * sd(temp_noise)
        }
        
        #Calculate the noise in each region
        
        for (r in 1:length(region_start)){
          temp_noise <- eie_df[region_start[r]:region_end[r], m + num_of_is + 1]
          noise_vec[r] <- noise_calculation(temp_noise)
        }
        
        #Define the noise as the 20th percentile noise region
        
        noise <- noise_vec %>%
          sort()
        
        noise <- noise[as.integer(length(noise)/5)]
        
        peak_area_df <- metabolite_peaks_df[,seq(7, ncol(metabolite_peaks_df), 7)]
        
        comment_df[,m] <- ifelse(peak_area_df[,m] < noise, "<LOD", comment_df[,m])
        
        ###Annotate injections that are not detected
        
        comment_df[,m] <- ifelse(peak_area_df[,m] == 0, "NPD", comment_df[,m])
        
      }
      
      timing$stage("Filter interfered peaks", data_file_names[d])
      ######3.10.9 Filter interfered peaks######
      
      #Build an interference data frame since some are metabolites and some are internal standards
      
      interference_df <- cbind(is_peaks_df[,seq(2, ncol(is_peaks_df), 7)],
                               metabolite_peaks_df[,seq(2, ncol(metabolite_peaks_df), 7)])
      
      for (m in 1:num_of_metabolites){
        
        #Skip metabolites with no reported interference
        
        if(is.na(mass_df$interference[m])){
          next
        }
        
        #Get the names of the interferences from mass_df
        
        interferences <- strsplit(mass_df$interference[m], ", ") %>%
          unlist()
        
        for (k in 1:length(interferences)){
          
          interference <- paste(interferences[k], ".apex.seconds", sep = "")
          
          #Check if interference appears as a metabolite or internal standard.
          #If not, provide a warning
          
          if(!(interferences[k] %in% name_vec)){
            print(paste("Error: ", "Interference ", interferences[k], " is not an analyte or internal standard", sep = ""))
          }
          
          #Check if a interference window was provided. If not, produce an error.
          
          if(is.na(mass_df$interference.comigration.threshold.seconds[m])){
            print(paste("Error: No interference.comigration.threshold.seconds provided for ", mass_df$name[m], sep = ""))
          }
          
          #See if there is any overlap between the metabolite peak and its interference 
          
          metabolite_name <- paste(mass_df$name[m], ".apex.seconds", sep = "")
          
          for (i in 1:num_of_injections){
            for (j in 1:num_of_injections){
              
              diff_temp <- (metabolite_peaks_df[i,metabolite_name] - interference_df[j,interference]) %>%
                abs() 
              
              comment_df[i,mass_df$name[m]] <- ifelse(diff_temp < mass_df$interference.comigration.threshold.seconds[m], "Interfered", comment_df[i,mass_df$name[m]])
            }
          }
        }
      }
       
      #Combine internal standard and metabolite data frames for plotting
      
      peaks_df <- cbind(is_peaks_df, metabolite_peaks_df)
      
      #update comment data frame account for internal standards
      
      is_comment_df <- matrix(nrow = num_of_injections, ncol = nrow(is_df), "") %>%
        as.data.frame()
      
      colnames(is_comment_df) <- is_df$name
      
      comment_df <- cbind(is_comment_df, comment_df)
      

      print("Peak Picking and Filtering for Analytes Complete")
      
      timing$stage("Generate and save plots", data_file_names[d])
      #####3.11 Plotting#####
      #Update progress bar
      incProgress(1/total_steps, detail = paste("Plotting & Exporting Electropherograms"))
      print("Plotting Electropherograms")
      
      #Make a list to save plots to
      
      plot_list <- vector("list", length(name_vec))
      names(plot_list) <- name_vec
      
      for (n in 1:length(name_vec)){
        
        ######3.11.1 Create annotation and peak fill data frame######
        
        peak_mt_df <- peaks_df[,seq(from = 2, to = ncol(peaks_df), by = 7)]
        
        ann_df <- data.frame("peak.number" = c(1:num_of_injections),
                             "comment" = comment_df[,n],
                             "peak.apex.seconds" = peak_mt_df[,n],
                             "peak.height.counts" = eie_df[which(eie_df$mt.seconds %in% (peak_mt_df[,n])),n+1])
        
        #Create peak fill data frame
        
        pf_df <- data.frame("peak.number" = 1,
                            "mt.seconds" = eie_df[,1],
                            "intensity" = eie_df[,n+1])
        
        mt_vec <- vector()
        start_df <- peaks_df[,seq(from = 1, to = ncol(peaks_df), by = 7)]
        end_df <- peaks_df[,seq(from = 3, to = ncol(peaks_df), by = 7)]
        
        for (i in 1:num_of_injections){
          
          #Create a migration time vector to track where peaks elute
          
          if(comment_df[i,n] == ""){
            mt_vec_temp <- eie_df$mt.seconds[between(eie_df$mt.seconds, start_df[i,n], end_df[i,n])]
            mt_vec <- append(mt_vec, mt_vec_temp)
            
            #Update peak.number in pf_df
            
            pf_df$peak.number <- ifelse(pf_df$mt.seconds >= start_df[i,n], i, pf_df$peak.number)
            pf_df$peak.number <- as.factor(pf_df$peak.number)
          }else{
            next
          }
        }
        
        pf_df$intensity <- ifelse(pf_df$mt.seconds %in% mt_vec == TRUE, pf_df$intensity , 0)
        
        ##Add baseline intensity
        
        pf_df$baseline <- 0
        
        for (i in 1:num_of_injections){
          if(comment_df[i,n] == ""){
            lower_intensity <- min(c(peaks_df[i, n * 7 - 3], peaks_df[i, n * 7 - 1]))
            pf_df$baseline <- ifelse(pf_df$mt.seconds >= start_df[i,n] & pf_df$mt.seconds <= end_df[i,n], lower_intensity, pf_df$baseline)
          }else{
            next
          }
        }
        
        #Only retain filling data required for plotting
        
        pf_df <- subset(pf_df, pf_df$intensity != 0)
        
        #Save variables to a list
        
        plot_list[[name_vec[n]]] <- list("eie_data" = eie_df[,c(1,(n + 1))], 
                                         "annotation_data" = ann_df, 
                                         "integration_data" = pf_df,
                                         "label_data" = c(name_vec[n],
                                                          mz_vec[n]),
                                         "x_axis_data" = c(start_df[1,n],
                                                           end_df[num_of_injections,n]),
                                         "y_axis_data" = c(max(ann_df$peak.height.counts), 
                                                           max(eie_df[,(n + 1)])))
        
      }

      
      ######3.11.2 Create plot function for ggplots######
      plot_function <- function(eie_data, annotation_data, integration_data, label_data, x_axis_data, y_axis_data){
        
        extra_space <- ifelse(x_axis_data[1] > 70, 1, 0)
        
        ggplot(data = eie_data) +
          geom_line(aes(x = mt.seconds / 60, y = eie_data[,2]), colour = "grey50") +
          theme_classic() +
          coord_cartesian(xlim = c(x_axis_data[1] / 60 - extra_space, x_axis_data[2] / 60 + extra_space),
                          ylim = c(0, 1.5 * y_axis_data[1])) +
          scale_y_continuous(name = "Ion Counts",
                             labels = function(x) format(x, scientific = TRUE),
                             expand = c(0,0),
                             breaks = scales::pretty_breaks(n = 10)) +
          scale_x_continuous(name = "Migration Time (Minutes)",
                             breaks = scales::pretty_breaks(n = 10))+
          ggtitle(paste(label_data[1], " EIE", " (m/z = ", label_data[2],")",sep = ""),
                  subtitle = paste("Data File: ", data_files[d])) +
          geom_ribbon(data = integration_data,
                      aes(x = mt.seconds/60, ymax = intensity, ymin = baseline, fill = peak.number),
                      alpha =0.4) +
          geom_text(data = annotation_data,
                    label = annotation_data$peak.number,
                    size  = font_size_1,
                    family = "sans",
                    aes(x = peak.apex.seconds/60,
                        y = peak.height.counts + 0.1 * y_axis_data[1])) +
          geom_text(data = annotation_data,
                    label = annotation_data$comment,
                    size  = font_size_1,
                    family = "sans",
                    aes(x = peak.apex.seconds/60,
                        y = peak.height.counts + 0.2 * y_axis_data[1])) +
          theme(legend.position = "none",
                text = element_text(size = font_size_2, family = "sans"))
        
      }
      
      ######3.11.3 Create Plot Function for making plotly plots######
      plot_function_plotly <- function(eie_data, annotation_data, integration_data, label_data, x_axis_data, y_axis_data) {
        plot_ly() %>%
          add_lines(
            x = ~eie_data[,1] / 60,
            y = ~eie_data[,2],
            name = 'Electropherogram',
            line = list(color = 'grey')
          ) %>%
          add_ribbons(
            x = ~integration_data$mt.seconds / 60,
            ymin = ~integration_data$baseline,
            ymax = ~integration_data$intensity,
            fillcolor = 'rgba(100,100,255,0.4)',
            line = list(color = 'rgba(100,100,255,0.4)'),
            name = 'Peak Integration'
          ) %>%
          add_text(
            x = ~annotation_data$peak.apex.seconds / 60,
            y = ~annotation_data$peak.height.counts + (0.1 * y_axis_data[1]),
            text = ~annotation_data$peak.number,
            showlegend = FALSE
          ) %>%
          add_text(
            x = ~annotation_data$peak.apex.seconds / 60,
            y = ~annotation_data$peak.height.counts + (0.2 * y_axis_data[1]),
            text = ~annotation_data$comment,
            showlegend = FALSE
          ) %>%
          layout(
            title = paste(label_data[1], " EIE", " (m/z =", label_data[2], ")"),
            xaxis = list(title = 'Migration Time (Minutes)'),
            yaxis = list(title = 'Ion Counts'),
            template = 'plotly_white',
            shapes = list(
              list(
                type = "line",
                x0 = 2,  #Initial position of the red line
                y0 = 0,
                x1 = 2,
                y1 = 1,
                xref = "x",
                yref = "paper",
                line = list(color = "Red", width = 2)
              ),
              list(
                type = "line",
                x0 = 4,  #Initial position of the blue line
                y0 = 0,
                x1 = 4,
                y1 = 1,
                xref = "x",
                yref = "paper",
                line = list(color = "Blue", width = 2)
              )
            )
          ) %>% config(editable = TRUE)
      }
      
      ######3.11.4 Initialize Saving Locations and Folders######
      
      if (parameters_df$plot.format == "Sample"){
        
        data_files_name <- list.files(path = "mzML Files")
        data_files_name <- gsub(".mz5", "", data_file_names, fixed = TRUE)
        
        #Create sub-folders
        
        dir.create(file.path(file_name, "Plots", "Analytes"), showWarnings = FALSE)
        dir.create(file.path(file_name, "Plots", "Internal Standards"), showWarnings = FALSE)
        
        dir.create(path = paste(file_name, "/Plots/", "Internal Standards", "/", data_files_name[d], sep = ""),
                   showWarnings = TRUE)
        
        dir.create(path = paste(file_name, "/Plots/", "Analytes", "/", data_files_name[d], sep = ""),
                   showWarnings = TRUE)
        
        ######3.11.5 Save plots as ggplots######
        #Save Internal Standard Plots
        
        plotly_plots <- list()
        
        for (n in 1:num_of_is){
          
          folder <- "Internal Standards"
          name <- name_vec[n]
          
          font_size_1 <- 7
          font_size_2 <- 25
          
          ggsave(filename = paste(n, "_", name,".png",sep=""),
                 width = 16,
                 height = 9,
                 plot = plot_function(eie_data = plot_list[[n]][[1]], 
                                      annotation_data = plot_list[[n]][[2]], 
                                      integration_data = plot_list[[n]][[3]],
                                      label_data = plot_list[[n]][[4]],
                                      x_axis_data = plot_list[[n]][[5]],
                                      y_axis_data = plot_list[[n]][[6]]),
                 path = paste(file_name, "/Plots/", folder, "/", data_files_name[d], sep = ""))
          
        }
        
        #Save Analyte Plots
        
        for (n in (num_of_is + 1):length(name_vec)){
          
          folder <- "Analytes"
          name <- name_vec[n]
          
          #Analytes using RMT
          
          if (mi_df$description[n - num_of_is] != "mi"){
            
            font_size_1 <- 4
            font_size_2 <- 12
            
            is_index <- which(name_vec == mi_df$description[(n - num_of_is)])
            
            figure <- ggarrange(plot_function(eie_data = plot_list[[n]][[1]], 
                                              annotation_data = plot_list[[n]][[2]], 
                                              integration_data = plot_list[[n]][[3]],
                                              label_data = plot_list[[n]][[4]],
                                              x_axis_data = c(min(eie_df$mt.seconds), max(eie_df$mt.seconds)),
                                              y_axis_data = plot_list[[n]][[6]][2]), 
                                plot_function(eie_data = plot_list[[n]][[1]], 
                                              annotation_data = plot_list[[n]][[2]], 
                                              integration_data = plot_list[[n]][[3]],
                                              label_data = plot_list[[n]][[4]],
                                              x_axis_data = plot_list[[n]][[5]],
                                              y_axis_data = plot_list[[n]][[6]]), 
                                plot_function(eie_data = plot_list[[is_index]][[1]], 
                                              annotation_data = plot_list[[is_index]][[2]], 
                                              integration_data = plot_list[[is_index]][[3]],
                                              label_data = plot_list[[is_index]][[4]],
                                              x_axis_data = plot_list[[n]][[5]],
                                              y_axis_data = plot_list[[is_index]][[6]]),
                                ncol = 1, nrow = 3)
            
            ggsave(filename = paste(n, "_", name,".png",sep=""),
                   width = 16,
                   height = 9,
                   plot = figure,
                   path = paste(file_name, "/Plots/", folder, "/", data_files_name[d], sep = ""))
          }
          
          #Analytes using MI
          
          if (mi_df$description[n - num_of_is] == "mi"){
            
            font_size_1 <- 4
            font_size_2 <- 12
            
            is_index <- which(name_vec == mi_df$description[(n - num_of_is)])
            
            figure <- ggarrange(plot_function(eie_data = plot_list[[n]][[1]], 
                                              annotation_data = plot_list[[n]][[2]], 
                                              integration_data = plot_list[[n]][[3]],
                                              label_data = plot_list[[n]][[4]],
                                              x_axis_data = c(min(eie_df$mt.seconds), max(eie_df$mt.seconds)),
                                              y_axis_data = plot_list[[n]][[6]][2]), 
                                plot_function(eie_data = plot_list[[n]][[1]], 
                                              annotation_data = plot_list[[n]][[2]], 
                                              integration_data = plot_list[[n]][[3]],
                                              label_data = plot_list[[n]][[4]],
                                              x_axis_data = plot_list[[n]][[5]],
                                              y_axis_data = plot_list[[n]][[6]]),
                                ncol = 1, nrow = 2)
            
            ggsave(filename = paste(n, "_", name,".png",sep=""),
                   width = 16,
                   height = 9,
                   plot = figure,
                   path = paste(file_name, "/Plots/", folder, "/", data_files_name[d], sep = ""))
          }
        }
      }
      
      if (parameters_df$plot.format == "Metabolite"){
        
        #Save plots to their respective folders within the "Plots" folder
        data_files_name <- list.files(path = "mzML Files")
        data_files_name <- gsub(".mz5", "", data_file_names, fixed = TRUE)
        
        #Save Plots
        for (n in 1:length(name_vec)){
          
          folder <- name_vec[n]
          
          font_size_1 <- 7
          font_size_2 <- 25
          
          ggsave(filename = paste(data_files_name[d],".png",sep=""),
                 width = 16,
                 height = 9,
                 plot = plot_function(eie_data = plot_list[[n]][[1]], 
                                      annotation_data = plot_list[[n]][[2]], 
                                      integration_data = plot_list[[n]][[3]],
                                      label_data = plot_list[[n]][[4]],
                                      x_axis_data = plot_list[[n]][[5]],
                                      y_axis_data = plot_list[[n]][[6]]),
                 path = paste(file_name, "/Plots/", folder, "/", sep = ""))
          
        }
      }
      
      ######3.11.6 Save plots as editable plotly plots######
      plot_list[[name_vec[n]]] <- list("eie_data" = eie_df[,c(1,(n + 1))], 
                                       "annotation_data" = ann_df, 
                                       "integration_data" = pf_df,
                                       "label_data" = c(name_vec[n],
                                                        mz_vec[n]),
                                       "x_axis_data" = c(start_df[1,n],
                                                         end_df[num_of_injections,n]),
                                       "y_axis_data" = c(max(ann_df$peak.height.counts), 
                                                         max(eie_df[,(n + 1)])))
      
      
      #Save Internal Standard Plots to a list
      plotly_plots <- lapply(1:num_of_is, function(n) {
        
        name <- name_vec[n]
        
        plot <- plot_function_plotly(eie_data = plot_list[[n]][[1]], 
                                     annotation_data = plot_list[[n]][[2]], 
                                     integration_data = plot_list[[n]][[3]],
                                     label_data = plot_list[[n]][[4]],
                                     x_axis_data = plot_list[[n]][[5]],
                                     y_axis_data = plot_list[[n]][[6]])
        return(plot)
      })
      
      names(plotly_plots) <- sapply(1:num_of_is, function(n) {
        paste0(data_files_name[d], "_", name_vec[n])
      })

      plotly_objects <- c(plotly_objects, plotly_plots)  
      
      #Save Analyte Plots as Plotly Objects 
      for (n in (num_of_is + 1):length(name_vec)) {
        
        name <- name_vec[n]
        
        #Analytes using RMT
        if (mi_df$description[n - num_of_is] != "mi") {
          
          is_index <- which(name_vec == mi_df$description[(n - num_of_is)])
          
          plotly_figure <- plot_function_plotly(
              eie_data = plot_list[[n]][[1]],
              annotation_data = plot_list[[n]][[2]],
              integration_data = plot_list[[n]][[3]],
              label_data = plot_list[[n]][[4]],
              x_axis_data = plot_list[[n]][[5]],
              y_axis_data = plot_list[[n]][[6]])
          
          #Build plotly object before saving them
          plotly_figure <- plotly::plotly_build(plotly_figure)
          
          #Name plots based on the file being processed and the metabolite
          plot_name <- paste0(data_files_name[d], "_", name_vec[n])
          
          #Save the plotly object
          plotly_objects[[plot_name]] <- plotly_figure
        }
        
        #Saving Analytes using MI as plotly objects
        if (mi_df$description[n - num_of_is] == "mi") {
          
          plotly3 <- plot_function_plotly(
              eie_data = plot_list[[n]][[1]],
              annotation_data = plot_list[[n]][[2]],
              integration_data = plot_list[[n]][[3]],
              label_data = plot_list[[n]][[4]],
              x_axis_data = plot_list[[n]][[5]],
              y_axis_data = plot_list[[n]][[6]])
          
          #Build plotly object before saving them
          plotly3 <- plotly::plotly_build(plotly3)
          
          #Name plots based on the file being processed and the metabolite
          plot_name <- paste0(data_files_name[d], "_", name_vec[n])
          
          #Save the plotly object
          plotly_objects[[plot_name]] <- plotly3
        }
      }
      
      #Save plotly_objects to reactive varibale plotly_data
      plotly_data(plotly_objects)
      
      #Clear workspace
      rm(list = c())
      
      print("Plotting Complete")
      
      
      timing$stage("Prepare result tables", data_file_names[d])
      #####3.12 Export Data#####
      
      ######3.12.1 Generate peak area data frame######
      
      peak_area_df <- cbind("file.name" = c(data_files_name[d], rep("", num_of_injections - 1)),
                            "peak.number" = c(1:num_of_injections),
                            peaks_df[,seq(from = 7, to = ncol(peaks_df), by = 7)])
      colnames(peak_area_df)[3:(length(name_vec) + 2)] <- name_vec
      
      #Update values to include <LOD and Interfered
      
      for (i in 1:num_of_injections){
        peak_area_df[i,3:ncol(peak_area_df)] <- ifelse(comment_df[i,] == "", peak_area_df[i, 3:ncol(peak_area_df)], comment_df[i,])
      }
      
      peak_area_report <- peak_area_df
      
      
      ######3.12.2 Generate peak migration time data frame######
      
      peak_mt_df <- cbind("file.name" = c(data_files_name[d], rep("", num_of_injections - 1)),
                          "peak.number" = c(1:num_of_injections),
                          peaks_df[,seq(from = 2, to = ncol(peaks_df), by = 7)] / 60)
      colnames(peak_mt_df)[3:(length(name_vec) + 2)] <- name_vec
      
      peak_mt_report <- peak_mt_df
      
      timing$stage("Save plot_list_data.RData", data_file_names[d])
      #Save the entire plot_list to a .RData file in the subfolder as well as peaks_df
      save(plot_list, peaks_df, mi_df, peak_mt_df,peak_area_df, comment_df, file = file.path(subfolder_path, "plot_list_data.RData"))
      
      timing$stage("Write SQLite results", data_file_names[d])
      #Create SQL database if one doesnt exist and save all data to it for each file.
      db_path_pipeline <- get_db_path(file_name)
      pipeline_con <- dbConnect(SQLite(), db_path_pipeline)
      dbExecute(pipeline_con, "PRAGMA journal_mode = WAL")
      
      initialize_db(pipeline_con)
      
      for (tbl in c("annotation_data", "integration_data","plot_axis_data", "peaks_df","peak_area_df","peak_mt_df")){
        dbExecute(pipeline_con,
                  paste0("DELETE FROM ", tbl, " WHERE file_name = ?"),
                  params = list(data_files_name[d])
        )
      }
      #Write data to SQL database
      write_file_to_db(
        con          = pipeline_con,
        file_name    = data_files_name[d],
        plot_list    = plot_list,
        peaks_df     = peaks_df,
        peak_area_df = peak_area_df,
        peak_mt_df   = peak_mt_df,
        comment_df   = comment_df,
        mi_df        = mi_df
      )
      
      dbDisconnect(pipeline_con)
      print(paste("database write complete for:", data_files_name[d]))
      
      
      timing$stage("Remove temporary mz5 file", data_file_names[d])
      #Delete temporary mz5 file
      file.remove(paste(file, "temp.mz5", sep = "_"))
  
      print(paste("Completed Analysis of Data File: ", data_file_names[d], sep = ""))
      print("")
      
      timing$stage("Export CSV reports", data_file_names[d])
      ######3.12.3 Export data######
      write.csv(peak_area_report,
                file = paste(file_name, "/", "Metabolite Peak Areas.csv", sep = ""),
                row.names = FALSE)
      
      write.csv(peak_mt_report,
                file = paste(file_name, "/", "Metabolite Migration Times.csv", sep = ""),
                row.names = FALSE)
      
      timing$stage("Save plotly_objects.RData", data_file_names[d])
      #Save plotly_objects to .RDA file
      save_plotly_objects(file_name)
      timing$finish()
      })
    }#End of loop
    
    timing$finish()
    }, error = function(e) {
      timing$finish("failed", conditionMessage(e))
      stop(e)
    }, interrupt = function(e) {
      timing$finish("interrupted", "Run interrupted")
      stop(e)
    })
  })#End of main button
  
  
  ####4. Visualization tab####
  
  #####4.1 Define reactive expressions and specific required functions#####
  ######4.1.1 Define Reactive expressions######
  plot_list_data_values <- reactiveValues(
    plot_list      = NULL,
    comment_df     = NULL,
    peaks_df       = NULL,
    peak_mt_df     = NULL,
    peak_area_df   = NULL,
    mi_df          = NULL,
    file_name      = NULL
  )
  
  #Define reactive values for storing information from plot_list_data.Rdata
  run_metadata <- reactiveValues(
    comment_df    = NULL,
    plot_list     = NULL,
    peaks_df      = NULL,
    peak_mt_df    = NULL,
    peak_area_df  = NULL,
    mi_df         = NULL,
    file_name     = NULL
  )
  
  #Define reactive expression to store unaltered peak information to allow users to undo a deletion if a peak is deleted accidentally 
  previous_metadata <- reactiveValues(
    comment_df    = NULL,
    plot_list     = NULL,
    peaks_df      = NULL,
    peak_mt_df    = NULL,
    peak_area_df  = NULL,
    mi_df         = NULL,
    file_name     = NULL
  )
  
  #Reactive expression to store the path to the selected results folder. Required to grab the correct annotation_data information from the subsequent runs and useful for processing data if you have closed the app.
  main_folder <- reactive({
    req(input$results_folder, input$results_folder %in% results_folders_reactive())
    input$results_folder
  })
  
  #Reactive expression to store the path to the selected results folder
  plot_list_data <- reactive({
    req(main_folder(), input$file_selector)
    file_name <- input$file_selector
    subfolder_path <- file.path(main_folder(), "Data", file_name)
    
    #Check if the file exists before attempting to load it
    if (file.exists(file.path(subfolder_path, "plot_list_data.RData"))) {
      load(file.path(subfolder_path, "plot_list_data.RData"))
      return(list(plot_list = plot_list, peaks_df = peaks_df, peak_mt_df = peak_mt_df, peak_area_df = peak_area_df, comment_df = comment_df, mi_df = mi_df))
    } else {
      warning("File not found: ", file.path(subfolder_path, "plot_list_data.RData"))
      return(NULL)
    }
  })

  
  #Reactive expression to get the plot names associated with a specified .mz5 file from the SQL database, based on the naming pattern.
  filtered_plot_names <- reactive({
    req(db_con(), input$file_selector, main_folder())
    con       <- db_con()
    file_name <- sub("\\.mz5$", "", input$file_selector)
    tryCatch({
      metabolites <- dbGetQuery(con,
                                "SELECT DISTINCT metabolite FROM plot_axis_data
         WHERE file_name = ?
         ORDER BY metabolite",
                                params = list(file_name)
      )$metabolite
      #Return in same format as before: "<file>_<metabolite>"
      paste0(file_name, "_", metabolites)
    }, error = function(e) character(0))
  })
  
  
  #Variable to save edited plot names for rerendering without rerendering ALL plots
  modified_peak_plots <- reactiveValues(names = character())
  
  #Variable to trigger plots reloading after editing them
  plot_render_trigger <- reactiveVal(0)
  
  
  ######4.1.2 SQL Database functions######
  
  #Database path helper 
  get_db_path <- function(results_folder) {
    file.path(results_folder, "results.sqlite")
  }
  
  #Check if database exists in results folder
  db_exists <- function(results_folder){
    file.exists(get_db_path(results_folder))
  }
  
  #Function for connection to the database
  open_db <- function(results_folder){
    con <- dbConnect(RSQLite::SQLite(), get_db_path(results_folder))
    #enable write-ahead logging for better performance
    dbExecute(con, "PRAGMA journal_mode = WAL")
    return(con)
  }
  
  #Function for initializing the SQL database, creating all tables and indexes for increased efficiency.
  initialize_db <- function(con){
    
  #Create tables for storing all components of the plot_list_data RData file/ allrequired files generated from running the main peakmeister loop  
    #Table for annotation data
    dbExecute(con, "CREATE TABLE IF NOT EXISTS annotation_data (
              file_name       TEXT NOT NULL,
              metabolite      TEXT NOT NULL,
              peak_number      INTEGER NOT NULL,
              comment         TEXT,
              peak_apex_secs  REAL,
              peak_height     
              )")
    #Table for integration data
    dbExecute(con, "CREATE TABLE IF NOT EXISTS integration_data (
              file_name     TEXT NOT NULL,
              metabolite    TEXT NOT NULL,
              peak_number   INTEGER NOT NULL,
              mt_seconds    REAL NOT NULL,
              intensity     REAL NOT NULL,
              baseline      REAL NOT NULL
    )")
    #Table for storing x_axis_data, y_axis_data, and label_data components from plot_list
    dbExecute(con, "CREATE TABLE IF NOT EXISTS plot_axis_data (
              file_name     TEXT NOT NULL,
              metabolite    TEXT NOT NULL,
              x_start       REAL,
              x_end         REAL,
              y_max_annot   REAL,
              y_max_eie     REAL,
              label_name    TEXT,
              label_mz      TEXT
    )")
    #Table for storing peak_df from peakmeister enginine loop
    dbExecute(con, "CREATE TABLE IF NOT EXISTS peaks_df (
              file_name       TEXT NOT NULL,
              metabolite      TEXT NOT NULL,
              peak_number     INTEGER NOT NULL,
              start_seconds   REAL,
              apex_seconds    REAL,
              end_seconds     REAL,
              start_intensity REAL,
              apex_intensity  REAL,
              end_intensity   REAL,
              peak_area       TEXT
    )")
    #Table for storing peak_area_df information
    dbExecute(con, "CREATE TABLE IF NOT EXISTS peak_area_df (
              file_name     TEXT NOT NULL,
              peak_number   INTEGER NOT NULL,
              metabolite    TEXT NOT NULL,
              peak_area     TEXT
    )")
    #Table for storing peak_mt_df information
    dbExecute(con, "CREATE TABLE IF NOT EXISTS peak_mt_df (
              file_name       TEXT NOT NULL,
              peak_number     INTEGER NOT NULL,
              metabolite      TEXT NOT NULL,
              migration_time  REAL
    )")
    #Table for storing mi_df information
    dbExecute(con, "CREATE TABLE IF NOT EXISTS mi_df (
              metabolite    TEXT NOT NULL,
              left_is       TEXT,
              right_is      TEXT,
              description   TEXT,
              peak_indices  TEXT
    )")
    #Table for tracking edits
    dbExecute(con, "CREATE TABLE IF NOT EXISTS edit_log(
              id            INTEGER PRIMARY KEY AUTOINCREMENT,
              file_name     TEXT NOT NULL,
              metabolite    TEXT_NOT_NULL,
              edit_type     TEXT NOT NULL,
              edited_at     TEXT NOT NULL
    )")
    
    #Create indexes on (filename, metabolite) for every table for optimal data access 
    dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_ann    ON annotation_data(file_name, metabolite)")
    dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_intg   ON integration_data(file_name, metabolite)")
    dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_axis   ON plot_axis_data(file_name, metabolite)")
    dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_peaks  ON peaks_df(file_name, metabolite)")
    dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_area   ON peak_area_df(file_name, metabolite)")
    dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_mt     ON peak_mt_df(file_name, metabolite)")
  }
  
  #Function for wide -> long pivitors
  safe_col <- function(df, col){
    if (col %in% colnames(df)) df[[col]] else rep(NA_real_, nrow(df))
  }
  
  #Function to write/add each datafiles data into the SQL database, one datafile at a time as its generated in the main PeakMeister loop OR when generating a SQL file from preexisting .RData files
  write_file_to_db <- function(con, file_name, plot_list, peaks_df, peak_area_df, peak_mt_df, comment_df, mi_df) {
    metabolite_names <- names(plot_list)
    
    dbBegin(con)
    tryCatch({
      
      #annotation_data
      ann_all <- do.call(rbind,lapply(metabolite_names, function(metab) {
       ann <- plot_list[[metab]]$annotation_data
       ann_comments <- if (!is.null(comment_df) && metab %in% colnames(comment_df)){
         comment_df[ann$peak.number,metab]
       } else {
         ann$comment
       }
         data.frame(
           file_name        = file_name,
           metabolite       = metab,
           peak_number      = ann$peak.number,
           comment          = as.character(ann_comments),
           peak_apex_secs   = ann$peak.apex.seconds,
           peak_height      = ann$peak.height.counts,
           stringsAsFactors = FALSE
         )
      }))
      dbAppendTable(con, "annotation_data", ann_all)
      
      #integration_data
      intg_list <- lapply(metabolite_names, function(metab) {
        intg <- plot_list[[metab]]$integration_data
        if (is.null(intg) || nrow(intg) == 0) return(NULL)
        intg$peak.number <- as.numeric(as.character(intg$peak.number))
        data.frame(
          file_name   = file_name,
          metabolite  = metab,
          peak_number = intg$peak.number,
          mt_seconds  = intg$mt.seconds,
          intensity   = intg$intensity,
          baseline    = intg$baseline,
          stringsAsFactors = FALSE
        )
      })
      intg_all <- do.call(rbind, Filter(Negate(is.null), intg_list))
      if (!is.null(intg_all) && nrow(intg_all) > 0) {
        dbAppendTable(con, "integration_data", intg_all)
      }
      
      #plot axis data
      axis_all <- do.call(rbind, lapply(metabolite_names, function(metab){
        data.frame(
          file_name        = file_name,
          metabolite       = metab,
          x_start          = plot_list[[metab]]$x_axis_data[1],
          x_end            = plot_list[[metab]]$x_axis_data[2],
          y_max_annot      = plot_list[[metab]]$y_axis_data[1],
          y_max_eie        = plot_list[[metab]]$y_axis_data[2],
          label_name       = as.character(plot_list[[metab]]$label_data[1]),
          label_mz         = as.character(plot_list[[metab]]$label_data[2]),
          stringsAsFactors = FALSE
        )
      }))
      dbAppendTable(con, "plot_axis_data", axis_all)
      
      #Peaks_df data -> Need to pivot this data from wide format to long format
      num_peaks  <- nrow(peaks_df)
      peaks_long <- do.call(rbind, lapply(metabolite_names, function(metab) {
        data.frame(
          file_name       = file_name,
          metabolite      = metab,
          peak_number     = 1:num_peaks,
          start_seconds   = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(metab, ".start.seconds")))),
          apex_seconds    = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(metab, ".apex.seconds")))),
          end_seconds     = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(metab, ".end.seconds")))),
          start_intensity = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(metab, ".start_intensity")))),
          apex_intensity  = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(metab, ".apex_intensity")))),
          end_intensity   = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(metab, ".end_intensity")))),
          peak_area       = as.character(safe_col(peaks_df, paste0(metab, ".peak.area"))),
          stringsAsFactors = FALSE
        )
      }))
      dbAppendTable(con, "peaks_df", peaks_long)
      
      #peak_area_df datatable -> also need to change wide - long
      area_meta_cols <- colnames(peak_area_df)[
        !(colnames(peak_area_df) %in% c("file.name", "peak.number"))
      ]
      peak_area_df_chr <- as.data.frame(peak_area_df)
      peak_area_df_chr[area_meta_cols] <- lapply(
        peak_area_df_chr[area_meta_cols], as.character
      )
      area_long <- pivot_longer(
        peak_area_df_chr,
        cols = all_of(area_meta_cols),
        names_to = "metabolite", values_to = "peak_area"
      )
      area_long$file_name <- file_name
      area_long$peak_area <- as.character(area_long$peak_area)
      area_long <- area_long[, c("file_name", "peak.number", "metabolite", "peak_area")]
      colnames(area_long)[2] <- "peak_number"
      dbAppendTable(con, "peak_area_df", area_long)
      
      #peak_mt_df
      mt_meta_cols <- colnames(peak_mt_df)[
        !(colnames(peak_mt_df) %in% c("file.name","peak.number"))
      ]
      mt_long <- pivot_longer(
        as.data.frame(peak_mt_df),
        cols = all_of(mt_meta_cols),
        names_to = "metabolite", values_to = "migration_time"
      )
      mt_long$file_name <- file_name
      mt_long <- mt_long[, c("file_name","peak.number","metabolite","migration_time")]
      colnames(mt_long)[2] <- "peak_number"
      
      dbAppendTable(con, "peak_mt_df", mt_long)
      
      #mi_df
      existing_mi <- dbGetQuery(con, "SELECT COUNT(*) FROM mi_df")[[1]]
      if (existing_mi ==0){
        fixed_cols <- c("name", "left_is", "right_is","description")
        peak_cols <- colnames(mi_df)[!(colnames(mi_df) %in% fixed_cols)]
        mi_all <- do.call(rbind, lapply(1:nrow(mi_df), function(i){
          data.frame(
            metabolite = as.character(mi_df$name[i]),
            left_is    = as.character(mi_df$left_is[i]),
            right_is   = as.character(mi_df$right_is[i]),
            description = as.character(mi_df$description[i]),
            peak_indices = toJSON(
              as.list(mi_df[i, peak_cols, drop = FALSE]),
              auto_unbox = TRUE
            ), 
            stringsAsFactors = FALSE
          )
        }))
        dbAppendTable(con, "mi_df", mi_all)
      }
      
      dbCommit(con)
    }, error = function(e){
      tryCatch(dbRollback(con), error = function(e2) NULL)
      stop(paste("Database write failed for", file_name, ":", conditionMessage(e)))
    })
  }
  
  #Read one metabolites plot data from database and eie_df from .Rdata file
  read_plot_data_from_db <- function(con, results_folder, file_name, metabolite) {
    
    #eie_data: load from .RData -> in .Rdata for now since eie_df is the largest dataset for each metabolite. Its never edited so, for now its in here.
    rdata_path <- file.path(results_folder, "Data",
                            paste0(file_name, ".mz5"),
                            "plot_list_data.RData")
    #Load into new environment 
    env <- new.env()
    load(rdata_path, envir = env)
    eie_raw <- env$plot_list[[metabolite]]$eie_data
    rm(env)
    
    eie_db <- data.frame(
      mt.seconds = eie_raw[[1]],
      intensity  = eie_raw[[2]],
      stringsAsFactors = FALSE
    )
    colnames(eie_db)[2] <- paste(metabolite, "intensity")
    
    #annotation_data
    ann_db <- dbGetQuery(con,
                         "SELECT peak_number, comment, peak_apex_secs, peak_height
     FROM annotation_data
     WHERE file_name = ? AND metabolite = ?
     ORDER BY peak_number",
                         params = list(file_name, metabolite)
    )
    colnames(ann_db) <- c("peak.number", "comment",
                          "peak.apex.seconds", "peak.height.counts")
    
    #integration_data
    intg_db <- dbGetQuery(con,
                          "SELECT peak_number, mt_seconds, intensity, baseline
     FROM integration_data
     WHERE file_name = ? AND metabolite = ?
     ORDER BY peak_number, mt_seconds",
                          params = list(file_name, metabolite)
    )
    colnames(intg_db) <- c("peak.number", "mt.seconds", "intensity", "baseline")
    
    #axis + label data
    ax_db <- dbGetQuery(con,
                        "SELECT x_start, x_end, y_max_annot, y_max_eie, label_name, label_mz
     FROM plot_axis_data
     WHERE file_name = ? AND metabolite = ?",
                        params = list(file_name, metabolite)
    )
    
    #peaks_df for this metabolite (used by integration functions)
    peaks_db <- dbGetQuery(con,
                           "SELECT peak_number, start_seconds, apex_seconds, end_seconds,
            start_intensity, apex_intensity, end_intensity, peak_area
     FROM peaks_df
     WHERE file_name = ? AND metabolite = ?
     ORDER BY peak_number",
                           params = list(file_name, metabolite)
    )
    
    list(
      eie_data         = eie_db,
      annotation_data  = ann_db,
      integration_data = intg_db,
      label_data       = c(ax_db$label_name, ax_db$label_mz),
      x_axis_data      = c(ax_db$x_start, ax_db$x_end),
      y_axis_data      = c(ax_db$y_max_annot, ax_db$y_max_eie),
      peaks_df_metab   = peaks_db   
    )
  }
  
  #Write edit log entry function -> maintains a audit trail...
  log_edit <- function(con, file_name, metabolite, edit_type){
    dbExecute(con,
              "INSERT INTO edit_log (file_name, metabolite, edit_type, edited_at) VALUES (?, ?, ?, ?)",
              params = list(file_name, metabolite, edit_type,
                            format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
              )
  }
  
  #Build SQL database from existing .Rdata files. Files already present are skipped.
  migrate_rdata_to_db <- function(results_folder, progress_callback = NULL) {
    
    db_path <- get_db_path(results_folder)
    con     <- dbConnect(RSQLite::SQLite(), db_path)
    initialize_db(con)
    
    data_dir    <- file.path(results_folder, "Data")
    subfolders  <- list.dirs(data_dir, full.names = TRUE, recursive = FALSE)
    rdata_files <- subfolders[file.exists(
      file.path(subfolders, "plot_list_data.RData")
    )]
    
    total   <- length(rdata_files)
    success <- 0
    failed  <- character(0)
    
    for (i in seq_along(rdata_files)) {
      subfolder  <- rdata_files[i]
      folder_name <- basename(subfolder)
      file_name   <- sub("\\.mz5$", "", folder_name)
      if (!is.null(progress_callback)) {
        progress_callback(i, total, file_name)
      }
      
      
      #Check if file is already in database and skip if need be.
      already <- dbGetQuery(con,
                            "SELECT COUNT(*) FROM plot_axis_data WHERE file_name = ?",
                            params = list(file_name)
      )[[1]]
      if (already > 0) {
        success <- success + 1
        next
      }
      
      env <- new.env()
      
      #Loop over files, marking files that failed to load and skipping them.
      tryCatch({
        load(file.path(subfolder, "plot_list_data.RData"), envir = env)
        write_file_to_db(
          con          = con,
          file_name    = file_name,
          plot_list    = env$plot_list,
          peaks_df     = env$peaks_df,
          peak_area_df = env$peak_area_df,
          peak_mt_df   = env$peak_mt_df,
          comment_df   = env$comment_df,
          mi_df        = env$mi_df
        )
        success <- success + 1
      }, error = function(e) {
        failed <<- c(failed, file_name)
        warning("Migration failed for ", file_name, ": ", conditionMessage(e))
      })
      rm(env); gc()
    }
    
    dbDisconnect(con)
    list(total = total, success = success, failed = failed)
  }
  
  
  #Load databases from SQL database for the selected file in wide format -> written to keep original logic
  load_wide_data_from_db <- function(con, file_name) {
    
    #peaks_df: pivot long -> wide. 7 columns per metabolite
    peaks_long <- dbGetQuery(con,
                             "SELECT * FROM peaks_df WHERE file_name = ? ORDER BY metabolite, peak_number",
                             params = list(file_name)
    )
    metabolites <- unique(peaks_long$metabolite)
    num_peaks   <- max(peaks_long$peak_number)
    
    #Make datamatix and fill column by column, matching column names to original peaks_df.
    peaks_wide <- matrix(NA_real_, nrow = num_peaks,
                         ncol = length(metabolites) * 7)
    col_names  <- character(length(metabolites) * 7)
    for (j in seq_along(metabolites)) {
      m   <- metabolites[j]
      sub <- peaks_long[peaks_long$metabolite == m, ]
      base <- (j - 1) * 7
      peaks_wide[sub$peak_number, base + 1] <- sub$start_seconds
      peaks_wide[sub$peak_number, base + 2] <- sub$apex_seconds
      peaks_wide[sub$peak_number, base + 3] <- sub$end_seconds
      peaks_wide[sub$peak_number, base + 4] <- sub$start_intensity
      peaks_wide[sub$peak_number, base + 5] <- sub$apex_intensity
      peaks_wide[sub$peak_number, base + 6] <- sub$end_intensity
      peaks_wide[sub$peak_number, base + 7] <- suppressWarnings(
        as.numeric(sub$peak_area)
      )
      col_names[base + (1:7)] <- paste0(m, c(".start.seconds", ".apex.seconds",
                                             ".end.seconds", ".start_intensity", ".apex_intensity",
                                             ".end_intensity", ".peak.area"))
    }
    peaks_wide_df           <- as.data.frame(peaks_wide)
    colnames(peaks_wide_df) <- col_names
    
    #peak_area_df: pivot long -> wide
    area_long <- dbGetQuery(con,
                            "SELECT peak_number, metabolite, peak_area FROM peak_area_df
       WHERE file_name = ? ORDER BY peak_number, metabolite",
                            params = list(file_name)
    )
    area_wide <- pivot_wider(
      area_long,
      names_from  = "metabolite",
      values_from = "peak_area"
    )
    area_wide <- cbind(
      "file.name"   = c(file_name, rep("", nrow(area_wide) - 1)),
      "peak.number" = area_wide$peak_number,
      area_wide[, -1]
    )
    
    #peak_mt_df: pivot long -> wide
    mt_long <- dbGetQuery(con,
                          "SELECT peak_number, metabolite, migration_time FROM peak_mt_df
       WHERE file_name = ? ORDER BY peak_number, metabolite",
                          params = list(file_name)
    )
    mt_wide <- pivot_wider(
      mt_long,
      names_from  = "metabolite",
      values_from = "migration_time"
    )
    mt_wide <- cbind(
      "file.name"   = c(file_name, rep("", nrow(mt_wide) - 1)),
      "peak.number" = mt_wide$peak_number,
      mt_wide[, -1]
    )
    
    #comment_df: derived from annotation_data
    ann_long <- dbGetQuery(con,
                           "SELECT peak_number, metabolite, comment FROM annotation_data
       WHERE file_name = ? ORDER BY peak_number, metabolite",
                           params = list(file_name))
    comment_wide <- pivot_wider(
      ann_long,
      names_from  = "metabolite",
      values_from = "comment",
      values_fill = "")
    comment_df <- as.data.frame(comment_wide[, -1])
    rownames(comment_df) <- comment_wide$peak_number
    
    #mi_df
    mi_db <- dbGetQuery(con, "SELECT * FROM mi_df")
    
    list(
      peaks_df     = peaks_wide_df,
      peak_area_df = area_wide,
      peak_mt_df   = mt_wide,
      comment_df   = comment_df,
      mi_df        = mi_db
    )
  }
  
  #Build full plot_list for current file from database + RData eie data. EIE data stored seperatley since we never need to change it, and it makes the stored SQL file smaller
  load_plot_list_from_db <- function(con, results_folder, file_name){
    metabolites <- dbGetQuery(con,
                              "SELECT DISTINCT metabolite FROM plot_axis_data WHERE file_name =?",
                              params = list(file_name)
                              )$metabolite
    plot_list <- lapply(metabolites, function(metab){
      read_plot_data_from_db(con, results_folder, file_name, metab)
    })
    names(plot_list) <- metabolites
    plot_list
  }
  
  
  #####4.2 Environmental initialization: File upload, and folder selection#####
  
  # Select filenames before looking for matching saved analyses.
  observe({
    files <- uploadedmz5()$FileName
    current <- isolate(input$file_selector)
    updateSelectInput(session, "file_selector", choices = files,
      selected = if (length(current) == 1L && current %in% files) current else head(files, 1L))
  })

  selected_results_filename <- reactive({
    files <- uploadedmz5()$FileName
    current <- input$file_selector
    if (!length(files)) return(character())
    if (length(current) == 1L && current %in% files) current else files[1]
  })

  # Inspect the matching filename parent and both earlier folder layouts only.
  # Checking Data/<filename>/plot_list_data.RData does not load raw mz5 contents.
  list_results_folders <- function(filename) {
    if (length(filename) != 1L || !valid_mz5_names(filename)) return(character())
    parent <- tryCatch(mz5_result_parent(filename), error = function(e) NULL)
    if (is.null(parent)) return(character())
    legacy_parent <- file.path(".", basename(parent))
    roots <- unique(c(".", parent, legacy_parent))
    roots <- roots[dir.exists(roots)]
    candidates <- unique(unlist(lapply(roots, function(path) {
      list.dirs(path, full.names = TRUE, recursive = FALSE)
    }), use.names = FALSE))
    candidates <- candidates[grepl("^Results", basename(candidates))]
    candidates <- candidates[file.exists(file.path(candidates, "Data", filename, "plot_list_data.RData"))]
    sort(candidates)
  }

  results_folders_reactive <- reactivePoll(5000, session,
    checkFunc = function() list_results_folders(selected_results_filename()),
    valueFunc = function() list_results_folders(selected_results_filename()))

  output$results_folder_status <- renderText({
    filename <- selected_results_filename()
    if (!length(filename)) return("Select an mz5 filename first in Engine.")
    if (!length(results_folders_reactive())) {
      paste("No saved analyses found for", filename)
    } else paste("Saved analyses matching", filename)
  })

  last_results_filename <- reactiveVal(character())
  observeEvent(list(selected_results_filename(), results_folders_reactive()), {
    filename <- selected_results_filename()
    folders <- results_folders_reactive()
    current <- input$results_folder
    same_file <- identical(filename, last_results_filename())
    last_results_filename(filename)
    choices <- c("Select a saved analysis" = "", stats::setNames(folders, folders))
    updateSelectInput(session, "results_folder", choices = choices,
      selected = if (same_file && length(current) == 1L && current %in% folders) current else "")
  }, ignoreNULL = FALSE)

  db_con <- reactiveVal(NULL)

  clear_loaded_analysis <- function() {
    for (state in list(plot_list_data_values, run_metadata, previous_metadata)) {
      for (key in c("plot_list", "comment_df", "peaks_df", "peak_mt_df", "peak_area_df", "mi_df", "file_name")) {
        state[[key]] <- NULL
      }
    }
    modified_peak_plots$names <- character()
    plot_render_trigger(plot_render_trigger() + 1L)
  }

  #Deal with changing results folder and then load database or make database
  observeEvent(input$results_folder, {
    #Close previous connection if open
    old_con <- db_con()
    if (!is.null(old_con)) {
      tryCatch(dbDisconnect(old_con), error = function(e) NULL)
      db_con(NULL)
    }
    
    clear_loaded_analysis()
    folder <- input$results_folder
    if (length(folder) != 1L || !nzchar(folder)) return()
    filename <- selected_results_filename()
    error <- tryCatch({ check_saved_analysis(folder, filename); NULL },
                      error = function(e) conditionMessage(e))
    if (!is.null(error)) {
      showNotification(error, type = "error", duration = 10)
      return()
    }
    updateSelectInput(session, "file_selector2", choices = uploadedmz5()$FileName, selected = filename)
    updateSelectInput(session, "results_folder2", choices = c("Select a saved analysis" = "",
                      stats::setNames(list_results_folders(filename), list_results_folders(filename))), selected = folder)

    if (db_exists(folder)) {
      #Database exists — open it and we're done
      con <- open_db(folder)
      db_con(con)
      cat("Database connected:", get_db_path(folder), "\n")
      
    } else {
      #No database found = count number of .Rdata files to be made into a SQL database.
      rdata_count <- length(list.dirs(
        file.path(folder, "Data"),
        full.names = TRUE, recursive = FALSE
      ))
      
      #Show modal and button options when no SQL database is found for migrating data over to SQL format.
      showModal(modalDialog(
        title = "Database Not Found",
        tags$p(
          "No database found for this Results folder.",
          tags$br(),
          paste0("Found ", rdata_count," processed data file(s) that can be migrated.")),
        tags$p(
          "Would you like to build the database now from your existing data?",
          tags$br(),
          tags$em("This only needs to be done once per Results folder.")),
        footer = tagList(
          modalButton("No — skip"),
          actionButton("confirm_migration", "Yes — build database", class = "btn-primary")),
        easyClose = FALSE
      ))
    }
  }, ignoreNULL = FALSE)

  # Clear stale editing/undo state only when the selected filename actually changes.
  last_visualization_filename <- reactiveVal(character())
  observeEvent(selected_results_filename(), {
    filename <- selected_results_filename()
    if (identical(filename, last_visualization_filename())) return()
    last_visualization_filename(filename)
    con <- isolate(db_con())
    if (!is.null(con)) tryCatch(dbDisconnect(con), error = function(e) NULL)
    db_con(NULL)
    clear_loaded_analysis()
  }, ignoreNULL = FALSE, priority = 100)

  #Handle migration confirmation
  observeEvent(input$confirm_migration, {
    removeModal()
    req(input$results_folder)
    folder <- input$results_folder
    error <- tryCatch({ check_saved_analysis(folder, selected_results_filename()); NULL },
                      error = function(e) conditionMessage(e))
    if (!is.null(error)) {
      showNotification(error, type = "error")
      return()
    }

    #Show progress modal
    showModal(modalDialog(
      title = "Building Database",
      tags$p("Migrating your existing data to the database..."),
      tags$p(id = "migration_progress", "Starting..."),
      footer = NULL,
      easyClose = FALSE
    ))
    
    withProgress(message = "Building database from existing data", value = 0, {
      result <- migrate_rdata_to_db(
        results_folder    = folder,
        progress_callback = function(i, total, file_name) {
          incProgress(
            amount  = 1 / total,
            detail  = paste("Processing:", file_name,
                            paste0("(", i, "/", total, ")"))
          )
        }
      )
    })
    
    removeModal()
    
    if (length(result$failed) == 0) {
      showNotification(
        paste0("Database built successfully! ",
               result$success, " file(s) migrated."),
        type     = "message",
        duration = 5
      )
    } else {
      showNotification(
        paste0("Migration complete. ", result$success, " succeeded, ",
               length(result$failed), " failed: ",
               paste(result$failed, collapse = ", ")),
        type     = "warning",
        duration = 10
      )
    }
    
    #Open new database with open_db and store connection in reactive variable.
    con <- open_db(folder)
    db_con(con)
    cat("Database built and connected:", get_db_path(folder), "\n")
  })
  
  #Clean up database connection when session ends
  session$onSessionEnded(function() {
    con <- isolate(db_con())
    if (!is.null(con)) {
      tryCatch(dbDisconnect(con), error = function(e) NULL)
    }
  })
  
  #When file selector changes, load all data from database into reactive values
  observeEvent(input$file_selector, {
    req(db_con(), input$file_selector)
    con       <- db_con()
    folder    <- main_folder()
    file_name <- sub("\\.mz5$", "", input$file_selector)

    #Checkto make sure file exists
    file_exists_in_db <- tryCatch({
      count <- dbGetQuery(con,
                          "SELECT COUNT(*) FROM plot_axis_data WHERE file_name = ?",
                          params = list(file_name))[[1]]
      count > 0
    }, error = function(e) FALSE)
    
    if (!file_exists_in_db) {
      cat("File not yet in database, skipping load:", file_name, "\n")
      return()
    }
    
    #Load plot_list (eie from .RData, rest from database)
    plot_list <- tryCatch(
      load_plot_list_from_db(con, folder, file_name),
      error = function(e) {
        showNotification(paste("Error loading data:", conditionMessage(e)),
                         type = "error")
        NULL
      }
    )
    req(!is.null(plot_list))

    #Load wide-format data for integration functions
    wide <- load_wide_data_from_db(con, file_name)

    #Populate reactive values
    plot_list_data_values$plot_list    <- plot_list
    plot_list_data_values$comment_df   <- wide$comment_df
    plot_list_data_values$peaks_df     <- wide$peaks_df
    plot_list_data_values$peak_mt_df   <- wide$peak_mt_df
    plot_list_data_values$peak_area_df <- wide$peak_area_df
    plot_list_data_values$mi_df        <- wide$mi_df
    plot_list_data_values$file_name    <- file_name

    run_metadata$plot_list    <- plot_list
    run_metadata$comment_df   <- wide$comment_df
    run_metadata$peaks_df     <- wide$peaks_df
    run_metadata$peak_mt_df   <- wide$peak_mt_df
    run_metadata$peak_area_df <- wide$peak_area_df
    run_metadata$mi_df        <- wide$mi_df
    run_metadata$file_name    <- input$file_selector

    cat("Loaded from database:", file_name,
        "| metabolites:", length(plot_list), "\n")

    plot_render_trigger(plot_render_trigger() + 1)
  })
  
#Obsereve block for loading data when app first launches when a results folder is selected and files havent been changed. Required for direct data acquisition after app startup and no file switch.
  observeEvent(db_con(), {
    req(db_con(), input$file_selector)
    con       <- db_con()
    folder    <- main_folder()
    file_name <- sub("\\.mz5$", "", input$file_selector)
    
    #Check for datafile
    file_exists_in_db <- tryCatch({
      count <- dbGetQuery(con,
                          "SELECT COUNT(*) FROM plot_axis_data WHERE file_name = ?",
                          params = list(file_name))[[1]]
      count > 0
    }, error = function(e) FALSE)
    
    if (!file_exists_in_db) {
      cat("File not yet in database, skipping load:", file_name, "\n")
      return()
    }
    
    #Load plot list eie_data and editable ddata from database
    plot_list <- tryCatch(
      load_plot_list_from_db(con, folder, file_name),
      error = function(e) {
        showNotification(paste("Error loading data:", conditionMessage(e)),
                         type = "error")
        NULL
      }
    )
    req(!is.null(plot_list))
    
    #load data in wide format for integration and editing functions
    wide <- load_wide_data_from_db(con, file_name)
    
    #Populate values
    plot_list_data_values$plot_list    <- plot_list
    plot_list_data_values$comment_df   <- wide$comment_df
    plot_list_data_values$peaks_df     <- wide$peaks_df
    plot_list_data_values$peak_mt_df   <- wide$peak_mt_df
    plot_list_data_values$peak_area_df <- wide$peak_area_df
    plot_list_data_values$mi_df        <- wide$mi_df
    plot_list_data_values$file_name    <- file_name
    
    run_metadata$plot_list    <- plot_list
    run_metadata$comment_df   <- wide$comment_df
    run_metadata$peaks_df     <- wide$peaks_df
    run_metadata$peak_mt_df   <- wide$peak_mt_df
    run_metadata$peak_area_df <- wide$peak_area_df
    run_metadata$mi_df        <- wide$mi_df
    run_metadata$file_name    <- input$file_selector
    
    cat("Loaded from database (on connection):", file_name,
        "| metabolites:", length(plot_list), "\n")
    
    #Update annotation table rendering with new data
    plot_render_trigger(plot_render_trigger() + 1)
  })
  
  
  #####4.3 Render and display plots and annotation information#####
  
  ######4.3.1 Populate plot table and render selected plot######
  #Populate the plot table
  output$plot_table <- DT::renderDataTable({
    data.frame(Plot = filtered_plot_names())
  }, selection = 'single')

  
  #Render and construct plot on demand from SQL databse data
  output$selected_plot <- renderPlotly({
    plot_render_trigger()
    req(input$plot_table_rows_selected, db_con())
    con       <- db_con()
    folder    <- main_folder()
    file_name <- sub("\\.mz5$", "", input$file_selector)
    
    selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
    metabolite         <- sub(paste0("^", file_name, "_"), "", selected_plot_name)
    
    plot_data <- tryCatch(
      read_plot_data_from_db(con, folder, file_name, metabolite),
      error = function(e) NULL
    )
    req(!is.null(plot_data))
    
    plot <- plot_function_plotly(
      eie_data         = plot_data$eie_data,
      annotation_data  = plot_data$annotation_data,
      integration_data = plot_data$integration_data,
      label_data       = plot_data$label_data,
      x_axis_data      = plot_data$x_axis_data,
      y_axis_data      = plot_data$y_axis_data
    )
    
    #Restore zoom state if one exists -> so that the plot doesnt reset every single time you make an edit.
    if (!is.null(zoom_state$xmin)) {
      plot <- plot %>% layout(
        xaxis = list(range = c(zoom_state$xmin, zoom_state$xmax)),
        yaxis = list(range = c(zoom_state$ymin, zoom_state$ymax))
      )
    }
    
    plot %>%
      layout(dragmode = "zoom") %>%
      config(modeBarButtonsToAdd = list("drawrect")) %>%
      event_register("plotly_click") %>%
      event_register("plotly_relayout")
  })
  
  
  
  ######4.3.2 Functions for dealing with moveable lines on plotlyplots######
  #Define logic for toggling between drag mode and draw mode for integrating.
  drag_mode <- reactiveVal("zoom")
  
  #Button toggle logic
  observeEvent(input$draw_mode, {
    
    #Mode for drawing boxes on plot to highlight regions to integrate between 
    if (input$draw_mode) {
      drag_mode("drawrect")
      print("Mode switched to: drawrect")
      plotlyProxy("selected_plot", session) %>%
        plotlyProxyInvoke("relayout", list(dragmode = "drawrect"))
    } else {
      
      #Regular zoom mode for zooming in on peaks
      drag_mode("zoom")
      print("Mode switched to: zoom")
      plotlyProxy("selected_plot", session) %>%
        plotlyProxyInvoke("relayout", list(dragmode = "zoom"))
    }
  })
  
  #Store zoom state reactively when user zooms so we can regenerate the plots without zooming out -> required for efficiency
  zoom_state <- reactiveValues(
    xmin = NULL, xmax = NULL,
    ymin = NULL, ymax = NULL
  )
  
  #Capture zoom from relayout events for regenerating plot.
  observeEvent(event_data("plotly_relayout"), {
    rd <- event_data("plotly_relayout")
    if (!is.null(rd[["xaxis.range[0]"]])) {
      zoom_state$xmin <- rd[["xaxis.range[0]"]]
      zoom_state$xmax <- rd[["xaxis.range[1]"]]
      zoom_state$ymin <- rd[["yaxis.range[0]"]]
      zoom_state$ymax <- rd[["yaxis.range[1]"]]
    }
    # Reset on double-click (autosize)
    if (!is.null(rd[["xaxis.autorange"]])) {
      zoom_state$xmin <- NULL
      zoom_state$xmax <- NULL
      zoom_state$ymin <- NULL
      zoom_state$ymax <- NULL
    }
  })
  
  #Reset zoom when switching metabolites
  observeEvent(input$plot_table_rows_selected, {
    zoom_state$xmin <- NULL
    zoom_state$xmax <- NULL
    zoom_state$ymin <- NULL
    zoom_state$ymax <- NULL
  }) 
  
  #Deal with events occuring on the plot (drawing boxes)
  observeEvent(event_data("plotly_relayout"), {
    relayout_data <- event_data("plotly_relayout")
    print(paste("Current drag_mode:", drag_mode()))
    print(relayout_data)
    
    #Lets user keep drawing after first rectangle drawn
    if (drag_mode() == "drawrect") {
      print("Re-asserting drawrect mode")
      plotlyProxy("selected_plot", session) %>%
        plotlyProxyInvoke("relayout", list(dragmode = "drawrect"))
    }
    
    #Deal with manual dragging of red and blue lines
    if (!is.null(relayout_data[["shapes[0].x0"]])) {
      line_positions$red <- relayout_data[["shapes[0].x0"]]
      print(paste("Red line updated to:", relayout_data[["shapes[0].x0"]]))
    }
    if (!is.null(relayout_data[["shapes[1].x0"]])) {
      line_positions$blue <- relayout_data[["shapes[1].x0"]]
      print(paste("Blue line updated to:", relayout_data[["shapes[1].x0"]]))
    }
    
    #Draw rectangle to dictate start and end integration points
    if (!is.null(relayout_data[["shapes"]]) && 
        any(relayout_data[["shapes"]]$type == "rect")) {
      
      shapes_df <- relayout_data[["shapes"]]
      rect_row <- shapes_df[shapes_df$type == "rect", ]
      
      x_start <- min(rect_row$x0, rect_row$x1)
      x_end   <- max(rect_row$x0, rect_row$x1)
      
      print(paste("Rectangle captured - x_start:", x_start, "x_end:", x_end))
      
      line_positions$red  <- x_start
      line_positions$blue <- x_end
      
      #Remove rectangle, update red/blue lines
      plotlyProxy("selected_plot", session) %>%
        plotlyProxyInvoke("relayout", list(
          shapes = list(
            list(
              type = "line",
              x0 = x_start, x1 = x_start,
              y0 = 0, y1 = 1,
              xref = "x", yref = "paper",
              line = list(color = "Red", width = 2)
            ),
            list(
              type = "line",
              x0 = x_end, x1 = x_end,
              y0 = 0, y1 = 1,
              xref = "x", yref = "paper",
              line = list(color = "Blue", width = 2)
            )
          )
        ))
    }
  })
  
  
  #Output red and blue line positions.
  output$red_line_position <- renderText({
    paste("Red Line Position (Left Boundary):", line_positions$red)
  })
  output$blue_line_position <- renderText({
    paste("Blue Line Position (Right Boundary):", line_positions$blue)
  })
  

  ######4.3.3 Render the annotation data table associated with selecetd plot######
  output$peak_info_table <- DT::renderDataTable({
    plot_render_trigger()
    req(input$plot_table_rows_selected, db_con())
    con       <- db_con()
    file_name <- sub("\\.mz5$", "", input$file_selector)
    
    #Extract metabolite name
    selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
    metabolite <- sub(paste0("^", file_name, "_"), "", selected_plot_name)
    
    #get annotation data for specified file and metabolite from database
    ann <- dbGetQuery(con,
                      "SELECT peak_number, comment, peak_apex_secs, peak_height
     FROM annotation_data
     WHERE file_name = ? AND metabolite = ?
     ORDER BY peak_number",
                      params = list(file_name, metabolite)
    )
    
    #rename columns
    colnames(ann) <- c("peak.number", "comment", "peak.apex.seconds", "peak.height.counts")
    ann
  }, options = list(pageLength = 13), selection = "multiple")
  
  
  #####4.4 Delete Peaks and Update Metadata#####
  ######4.4.1 Initialize environment, define reactive expressions, and define functions######
  
  #Define reactive expression to load plot_list_data
  loaded_plot_list_data <- reactive({
    req(main_folder(), input$file_selector)
    file_name <- input$file_selector
    subfolder_path <- file.path(main_folder(), "Data", file_name)

    #Debug print statements for paths
    print(paste("Main folder:", main_folder()))
    print(paste("File selector:", file_name))
    print(paste("Subfolder path:", subfolder_path))

    #Check if the file exists before attempting to load it
    if (file.exists(file.path(subfolder_path, "plot_list_data.RData"))) {
      load(file.path(subfolder_path, "plot_list_data.RData"))
      return(list(plot_list = plot_list, peaks_df = peaks_df, peak_mt_df = peak_mt_df, peak_area_df = peak_area_df, comment_df = comment_df, mi_df = mi_df))
    } else {
      warning("File not found: ", file.path(subfolder_path, "plot_list_data.RData"))
      return(NULL)
    }
  })

  #Function for saving and updating raw data aftyer any edit to a plot. Writes edited data for a single metabolite into the SQL database and then rerenders plot with new data. Also keeps writing .Rdata as a backup.
  save_plot_data <- function(plotname) {
    req(db_con())
    con       <- db_con()
    file_name <- plot_list_data_values$file_name
    folder    <- main_folder()
    
    #Update existing annotation data
    ann <- plot_list_data_values$plot_list[[plotname]]$annotation_data
    dbExecute(con,
              "DELETE FROM annotation_data WHERE file_name = ? AND metabolite = ?",
              params = list(file_name, plotname)
    )
    ann_db <- data.frame(
      file_name      = file_name,
      metabolite     = plotname,
      peak_number    = ann$peak.number,
      comment        = as.character(ann$comment),
      peak_apex_secs = ann$peak.apex.seconds,
      peak_height    = ann$peak.height.counts,
      stringsAsFactors = FALSE
    )
    dbAppendTable(con, "annotation_data", ann_db)
    
    #Update integration data
    intg <- plot_list_data_values$plot_list[[plotname]]$integration_data
    dbExecute(con,
              "DELETE FROM integration_data WHERE file_name = ? AND metabolite = ?",
              params = list(file_name, plotname)
    )
    if (!is.null(intg) && nrow(intg) > 0) {
      intg$peak.number <- as.numeric(as.character(intg$peak.number))
      intg_db <- data.frame(
        file_name   = file_name,
        metabolite  = plotname,
        peak_number = intg$peak.number,
        mt_seconds  = intg$mt.seconds,
        intensity   = intg$intensity,
        baseline    = intg$baseline,
        stringsAsFactors = FALSE
      )
      dbAppendTable(con, "integration_data", intg_db)
    }
    
    #Update peaks_df
    peaks_df  <- plot_list_data_values$peaks_df
    num_peaks <- nrow(peaks_df)
    dbExecute(con,
              "DELETE FROM peaks_df WHERE file_name = ? AND metabolite = ?",
              params = list(file_name, plotname)
    )
    peaks_db <- data.frame(
      file_name       = file_name,
      metabolite      = plotname,
      peak_number     = 1:num_peaks,
      start_seconds   = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(plotname, ".start.seconds")))),
      apex_seconds    = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(plotname, ".apex.seconds")))),
      end_seconds     = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(plotname, ".end.seconds")))),
      start_intensity = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(plotname, ".start_intensity")))),
      apex_intensity  = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(plotname, ".apex_intensity")))),
      end_intensity   = suppressWarnings(as.numeric(safe_col(peaks_df, paste0(plotname, ".end_intensity")))),
      peak_area       = as.character(safe_col(peaks_df, paste0(plotname, ".peak.area"))),
      stringsAsFactors = FALSE
    )
    dbAppendTable(con, "peaks_df", peaks_db)
    
    #Update peak area entries
    peak_area_df <- plot_list_data_values$peak_area_df
    dbExecute(con,"DELETE FROM peak_area_df WHERE file_name = ? AND metabolite = ?",params = list(file_name, plotname))
    
    area_db <- data.frame(
      file_name   = file_name,
      peak_number = 1:num_peaks,
      metabolite  = plotname,
      peak_area   = as.character(peak_area_df[, plotname]),
      stringsAsFactors = FALSE
    )
    dbAppendTable(con, "peak_area_df", area_db)
    
    #Record data edit in edit_log trail.
    log_edit(con, file_name, plotname, "save")
    
    #Write and update .Rdata in the datafolder as backup
    subfolder_path <- file.path(folder, "Data", paste0(file_name, ".mz5"))
    if (dir.exists(subfolder_path)) {
      plot_list    <- plot_list_data_values$plot_list
      comment_df   <- plot_list_data_values$comment_df
      peaks_df     <- plot_list_data_values$peaks_df
      peak_mt_df   <- plot_list_data_values$peak_mt_df
      peak_area_df <- plot_list_data_values$peak_area_df
      mi_df        <- plot_list_data_values$mi_df
      save(plot_list, comment_df, peaks_df, peak_mt_df, mi_df, peak_area_df,
           file = file.path(subfolder_path, "plot_list_data.RData"))
    }
    
    #Update reactive values with new metabolite data.
    updated_plot_data <- tryCatch(
      read_plot_data_from_db(con, folder, file_name, plotname),
      error = function(e) NULL)
    if (!is.null(updated_plot_data)) {
      plot_list_data_values$plot_list[[plotname]] <- updated_plot_data
      run_metadata$plot_list[[plotname]]          <- updated_plot_data
    }
    
    #Force plot and annotation table to re-render
    plot_render_trigger(plot_render_trigger() + 1)
  }
  
  
  #reload_plot_data: now just re-reads from database into reactive values. Kept for backward compatibility with any code that calls it
  reload_plot_data <- function(subfolder_path = NULL) {
    req(db_con())
    con       <- db_con()
    folder    <- main_folder()
    file_name <- plot_list_data_values$file_name
    req(!is.null(file_name))
    
    wide <- load_wide_data_from_db(con, file_name)
    plot_list_data_values$comment_df   <- wide$comment_df
    plot_list_data_values$peaks_df     <- wide$peaks_df
    plot_list_data_values$peak_mt_df   <- wide$peak_mt_df
    plot_list_data_values$peak_area_df <- wide$peak_area_df
    run_metadata$comment_df            <- wide$comment_df
    run_metadata$peaks_df              <- wide$peaks_df
    run_metadata$peak_mt_df            <- wide$peak_mt_df
    run_metadata$peak_area_df          <- wide$peak_area_df
  }
  
  
  
  #Intiialize reactive values 
  observeEvent(input$file_selector, {
    req(loaded_plot_list_data())
    
    #Initialize reactive values only once when a new file is selected
    plot_list_data_values$plot_list <- loaded_plot_list_data()$plot_list
    plot_list_data_values$comment_df <- loaded_plot_list_data()$comment_df
    plot_list_data_values$peaks_df <- loaded_plot_list_data()$peaks_df
    plot_list_data_values$peak_mt_df <- loaded_plot_list_data()$peak_mt_df
    plot_list_data_values$peak_area_df <- loaded_plot_list_data()$peak_area_df
    plot_list_data_values$mi_df <- loaded_plot_list_data()$mi_df
    
    #Store metadata in a separate reactive structure to allow users to access/undo peak deletions
    run_metadata$comment_df <- plot_list_data_values$comment_df
    run_metadata$plot_list <- plot_list_data_values$plot_list
    run_metadata$file_name <- input$file_selector
    run_metadata$peaks_df <- plot_list_data_values$peaks_df
    run_metadata$peak_mt_df <- plot_list_data_values$peak_mt_df
    run_metadata$peak_area_df <- plot_list_data_values$peak_area_df
    run_metadata$mi_df <- plot_list_data_values$mi_df
  })
  
  ######4.4.2 Server logic for deleting peaks in plots######
  observeEvent(input$delete_peak, {
    req(input$peak_info_table_rows_selected)
    
    #Store the current state before making changes
    previous_metadata$comment_df <- run_metadata$comment_df
    previous_metadata$plot_list <- run_metadata$plot_list
    previous_metadata$peaks_df <- run_metadata$peaks_df
    previous_metadata$peak_mt_df <- run_metadata$peak_mt_df
    previous_metadata$peak_area_df <- run_metadata$peak_area_df
    previous_metadata$mi_df <- run_metadata$mi_df
    
    selected_rows <- input$peak_info_table_rows_selected
    
    #Extract the current annotation_data
    selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
    base_file_name <- sub("\\.mz5$", "", input$file_selector)
    plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
    plot_list <- plot_list_data_values$plot_list
    
    #Mark the plot as edited so only plots that have been edited can be regenerated later
    modified_peak_plots$names <- unique(c(modified_peak_plots$names, selected_plot_name))
    
    if (plotname %in% names(plot_list)) {
      annotation_data <- plot_list[[plotname]]$annotation_data
      integration_data <- plot_list[[plotname]]$integration_data
      
      #Loop through each selected row
      for (selected_row in selected_rows) {
        #Get the peak number from the selected row to allow matching accross various acccesible data streams
        peak_number <- annotation_data$peak.number[selected_row]
        
        #Update the comment column for the selected peak
        annotation_data$comment[selected_row] <- "NPD"
        
        #Update comment_df so peak areas can be recauculated/reintegrated later on
        if (plotname %in% colnames(plot_list_data_values$comment_df)) {
          plot_list_data_values$comment_df[peak_number, plotname] <- "NPD"
        }
        
        if (plotname %in% colnames(plot_list_data_values$peak_area_df)) {
          plot_list_data_values$peak_area_df[peak_number, plotname] <- "NPD"
        }
        
        #Remove all integration_data associated with the relevant peak number so regenerated plots will update with the deleted peaks once regenerated
        integration_data <- integration_data[integration_data$peak.number != peak_number, ]
      }
      
      #Update the plot_list with the modified annotation_data
      plot_list[[plotname]]$annotation_data <- annotation_data
      plot_list[[plotname]]$integration_data <- integration_data
      
      #Update the reactive plot_list_data
      plot_list_data_values$plot_list <- plot_list
    }
    
    #Save changes 
    save_plot_data(plotname)
    
    #Reload the updated data
    subfolder_path <- file.path(input$results_folder, "Data", run_metadata$file_name)
    reload_plot_data(subfolder_path)
    
    #Update reactive values directly
    plot_list_data_values$plot_list <- plot_list
    plot_list_data_values$comment_df <- plot_list_data_values$comment_df
    
    
    #Update the peak_info_table
    output$peak_info_table <- DT::renderDataTable({
      req(input$plot_table_rows_selected)
      selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
      
      base_file_name <- sub("\\.mz5$", "", input$file_selector)
      plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
      plot_list <- plot_list_data_values$plot_list
      
      #Ensure the plotname exists in the plot_list
      if (plotname %in% names(plot_list)) {
        annotation_data <- plot_list[[plotname]]$annotation_data
        annotation_data
      } else {
        data.frame()
      }
    }, options = list(pageLength = 13), selection = 'multiple')
})
  

  ######4.4.3 Button for deleting all but the first and last peaks######
  observeEvent(input$deletemiddlepeaks, {
    req(input$plot_table_rows_selected)
    
    #Store the current state before making changes
    previous_metadata$comment_df <- run_metadata$comment_df
    previous_metadata$plot_list <- run_metadata$plot_list
    previous_metadata$peaks_df <- run_metadata$peaks_df
    previous_metadata$peak_mt_df <- run_metadata$peak_mt_df
    previous_metadata$peak_area_df <- run_metadata$peak_area_df
    previous_metadata$mi_df <- run_metadata$mi_df
    
    #Extract the current plot name
    selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
    base_file_name <- sub("\\.mz5$", "", input$file_selector)
    plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
    plot_list <- plot_list_data_values$plot_list
    
    #Mark the plot as edited so only plots that have been edited can be regenerated later
    modified_peak_plots$names <- unique(c(modified_peak_plots$names, selected_plot_name))
    
    if (plotname %in% names(plot_list)) {
      annotation_data <- plot_list[[plotname]]$annotation_data
      integration_data <- plot_list[[plotname]]$integration_data
      
      #Identify the first and last peak numbers
      first_peak <- annotation_data$peak.number[1]
      last_peak <- annotation_data$peak.number[nrow(annotation_data)]
      
      #Update the comment column for all peaks except the first and last
      annotation_data$comment[-c(1, nrow(annotation_data))] <- "NPD"
      
      #Update comment_df so peak areas can be recalculated/reintegrated later on
      if (plotname %in% colnames(plot_list_data_values$comment_df)) {
        plot_list_data_values$comment_df[-c(first_peak, last_peak), plotname] <- "NPD"
      }
      
      if (plotname %in% colnames(plot_list_data_values$peak_area_df)) {
        plot_list_data_values$peak_area_df[-c(first_peak, last_peak), plotname] <- "NPD"
      }
      
      #Remove all integration_data associated with peaks except the first and last
      integration_data <- integration_data[integration_data$peak.number %in% c(first_peak, last_peak), ]
      
      #Update the plot_list with the modified annotation_data
      plot_list[[plotname]]$annotation_data <- annotation_data
      plot_list[[plotname]]$integration_data <- integration_data
      
      #Update the reactive plot_list_data
      plot_list_data_values$plot_list <- plot_list
      
      #Save changes
      save_plot_data(plotname)
      
      #Reload the updated data
      subfolder_path <- file.path(input$results_folder, "Data", run_metadata$file_name)
      reload_plot_data(subfolder_path)
      
      #Update reactive values directly
      plot_list_data_values$plot_list <- plot_list
      plot_list_data_values$comment_df <- plot_list_data_values$comment_df
      
      #Update the peak_info_table
      output$peak_info_table <- DT::renderDataTable({
        req(input$plot_table_rows_selected)
        selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
        
        base_file_name <- sub("\\.mz5$", "", input$file_selector)
        plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
        plot_list <- plot_list_data_values$plot_list
        
        #Ensure the plotname exists in the plot_list
        if (plotname %in% names(plot_list)) {
          annotation_data <- plot_list[[plotname]]$annotation_data
          annotation_data
        } else {
          data.frame()
        }
      }, options = list(pageLength = 13), selection = 'multiple')
    }
  })
  
  ######4.4.4 Button for deleting ALL peaks in a plot######
  observeEvent(input$deleteallpeaks, {
    req(input$plot_table_rows_selected)
    
    #Store the current state before making changes
    previous_metadata$comment_df <- run_metadata$comment_df
    previous_metadata$plot_list <- run_metadata$plot_list
    previous_metadata$peaks_df <- run_metadata$peaks_df
    previous_metadata$peak_mt_df <- run_metadata$peak_mt_df
    previous_metadata$peak_area_df <- run_metadata$peak_area_df
    previous_metadata$mi_df <- run_metadata$mi_df
    
    #Extract the current plot name
    selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
    base_file_name <- sub("\\.mz5$", "", input$file_selector)
    plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
    plot_list <- plot_list_data_values$plot_list
    
    #Mark the plot as edited so only plots that have been edited can be regenerated later
    modified_peak_plots$names <- unique(c(modified_peak_plots$names, selected_plot_name))
    
    if (plotname %in% names(plot_list)) {
      annotation_data <- plot_list[[plotname]]$annotation_data
      integration_data <- plot_list[[plotname]]$integration_data
      
      #Update the comment column for all peaks
      annotation_data$comment <- "NPD"
      
      #Update comment_df so peak areas can be recalculated/reintegrated later on
      if (plotname %in% colnames(plot_list_data_values$comment_df)) {
        plot_list_data_values$comment_df[, plotname] <- "NPD"
      }
      
      if (plotname %in% colnames(plot_list_data_values$peak_area_df)) {
        plot_list_data_values$peak_area_df[, plotname] <- "NPD"
      }
      
      #Remove all integration_data
      integration_data <- integration_data[0, ]
      
      #Update the plot_list with the modified annotation_data
      plot_list[[plotname]]$annotation_data <- annotation_data
      plot_list[[plotname]]$integration_data <- integration_data
      
      #Update the reactive plot_list_data
      plot_list_data_values$plot_list <- plot_list
      
      #Save changes
      save_plot_data(plotname)
      
      #Reload the updated data
      subfolder_path <- file.path(input$results_folder, "Data", run_metadata$file_name)
      reload_plot_data(subfolder_path)
      
      #Update reactive values directly
      plot_list_data_values$plot_list <- plot_list
      plot_list_data_values$comment_df <- plot_list_data_values$comment_df
      
      #Update the peak_info_table
      output$peak_info_table <- DT::renderDataTable({
        req(input$plot_table_rows_selected)
        selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
        
        base_file_name <- sub("\\.mz5$", "", input$file_selector)
        plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
        plot_list <- plot_list_data_values$plot_list
        
        #Ensure the plotname exists in the plot_list
        if (plotname %in% names(plot_list)) {
          annotation_data <- plot_list[[plotname]]$annotation_data
          annotation_data
        } else {
          data.frame()
        }
      }, options = list(pageLength = 13), selection = 'multiple')
    }
  })
  
     
  #####4.5 Undo peak deletion#####
  observeEvent(input$undo, {
    #Revert to the previous state
    run_metadata$comment_df <- previous_metadata$comment_df
    run_metadata$plot_list <- previous_metadata$plot_list
    run_metadata$peaks_df <- previous_metadata$peaks_df
    run_metadata$peak_mt_df <- previous_metadata$peak_mt_df
    run_metadata$peak_area_df <- previous_metadata$peak_area_df
    run_metadata$mi_df <- previous_metadata$mi_df

    #Store metadata in a separate reactive structure to allow users to access/undo peak deletions
    plot_list_data_values$comment_df <- run_metadata$comment_df 
    plot_list_data_values$plot_list <- run_metadata$plot_list
    plot_list_data_values$peaks_df <- run_metadata$peaks_df
    plot_list_data_values$peak_mt_df <- run_metadata$peak_mt_df
    plot_list_data_values$peak_area_df <- run_metadata$peak_area_df
    plot_list_data_values$mi_df <- run_metadata$mi_df
    
    #Extract objects from run_metadata
    plot_list <- run_metadata$plot_list
    comment_df <- run_metadata$comment_df
    peaks_df <- run_metadata$peaks_df
    peak_mt_df <- run_metadata$peak_mt_df
    peak_area_df <- run_metadata$peak_area_df
    mi_df <- run_metadata$mi_df
    
    #Save the reverted state back to the file
    subfolder_path <- file.path(input$results_folder, "Data", run_metadata$file_name)
    save(plot_list, comment_df, mi_df,  peaks_df, peak_mt_df, peak_area_df, mi_df, file = file.path(subfolder_path, "plot_list_data.RData"))
    
    #Update the UI to reflect the reverted state
    output$peak_info_table <- DT::renderDataTable({
      req(input$plot_table_rows_selected)
      selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
      base_file_name <- sub("\\.mz5$", "", run_metadata$file_name)
      plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
      plot_list <- run_metadata$plot_list
      
      if (plotname %in% names(plot_list)) {
        annotation_data <- plot_list[[plotname]]$annotation_data
        annotation_data
      } else {
        data.frame()  
      }
    }, options = list(pageLength = 13), selection = 'multiple')
  })
  
  #####4.6 Regenerate plots#####

  ######4.6.1 Defining Plotting Functions######
   #ggplot function
   plot_function <- function(eie_data, annotation_data, integration_data, label_data, x_axis_data, y_axis_data, base_file_name, font_size_1, font_size_2){
  
     extra_space <- ifelse(x_axis_data[1] > 70, 1, 0)
  
     #Convert peak.number to factor
     integration_data$peak.number <- as.factor(integration_data$peak.number)
     
     ggplot(data = eie_data) +
      geom_line(aes(x = mt.seconds / 60, y = eie_data[,2]), colour = "grey50") +
       theme_classic() +
       coord_cartesian(xlim = c(x_axis_data[1] / 60 - extra_space, x_axis_data[2] / 60 + extra_space),
                       ylim = c(0, 1.5 * y_axis_data[1])) +
       scale_y_continuous(name = "Ion Counts",
                          labels = function(x) format(x, scientific = TRUE),
                          expand = c(0,0),
                          breaks = scales::pretty_breaks(n = 10)) +
       scale_x_continuous(name = "Migration Time (Minutes)",
                          breaks = scales::pretty_breaks(n = 10))+
       ggtitle(paste(label_data[1], " EIE", " (m/z = ", label_data[2],")",sep = ""),
              subtitle = paste("Data File: ", base_file_name)) +
       geom_ribbon(data = integration_data,
                   aes(x = mt.seconds/60, ymax = intensity, ymin = baseline, fill = peak.number),
                   alpha =0.4) +
       geom_text(data = annotation_data,
                 label = annotation_data$peak.number,
                 size  = font_size_1,
                 family = "sans",
                 aes(x = peak.apex.seconds/60,
                     y = peak.height.counts + 0.1 * y_axis_data[1])) +
       geom_text(data = annotation_data,
                 label = annotation_data$comment,
                 size  = font_size_1,
                 family = "sans",
                 aes(x = peak.apex.seconds/60,
                     y = peak.height.counts + 0.2 * y_axis_data[1])) +
       theme(legend.position = "none",
             text = element_text(size = font_size_2, family = "sans"))
   }
  
   #plotly function
   plot_function_plotly <- function(eie_data, annotation_data, integration_data, label_data, x_axis_data, y_axis_data) {
     plot_ly() %>%
       add_lines(
         x = ~eie_data[,1] / 60,
         y = ~eie_data[,2],
         name = 'Electropherogram',
         line = list(color = 'grey')
       ) %>%
       add_ribbons(
         x = ~integration_data$mt.seconds / 60,
         ymin = ~integration_data$baseline,
         ymax = ~integration_data$intensity,
         fillcolor = 'rgba(100,100,255,0.4)',
         line = list(color = 'rgba(100,100,255,0.4)'),
         name = 'Peak Integration'
       ) %>%
       add_text(
         x = ~annotation_data$peak.apex.seconds / 60,
         y = ~annotation_data$peak.height.counts + (0.1 * y_axis_data[1]),
         text = ~annotation_data$peak.number,
         showlegend = FALSE
       ) %>%
       add_text(
         x = ~annotation_data$peak.apex.seconds / 60,
         y = ~annotation_data$peak.height.counts + (0.2 * y_axis_data[1]),
         text = ~annotation_data$comment,
         showlegend = FALSE
       ) %>%
       layout(
         title = paste(label_data[1], " EIE", " (m/z =", label_data[2], ")"),
         xaxis = list(title = 'Migration Time (Minutes)'),
         yaxis = list(title = 'Ion Counts'),
         template = 'plotly_white',
         shapes = list(
           list(
             type = "line",
             x0 = 2,  #Initial position of the red line
             y0 = 0,
             x1 = 2,
             y1 = 1,
             xref = "x",
             yref = "paper",
             line = list(color = "Red", width = 2)
           ),
           list(
             type = "line",
             x0 = 4,  #Initial position of the blue line
             y0 = 0,
             x1 = 4,
             y1 = 1,
             xref = "x",
             yref = "paper",
             line = list(color = "Blue", width = 2)
           )
         )
       ) %>% config(editable = TRUE)
   }
  
   ######4.6.2 Defining functions for modifying saved peak area information######
   regenerate_metabolite_peak_areas <- function() {
     req(db_con())
     con       <- db_con()
     file_name <- plot_list_data_values$file_name
     
     #Pull all peak areas from database across all files in this Results folder
     area_long <- dbGetQuery(con,
                             "SELECT file_name, peak_number, metabolite, peak_area FROM peak_area_df
       ORDER BY file_name, peak_number"
     )
     
     #Pivot to wide format matching original CSV structure
     area_wide <- tidyr::pivot_wider(
       area_long,
       names_from  = "metabolite",
       values_from = "peak_area"
     )
     area_wide <- cbind(
       "file.name"   = area_wide$file_name,
       "peak.number" = area_wide$peak_number,
       area_wide[, !(colnames(area_wide) %in% c("file_name", "peak_number"))]
     )
     
     write.csv(area_wide,
               file      = file.path(main_folder(), "Metabolite Peak Areas.csv"),
               row.names = FALSE)
   }
   
  ######4.6.3 Function for regenerating plotly plots######
   regenerate_plots <- function() {
     #Define font sizes for plotting functions
     font_size_1 <- 7
     font_size_2 <- 25
     
     #Get the list of all folders in the Data directory
     data_subfolders <- list.dirs(path = file.path(input$results_folder, "Data"), full.names = TRUE, recursive = FALSE)
     
     #Initialize an empty list to store all plotly objects
     plotly_objects <- list()
     
     for (subfolder_path in data_subfolders) {
       #Extract the base file name from the subfolder name
       base_file_name <- sub("\\.mz5$", "", basename(subfolder_path))
       
       #Load the plot_list_data.RData file
       if (file.exists(file.path(subfolder_path, "plot_list_data.RData"))) {
         load(file.path(subfolder_path, "plot_list_data.RData"))
         
         for (plotname in names(plot_list)) {
           #Print the name of the plot being worked on
           print(paste("Regenerating plot:", plotname))
           
           #Retrieve the plot data for the current plotname
           plot_data <- plot_list[[plotname]]
           
           plot <- plot_function_plotly(
             eie_data = plot_data$eie_data,
             annotation_data = plot_data$annotation_data,
             integration_data = plot_data$integration_data,
             label_data = plot_data$label_data,
             x_axis_data = plot_data$x_axis_data,
             y_axis_data = plot_data$y_axis_data
           )
           
           plot <- event_register(plot, "plotly_click")
           
           #This is necessary to save plots correctly
           plotly_copy <- plotly::plotly_build(plot)
           
           #Name the plot to match with existing plots
           plot_name <- paste0(base_file_name, "_", plotname)
           
           #Update or add the regenerated plot in the plotly_objects list
           plotly_objects[[plot_name]] <- plotly_copy
         }
         plotly_data(plotly_objects)
       }
     }
     
     
     #Save all plotly objects to a single .RData file
     save(plotly_objects, file = file.path(input$results_folder, "plotly_objects.RData"), compress = "xz")
     print("Finished regenerating and compressing all plots!")
   }

   ######4.6.4 Function for regenerating ggplots plots and saving them to specified locations######
   regenerate_ggplot <- function() {
     #Define plotting housekeeping variables
     plot_format <- parametersData()$plot.format
     name_vec <- c(refMassListData()$name, massData()$name)
     num_of_metabolites <- nrow(massData())
     num_of_is <- nrow(refMassListData())
     
     #Get the list of all folders in the Data directory
     data_subfolders <- list.dirs(path = file.path(input$results_folder, "Data"), full.names = TRUE, recursive = FALSE)
     
     for (subfolder_path in data_subfolders) {
       #Extract the base file name from the subfolder name and remove the .mz5 extension
       base_file_name <- sub("\\.mz5$", "", basename(subfolder_path))
       
       #Load the plot_list_data.RData file
       if (file.exists(file.path(subfolder_path, "plot_list_data.RData"))) {
         load(file.path(subfolder_path, "plot_list_data.RData"))
         
         #Load data for plotting along with mi_df to determine which metabolites use MI and which use RMTs
         plot_list <- get("plot_list")
         mi_df <- get("mi_df")
         
         modified_plot_names <- sapply(modified_peak_plots$names, function(name) {
           sub(paste0("^", base_file_name, "_"), "", name)
         })
         
         modified_plots <- intersect(names(plot_list), modified_plot_names)
         
         
         #Regenerate plots according to "Sample" format
         if (plot_format == "Sample") {
           for (plotname in modified_plots) {
             
             print(paste("Regenerating plot:", plotname, "from file:", base_file_name))
             
             #Load plot data
             plot_data <- plot_list[[plotname]]
             
             #Define index for correctly naming saved files
             plot_index <- which(name_vec == plotname)
             
             #Plot both internal standard and analyte plots
             if (plotname %in% name_vec[1:num_of_is]) {
               font_size_1 <- 7
               font_size_2 <- 25
               
               #Save internal standard plots
               folder <- "Internal Standards"
               ggsave(filename = paste(plot_index, "_", plotname, ".png", sep = ""),
                      width = 16,
                      height = 9,
                      plot = plot_function(eie_data = plot_data$eie_data,
                                           annotation_data = plot_data$annotation_data,
                                           integration_data = plot_data$integration_data,
                                           label_data = plot_data$label_data,
                                           x_axis_data = plot_data$x_axis_data,
                                           y_axis_data = plot_data$y_axis_data,
                                           base_file_name = base_file_name,
                                           font_size_1 = font_size_1,
                                           font_size_2 = font_size_2),
                      path = paste(input$results_folder, "/Plots/", folder, "/", base_file_name, sep = ""))
             } else {
               folder <- "Analytes"
               
               #Analytes using RMT
               if (mi_df$description[which(name_vec == plotname) - num_of_is] != "mi") {
                 font_size_1 <- 4
                 font_size_2 <- 12
                 is_index <- which(name_vec == mi_df$description[which(name_vec == plotname) - num_of_is])
                 
                 figure <- ggarrange(plot_function(eie_data = plot_data$eie_data,
                                                   annotation_data = plot_data$annotation_data,
                                                   integration_data = plot_data$integration_data,
                                                   label_data = plot_data$label_data,
                                                   x_axis_data = c(min(plot_data$eie_data$mt.seconds), max(plot_data$eie_data$mt.seconds)),
                                                   y_axis_data = plot_data$y_axis_data,
                                                   base_file_name = base_file_name,
                                                   font_size_1 = font_size_1,
                                                   font_size_2 = font_size_2),
                                     plot_function(eie_data = plot_data$eie_data,
                                                   annotation_data = plot_data$annotation_data,
                                                   integration_data = plot_data$integration_data,
                                                   label_data = plot_data$label_data,
                                                   x_axis_data = plot_data$x_axis_data,
                                                   y_axis_data = plot_data$y_axis_data,
                                                   base_file_name = base_file_name,
                                                   font_size_1 = font_size_1,
                                                   font_size_2 = font_size_2),
                                     plot_function(eie_data = plot_list[[is_index]]$eie_data,
                                                   annotation_data = plot_list[[is_index]]$annotation_data,
                                                   integration_data = plot_list[[is_index]]$integration_data,
                                                   label_data = plot_list[[is_index]]$label_data,
                                                   x_axis_data = plot_list[[is_index]]$x_axis_data,
                                                   y_axis_data = plot_list[[is_index]]$y_axis_data,
                                                   base_file_name = base_file_name,
                                                   font_size_1 = font_size_1,
                                                   font_size_2 = font_size_2),
                                     ncol = 1, nrow = 3)
                 ggsave(filename = paste(plot_index, "_", plotname, ".png", sep = ""),
                        width = 16,
                        height = 9,
                        plot = figure,
                        path = paste(input$results_folder, "/Plots/", folder, "/", base_file_name, sep = ""))
               } else {
                 font_size_1 <- 4
                 font_size_2 <- 12
                 figure <- ggarrange(plot_function(eie_data = plot_data$eie_data,
                                                   annotation_data = plot_data$annotation_data,
                                                   integration_data = plot_data$integration_data,
                                                   label_data = plot_data$label_data,
                                                   x_axis_data = c(min(plot_data$eie_data$mt.seconds), max(plot_data$eie_data$mt.seconds)),
                                                   y_axis_data = plot_data$y_axis_data,
                                                   base_file_name = base_file_name,
                                                   font_size_1 = font_size_1,
                                                   font_size_2 = font_size_2),
                                     plot_function(eie_data = plot_data$eie_data,
                                                   annotation_data = plot_data$annotation_data,
                                                   integration_data = plot_data$integration_data,
                                                   label_data = plot_data$label_data,
                                                   x_axis_data = plot_data$x_axis_data,
                                                   y_axis_data = plot_data$y_axis_data,
                                                   base_file_name = base_file_name,
                                                   font_size_1 = font_size_1,
                                                   font_size_2 = font_size_2),
                                     ncol = 1, nrow = 2)
                 ggsave(filename = paste(plot_index, "_", plotname, ".png", sep = ""),
                        width = 16,
                        height = 9,
                        plot = figure,
                        path = paste(input$results_folder, "/Plots/", folder, "/", base_file_name, sep = ""))
               }
             }
           }
         }
         
         #Regenerate plots according to "Metabolite" format
         if (plot_format == "Metabolite") {
           for (n in 1:length(name_vec)) {
             plotname <- name_vec[n]
             if (plotname %in% modified_plot_names) {
               folder <- plotname
               ggsave(filename = paste(base_file_name, ".png", sep = ""),
                      width = 16,
                      height = 9,
                      plot = plot_function(eie_data = plot_list[[n]]$eie_data,
                                           annotation_data = plot_list[[n]]$annotation_data,
                                           integration_data = plot_list[[n]]$integration_data,
                                           label_data = plot_list[[n]]$label_data,
                                           x_axis_data = plot_list[[n]]$x_axis_data,
                                           y_axis_data = plot_list[[n]]$y_axis_data),
                      path = paste(input$results_folder, "/Plots/", folder, "/", sep = ""))
             }
           }
         }
       }
     }
   }
   
   ######4.6.5 Function for Regenerating Migration time file#####
   regenerate_metabolite_migration_times <- function() {
     req(db_con())
     con <- db_con()
     
     # Pull all apex times from annotation_data across all files
     mt_long <- dbGetQuery(con,
                           "SELECT file_name, peak_number, metabolite, peak_apex_secs
     FROM annotation_data
     ORDER BY file_name, peak_number"
     )
     
     # Pivot to wide format
     mt_wide <- tidyr::pivot_wider(
       mt_long,
       names_from  = "metabolite",
       values_from = "peak_apex_secs"
     )
     mt_wide <- cbind(
       "file.name"   = mt_wide$file_name,
       "peak.number" = mt_wide$peak_number,
       mt_wide[, !(colnames(mt_wide) %in% c("file_name", "peak_number"))]
     )
     
     write.csv(mt_wide,
               file      = file.path(main_folder(), "Metabolite Migration Times.csv"),
               row.names = FALSE)
   }


   ######4.6.6 Applying functions######
#Loading messages for loading screen
   loading_messages <- c(
     "Training squirrels to redraw plots...", 
     "Bribing CE ghosts to cooperate...",
     "Praying to the R gods that the code will work...",
     "Plotting revenge against noisy data...",
     "Summoning peak integration spirits...",
     "Sacrificing the old plots for better results...",
     "Convincing peaks they really do exist...",
     "Making a blood sacrifice to the R gods...",
     "Offering a peace treaty to signal noise...",
     "Selling soul for a perfectly integrated peak...",
     "Syncing plotting preferences with your Spotify...",
     "Telling algorithm to 'just work'...",
     "Praying the plots havent unionized...",
     "Debugging... or just staring at the screen intensely...",
     "Gaslighting peaks into aligning properly...",
     "Performing an exorcism on the previous plots...",
     "Pleading with the R gods for mercy...",
     "Recalculating everything... again...",
     "Casting a circle of salt around the data...",
     "Begging...",
     "Hoping that everything is still there at the end",
     "I feel like we are suposed to be loading something",
     "Integrating something completly random...",
     "Woah! Look at it go",
     "Our premium plan is faster...",
     "Trying not to crash, burn, and die",
     "Avoiding a crashout...",
     "Loading sounds...",
     "Discovering new ways to make you wait",
     "Try holding your breath through this"
   )
   
   #Modal for other button that creates modal with three other buttons
   observeEvent(input$other_options, {
     showModal(modalDialog(
       title = "Additional Options",
       h5("Select an action:"),
       br(),
       fluidRow(
         column(12,
                actionButton("regenerateplots", "Regenerate Changed Plots"),
                actionButton("rebuild_peak_areas_modal", "Rebuild Peak Area CSV"),
                actionButton("rebuild_mig_times_modal", "Rebuild Peak Migration Time CSV")
         )
       ),
       footer = modalButton("Close"),
       easyClose = TRUE
     ))
   })

  #Action button for regenerating plots
  observeEvent(input$regenerateplots, {
    
    #remove modal
    removeModal()
    
    #Initialize waiter loading screen
    initial_message <- sample(loading_messages, 1)
    waiter_show(
      html = tagList(
        spin_flower(),
        h4(id = "loading-message", initial_message, style = "color:black;")
      ),
      color = "rgba(210,210,210,0.7)"
    )
    
    #Java script for switching between different loading messages
    js_script <- sprintf("
    setTimeout(function() {
      var messages = %s;
      function updateMessage() {
        var elem = document.getElementById('loading-message');
        if (elem) {
          var newMessage = messages[Math.floor(Math.random() * messages.length)];
          elem.innerText = newMessage;
        }
      }
      setInterval(updateMessage, 4000); // Update every 5 seconds
    }, 2000);
  ", jsonlite::toJSON(loading_messages, auto_unbox = TRUE))
    runjs(js_script)
    
    #Apply plotting and area functions
    regenerate_metabolite_peak_areas()
    print("Finished applying peak area function")
    
    #Check if required data is populated before running regenerate_ggplot
    if (nrow(massData()) == 0 || nrow(refMassListData()) == 0 || nrow(parametersData()) == 0) {
      showNotification("Upload Mass List and Parameters Information!", type = "error")
    } else {
      regenerate_ggplot()
      showNotification("Successfully regenerated plots and updated metadata", type = "message")
    }
    
    #Reset modified peaks table 
    modified_peak_plots$names <- character(0)
    waiter_hide()
  })
  
  #Table to list modified plots that will be regenerated
  output$modified_peak_plots_table <- DT::renderDataTable({
    data.frame(Modified_Plots = modified_peak_plots$names)
  }, selection = 'single')
  
  #Rebuild peak area csv
  observeEvent(input$rebuild_peak_areas_modal, {
    removeModal()
    tryCatch({
      regenerate_metabolite_peak_areas()
      showNotification("Metabolite Peak Areas CSV rebuilt successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error rebuilding peak areas:", e$message), type = "error")
    })
  })
  
  #Rebuild migration CSV
  observeEvent(input$rebuild_mig_times_modal, {
    removeModal()
    tryCatch({
      regenerate_metabolite_migration_times()
      showNotification("Metabolite Migration Times CSV rebuilt successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error rebuilding migration times:", e$message), type = "error")
    })
  })
  
  #####4.7 Adjusting integration and adding peaks#####
  ######4.7.1 Function for manual peak number input######
  peak_input_modal <- function() {
    showModal(modalDialog(
      title = "Manual Peak Input",
      textInput("manual_peak_number", "Enter peak number of peak you are modifying:", ""),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_peak_number", "Confirm")
      )
    ))
  }

  
  ######4.7.2 Integration function######
  manual_peak_integration <- function(peak_number){
    
    #define variables required for plotting 
    selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
    base_file_name <- sub("\\.mz5$", "", input$file_selector)
    plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
    
    #load plot data
    plot_data <- plot_list_data_values$plot_list[[plotname]]
    eie_df <- plot_data$eie_data
    peaks_df <- plot_list_data_values$peaks_df
    peak_area_df <- plot_list_data_values$peak_area_df
    
    #Get the start and end positions from the reactive line_positions
    red_line <- line_positions$red #start
    blue_line <- line_positions$blue #end
    
    #Convert from minutes to seconds
    red_line <- as.numeric(red_line) * 60
    blue_line <- as.numeric(blue_line) *60
    
    #Filter the data between the red and blue lines. This defines the integration boundaries 
    integrated_data <- eie_df[eie_df$mt.seconds >= red_line & eie_df$mt.seconds <= blue_line, ]
    
    #Define name that will macth to existing metadata.
    intensity_column <- paste(plotname, "intensity", sep = " ")
    print("starting integration")
    
    #Load variables used for integrating 
    time <- integrated_data$mt.seconds
    intensity <- integrated_data[[intensity_column]]
     
     #Find the peak apex within the range
     apex_index <- which.max(intensity)
     apex_time <- time[apex_index]
     apex_intensity <- intensity[apex_index]
     
     #Calculate the area using the trapezoidal rule
     area <- tryCatch({
       AUC(time,
           intensity,
           method = "trapezoid",
           from = min(time),
           to = max(time),
           absolutearea = FALSE,
           na.rm = FALSE)
     }, error = function(e) {
       showNotification(paste("Error in AUC calculation:", e), type = "error")
       return(NA)
     })
     
     #Check is error is NA or not
     if (is.na(area)) {
       showNotification("Error: Area calculation returned NA", type = "error")
       return(NA)
     }
     
     #Adjust for baseline
     baseline_adjustment <- (blue_line - red_line) * min(intensity)
     area <- area - baseline_adjustment
     
     #Determine the peak number to use
     selected_peak_index <- if (!is.null(peak_number)) {
       as.numeric(peak_number)
     } else {
       input$peak_info_table_rows_selected
     }
     
     #Update the specific row in peaks_df with the new integration data
     peaks_df[selected_peak_index, paste0(plotname, ".start.seconds")] <- red_line
     peaks_df[selected_peak_index, paste0(plotname, ".apex.seconds")] <- apex_time
     peaks_df[selected_peak_index, paste0(plotname, ".end.seconds")] <- blue_line
     peaks_df[selected_peak_index, paste0(plotname, ".start_intensity")] <- min(intensity)
     peaks_df[selected_peak_index, paste0(plotname, ".apex_intensity")] <- apex_intensity
     peaks_df[selected_peak_index, paste0(plotname, ".end_intensity")] <- min(intensity)
     peaks_df[selected_peak_index, paste0(plotname, ".peak.area")] <- area
     
     #Update comment_df to remove any comments for the metabolite
     comment_df <- plot_list_data_values$comment_df
     if (plotname %in% colnames(comment_df)) {
       comment_df[selected_peak_index, plotname] <- ""}
     
     #Update peak_area_df with the new peak area
     peak_area_df[selected_peak_index, paste0(plotname)] <- area
     
     #Create temporary integration_data for the manually integrated peak to add to existing integration-data later
     temp_integration_data <- data.frame(
       "peak.number" = selected_peak_index,
       "mt.seconds" = time,
       "intensity" = intensity,
       "baseline" = min(intensity)
     )
     
     #Load existing integration_data. Convert the peak number from factor to character 
     integration_data <- plot_data$integration_data
       integration_data <- as.data.frame(lapply(integration_data, function(x) if(is.factor(x)) as.character(x) else x))
     
     #Check if the peak number exists in the existing integration_data and replace existing data if peak is already there
     if (selected_peak_index %in% integration_data$peak.number) {
       integration_data <- integration_data[integration_data$peak.number != selected_peak_index, ]
     }
     
     #Add the temporary integration_data to the existing integration_data
     integration_data <- rbind(integration_data, temp_integration_data)
     
     #Update annotation_data with the new peak.apex.seconds and peak.height.counts
     annotation_data <- plot_data$annotation_data
     annotation_data$peak.apex.seconds[selected_peak_index] <- apex_time
     annotation_data$peak.height.counts[selected_peak_index] <- apex_intensity
     annotation_data$comment[selected_peak_index] <- ""
     
     #Update values with new data
     plot_list_data_values$plot_list[[plotname]]$integration_data <- integration_data
     plot_list_data_values$plot_list[[plotname]]$annotation_data <- annotation_data
     plot_list_data_values$comment_df <- comment_df
     plot_list_data_values$peak_area_df <- peak_area_df
     plot_list_data_values$peaks_df <- peaks_df
     
     #Mark the plot as modified
     modified_peak_plots$names <- unique(c(modified_peak_plots$names, selected_plot_name))
     
     #Save the updated peaks_df back to the plot_list_data.RData file
     save_plot_data(plotname)
     
     #Clear workspace
     rm(list = c())
     gc()
     
     showNotification("Integration Successful. Metadata updated", type = "message")
  }
  
  
  ######4.7.3 Applying Integration Functions######
  observeEvent(input$manual_integrate, {
    #Check to see if the user has selected a row in the peak info table. If not, launch modal function.
    if (is.null(input$peak_info_table_rows_selected) || length(input$peak_info_table_rows_selected) == 0) {peak_input_modal()} else {
      manual_peak_integration(input$peak_info_table_rows_selected)
      regenerate_metabolite_peak_areas()
    }
    
    #Update the UI to reflect new data
    output$peak_info_table <- DT::renderDataTable({
      req(input$plot_table_rows_selected)
      selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
      base_file_name <- sub("\\.mz5$", "", run_metadata$file_name)
      plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
      plot_list <- run_metadata$plot_list
      
      if (plotname %in% names(plot_list)) {
        annotation_data <- plot_list[[plotname]]$annotation_data
        annotation_data
      } else {
        data.frame()  
      }
    }, options = list(pageLength = 13), selection = 'multiple')
  })
  
  #Button for if user uses modal peak selection
  observeEvent(input$confirm_peak_number, {
    removeModal()
    manual_peak_integration(input$manual_peak_number)
    regenerate_metabolite_peak_areas()
    
    #Update the UI to reflect new data
    output$peak_info_table <- DT::renderDataTable({
      req(input$plot_table_rows_selected)
      selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
      base_file_name <- sub("\\.mz5$", "", run_metadata$file_name)
      plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
      plot_list <- run_metadata$plot_list
      
      if (plotname %in% names(plot_list)) {
        annotation_data <- plot_list[[plotname]]$annotation_data
        annotation_data
      } else {
        data.frame()  
      }
    }, options = list(pageLength = 13), selection = 'multiple')
  })
  
  ######4.7.4 Logic for Adjusting Individual Integration baselines######
  #Reactive value for turning the baseline adjustment mode on or off for individual adjustment.
  baseline_adjustment_mode <- reactiveVal(FALSE)
  plotbaseline_adjustment_mode <- reactiveVal(FALSE)
  
  #Modal Function for specifying which peak you are editing if user fails to select a peak.
  peak_input_modal_adjustment <- function() {
    showModal(modalDialog(
      title = "Manual Peak Input",
      textInput("baseline_peak_number", "Enter peak number of peak you are adjusting:", ""),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("baseline_peak_number_adjustment", "Confirm")
      )))}
  
  #Function for adjusting individual baseline
  adjust_individual_baseline <- function(peak_index = NULL){
    
    #define variables required for plotting 
    selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
    base_file_name <- sub("\\.mz5$", "", input$file_selector)
    plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
    
    #Mark the plot as modified
    modified_peak_plots$names <- unique(c(modified_peak_plots$names, selected_plot_name))
    
    #Use clicked Y value as the new baseline
    new_baseline <- session$userData$clicked_position$y
    print(paste("New baseline:", new_baseline))
    
    #load plot data
    plot_data <- plot_list_data_values$plot_list[[plotname]]
    eie_df <- plot_data$eie_data
    peaks_df <- plot_list_data_values$peaks_df
    peak_area_df <- plot_list_data_values$peak_area_df
    
    #Determine the peak from the selected row in the peak info table or from the modal input
    if (is.null(peak_index)) {
      peak_index <- input$peak_info_table_rows_selected
    } else {
      peak_index <- as.numeric(peak_index)
    }
    
    #Get the integration boundaries from peaks_df
    left_boundary <- peaks_df[peak_index, paste0(plotname, ".start.seconds")]
    right_boundary <- peaks_df[peak_index, paste0(plotname, ".end.seconds")]
    
    #Filter the data between the left and right boundaries
    integrated_data <- eie_df[eie_df$mt.seconds >= left_boundary & eie_df$mt.seconds <= right_boundary, ]
    
    #Define name that will match to existing metadata
    intensity_column <- paste(plotname, "intensity", sep = " ")
    
    #Load variables used for integrating
    time <- integrated_data$mt.seconds
    intensity <- integrated_data[[intensity_column]]
    
    #Calculate the area using the trapezoidal rule with the new baseline
    area <- tryCatch({
      AUC(time,
          intensity,
          method = "trapezoid",
          from = min(time),
          to = max(time),
          absolutearea = FALSE,
          na.rm = FALSE)
    }, error = function(e) {
      showNotification(paste("Error in AUC calculation:", e), type = "error")
      return(NA)
    })
    
    #Check if error is NA or not
    if (is.na(area)) {
      showNotification("Error: Area calculation returned NA", type = "error")
      return(NA)
    }
    
    #Adjust for baseline
    baseline_adjustment <- (right_boundary - left_boundary) * new_baseline
    area <- area - baseline_adjustment
    
    #Update the specific row in peaks_df with the new integration data
    peaks_df[peak_index, paste0(plotname, ".start_intensity")] <- min(intensity)
    peaks_df[peak_index, paste0(plotname, ".end_intensity")] <- min(intensity)
    peaks_df[peak_index, paste0(plotname, ".peak.area")] <- area
    
    #Update peak_area_df with the new peak area
    peak_area_df[peak_index, paste0(plotname)] <- area
    
    #Create temporary integration_data for the manually integrated peak to add to existing integration-data later
    temp_integration_data <- data.frame(
      "peak.number" = peak_index,
      "mt.seconds" = time,
      "intensity" = intensity,
      "baseline" = new_baseline)
    
    #Load existing integration_data. Convert the peak number from factor to character 
    integration_data <- plot_data$integration_data
    
    integration_data <- as.data.frame(lapply(integration_data, function(x) if(is.factor(x)) as.character(x) else x))
    
    #Check if the peak number exists in the existing integration_data and replace existing data if peak is already there
    if (peak_index %in% integration_data$peak.number) {
      integration_data <- integration_data[integration_data$peak.number != peak_index, ]
    }
    
    #Add the temporary integration_data to the existing integration_data with the new baseline
    integration_data <- rbind(integration_data, temp_integration_data)
    
    #Update values with new data
    plot_list_data_values$plot_list[[plotname]]$integration_data <- integration_data
    plot_list_data_values$peak_area_df <- peak_area_df
    plot_list_data_values$peaks_df <- peaks_df
    
    
    #Save the updated peaks_df back to the plot_list_data.RData file
    save_plot_data(plotname)
    
    #Clear workspace
    rm(list = c())
    gc()
    
    showNotification("Baseline adjusted for individual peak. Area Updated. Metadata updated", type = "message")
  }
  
  ######4.7.5 Applying manual basline adjustment on individual peaks######
  #Toggle baseline adjustment mode for individual peaks
  observeEvent(input$adjust_indiv_baseline, {
    baseline_adjustment_mode(TRUE)
    plotbaseline_adjustment_mode(FALSE)
    showNotification("Baseline adjustment mode activated", type = "message")
    showNotification("Click on the plot to define the baseline for the individual peak", type = "message")
  })
  
  #Observe click events on the plot only when baseline adjustment mode is active
  observeEvent(event_data("plotly_click"), {
    if (baseline_adjustment_mode()) {
      click_data <- event_data("plotly_click")
      if (!is.null(click_data)) {
        clicked_x <- click_data$x
        clicked_y <- click_data$y
        
        print(paste("Clicked at:", clicked_x, clicked_y))
        #Store the clicked position for baseline adjustment
        session$userData$clicked_position <- list(x = clicked_x, y = clicked_y)
        
        #Check if a peak is selected
        if (is.null(input$peak_info_table_rows_selected) || length(input$peak_info_table_rows_selected) == 0) {peak_input_modal_adjustment()} else {adjust_individual_baseline(input$peak_info_table_rows_selected)
        }
        
        #Deactivate baseline adjustment mode after handling the click
        baseline_adjustment_mode(FALSE)
      } else {
        showNotification("No valid click data received.", type = "error")
        showNotification("Error: No click detected on the plot.", type = "error")
      }
    } else {
      showNotification("Baseline adjustment mode is inactive.", type = "message")  
    }
  })
  
  #Button for if user uses modal peak selection for baseline adjustment
  observeEvent(input$baseline_peak_number_adjustment, {
    removeModal()
    adjust_individual_baseline(input$baseline_peak_number)
  })
  
  
  ######4.7.6 Logic for Adjusting baseline across all peaks in selcetd plot######
  #Modal to render plot to adjust baseline on
  show_adjust_all_baselines_modal <- function() {
    showModal(modalDialog(
      title = "Adjust All Baselines",
      plotlyOutput("adjust_baseline_plot"),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_adjust_all_baselines", "Confirm")
      )
    ))
    plotbaseline_adjustment_mode(TRUE)
  }
  
  #Function to adjust all baselines in the plot
  adjust_plot_baseline <- function(plot_baseline) {
    
    #Initialize variables and load data
    selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
    base_file_name <- sub("\\.mz5$", "", input$file_selector)
    plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
    
    #Mark the plot as modified
    modified_peak_plots$names <- unique(c(modified_peak_plots$names, selected_plot_name))
    
    plot_data <- plot_list_data_values$plot_list[[plotname]]
    eie_df <- plot_data$eie_data
    peaks_df <- plot_list_data_values$peaks_df
    peak_area_df <- plot_list_data_values$peak_area_df
    integration_data <- plot_data$integration_data
    
    #Loop through all peaks 
    for (peak_index in 1:nrow(peaks_df)) {
      
      #Get boundaries of the individual peaks to integrate between 
      left_boundary <- peaks_df[peak_index, paste0(plotname, ".start.seconds")]
      right_boundary <- peaks_df[peak_index, paste0(plotname, ".end.seconds")]
      
      integrated_data <- eie_df[eie_df$mt.seconds >= left_boundary & eie_df$mt.seconds <= right_boundary, ]
      intensity_column <- paste(plotname, "intensity", sep = " ")
      time <- integrated_data$mt.seconds
      intensity <- integrated_data[[intensity_column]]
      
      #Area function
      area <- tryCatch({
        AUC(time,
            intensity,
            method = "trapezoid",
            from = min(time),
            to = max(time),
            absolutearea = FALSE,
            na.rm = FALSE)
      }, error = function(e) {
        showNotification(paste("Error in AUC calculation:", e), type = "error")
        return(NA)})
      
      #Check if area is NA
      if (is.na(area)) {
        showNotification("Error: Area calculation returned NA", type = "error")
        next
      }
      
      #Adjust baseline and update area
      baseline_adjustment <- (right_boundary - left_boundary) * plot_baseline
      area <- area - baseline_adjustment

      peaks_df[peak_index, paste0(plotname, ".peak.area")] <- area
      peak_area_df[peak_index, paste0(plotname)] <- area
      
      #Update integration data
      temp_integration_data <- data.frame(
        "peak.number" = peak_index,
        "mt.seconds" = time,
        "intensity" = intensity,
        "baseline" = plot_baseline)
      
      integration_data <- as.data.frame(lapply(integration_data, function(x) if(is.factor(x)) as.character(x) else x))
      if (peak_index %in% integration_data$peak.number) {
        integration_data <- integration_data[integration_data$peak.number != peak_index, ]
      }
      integration_data <- rbind(integration_data, temp_integration_data)
    }
    
    #Update metadata 
    plot_list_data_values$plot_list[[plotname]]$integration_data <- integration_data
    plot_list_data_values$peak_area_df <- peak_area_df
    plot_list_data_values$peaks_df <- peaks_df
    
    save_plot_data(plotname)
    
    showNotification("All baselines adjusted. Areas and metadata updated.", type = "message")
  }
  
  ######4.7.7 Applying baseline correction across all peaks in selected plot######
  #Button for the "adjust_all_baselines" button
  observeEvent(input$adjust_all_baselines, {
    plotbaseline_adjustment_mode(TRUE)
    baseline_adjustment_mode(FALSE) 
    show_adjust_all_baselines_modal()
    showNotification("Baseline adjustment mode for all peaks activated", type = "message")
    showNotification("Click on the plot to define the baseline for all peaks", type = "message")
  })
  
  #Function for rendering the selected plot in the modal
  output$adjust_baseline_plot <- renderPlotly({
    req(input$plot_table_rows_selected, db_con())
    con       <- db_con()
    folder    <- main_folder()
    file_name <- sub("\\.mz5$", "", input$file_selector)
    selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
    metabolite <- sub(paste0("^", file_name, "_"), "", selected_plot_name)
    
    plot_data <- tryCatch(
      read_plot_data_from_db(con, folder, file_name, metabolite),
      error = function(e) NULL
    )
    req(!is.null(plot_data))
    
    #generate plot
    plot <- plot_function_plotly(
      eie_data         = plot_data$eie_data,
      annotation_data  = plot_data$annotation_data,
      integration_data = plot_data$integration_data,
      label_data       = plot_data$label_data,
      x_axis_data      = plot_data$x_axis_data,
      y_axis_data      = plot_data$y_axis_data
    )
    
    plot %>%
      layout(dragmode = "select") %>%
      event_register("plotly_click")
  })
  
  #Observer for handling click data events on the plotly plot rendered in the modal to acquire the new baseline
  observeEvent(event_data("plotly_click"), {
    if (plotbaseline_adjustment_mode()) {
      click_data <- event_data("plotly_click")
      if (!is.null(click_data)) {
        plot_baseline <- click_data$y
        
        #Update all baselines with the new value and apply the baseline fun ction
        adjust_plot_baseline(plot_baseline)
        plotbaseline_adjustment_mode(FALSE)
        
        #Performance notifications 
        showNotification("New baseline set. All baselines adjusted.", type = "message")
      } else {
        showNotification("No valid click data received.", type = "error")
        showNotification("Error: No click detected on the plot.", type = "error")}
    } else {
      showNotification("Baseline adjustment mode is inactive.", type = "message")}
  })
  
  #Observer for the "confirm_adjust_all_baselines" button
  observeEvent(input$confirm_adjust_all_baselines, {
    removeModal()
  })
  
  #####4.8 Server logic for dual EIE viewing.#####
  #Populate plot tables of plots you wish you display
  output$plottable1 <- renderDataTable({
    data.frame(Plot = filtered_plot_names())
  }, selection = 'single')
  
  output$plottable2 <- renderDataTable({
    data.frame(Plot = filtered_plot_names())
  }, selection = 'single')
  
  output$plottable3 <- renderDataTable({
    data.frame(Plot = filtered_plot_names())
  }, selection = 'single')

  
  #Observe selected plots and render a combined subplot of these functions to zoom in on peaks
    observe({
      sel1 <- input$plottable1_rows_selected
      sel2 <- input$plottable2_rows_selected
      sel3 <- input$plottable3_rows_selected
      
      selected_indices <- unique(na.omit(c(sel1, sel2, sel3)))
      
      #generate combined plot depending on if 2 OR 3 plots are selected.
      if (length(selected_indices) >= 2) {
        req(db_con())
        con            <- db_con()
        folder         <- main_folder()
        file_name      <- sub("\\.mz5$", "", input$file_selector)
        selected_names <- filtered_plot_names()[selected_indices]
        
        #Generate each plot from SQL database on demand
        plots <- lapply(selected_names, function(name) {
          metabolite <- sub(paste0("^", file_name, "_"), "", name)
          plot_data  <- tryCatch(
            read_plot_data_from_db(con, folder, file_name, metabolite),
            error = function(e) NULL
          )
          if (is.null(plot_data)) return(NULL)
          plot_function_plotly(
            eie_data         = plot_data$eie_data,
            annotation_data  = plot_data$annotation_data,
            integration_data = plot_data$integration_data,
            label_data       = plot_data$label_data,
            x_axis_data      = plot_data$x_axis_data,
            y_axis_data      = plot_data$y_axis_data
          )
        })
        
        #Remove any plots that failed to load
        plots <- Filter(Negate(is.null), plots)
        
        #render if successful and combine plot names into one title 
        if (length(plots) == length(selected_names)) {
          cleaned_names  <- sub(paste0("^", file_name, "_"), "", selected_names)
          combined_title <- paste("EIE of", paste(cleaned_names, collapse = ", "))
          
          #combine individual plots into stacked subplots with shared x axis
          output$combined_plot <- renderPlotly({
            subplot(plots, nrows = length(plots), 
                    shareX = TRUE, titleX = TRUE, titleY = TRUE) %>%
              layout(title = combined_title)
          })
        }
      }
    })    
  
  #####4.9 Marking Interferences#####
  observeEvent(input$interferences, {
    req(input$peak_info_table_rows_selected)
    
    #Store the current state before making changes
    previous_metadata$comment_df <- run_metadata$comment_df
    previous_metadata$plot_list <- run_metadata$plot_list
    previous_metadata$peaks_df <- run_metadata$peaks_df
    previous_metadata$peak_mt_df <- run_metadata$peak_mt_df
    previous_metadata$peak_area_df <- run_metadata$peak_area_df
    previous_metadata$mi_df <- run_metadata$mi_df
    
    selected_rows <- input$peak_info_table_rows_selected
    
    #Extract the current annotation_data
    selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
    base_file_name <- sub("\\.mz5$", "", input$file_selector)
    plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
    plot_list <- plot_list_data_values$plot_list
    
    #Mark the plot as edited so only plots that have been edited can be regenerated later
    modified_peak_plots$names <- unique(c(modified_peak_plots$names, selected_plot_name))
    
    if (plotname %in% names(plot_list)) {
      annotation_data <- plot_list[[plotname]]$annotation_data
      integration_data <- plot_list[[plotname]]$integration_data
      
      #Loop through each selected row
      for (selected_row in selected_rows) {
        #Get the peak number from the selected row to allow matching accross various acccesible data streams
        peak_number <- annotation_data$peak.number[selected_row]
        
        #Update the comment column for the selected peak
        annotation_data$comment[selected_row] <- "Interferred"
        
        #Update comment_df so peak areas can be recauculated/reintegrated later on
        if (plotname %in% colnames(plot_list_data_values$comment_df)) {
          plot_list_data_values$comment_df[peak_number, plotname] <- "Interferred"
        }
        
        #Remove all integration_data associated with the relevant peak number so regenerated plots will update with the interferred peaks once regenerated
        integration_data <- integration_data[integration_data$peak.number != peak_number, ]
      }
      
      #Update the plot_list with the modified annotation_data
      plot_list[[plotname]]$annotation_data <- annotation_data
      plot_list[[plotname]]$integration_data <- integration_data
      
      #Update the reactive plot_list_data
      plot_list_data_values$plot_list <- plot_list
    }

    save_plot_data(plotname)
    
    #Reload the updated data
    subfolder_path <- file.path(input$results_folder, "Data", run_metadata$file_name)
    reload_plot_data(subfolder_path)
    
    #Update reactive values directly
    plot_list_data_values$plot_list <- plot_list
    plot_list_data_values$comment_df <- plot_list_data_values$comment_df
    
    
    #Update the peak_info_table
    output$peak_info_table <- DT::renderDataTable({
      req(input$plot_table_rows_selected)
      selected_plot_name <- filtered_plot_names()[input$plot_table_rows_selected]
      
      base_file_name <- sub("\\.mz5$", "", input$file_selector)
      plotname <- sub(paste0("^", base_file_name, "_"), "", selected_plot_name)
      plot_list <- plot_list_data_values$plot_list
      
      #Ensure the plotname exists in the plot_list
      if (plotname %in% names(plot_list)) {
        annotation_data <- plot_list[[plotname]]$annotation_data
        annotation_data
      } else {
        data.frame()
      }
    }, options = list(pageLength = 13), selection = 'multiple')
  })
  
####5. Downstream Processing ####   
  #####5.1 Initializing Environment for loading peak area and migration time data#####
  #Define the function to load files
  grab_metabolite_files <- function(results_folder2) {
    peak_areas_path <- file.path(results_folder2, "Metabolite Peak Areas.csv")
    migration_times_path <- file.path(results_folder2, "Metabolite Migration Times.csv")
    
    peak_areas <- if (file.exists(peak_areas_path)) {
      read.csv(peak_areas_path, stringsAsFactors = FALSE, check.names = FALSE)} else {
      warning("File not found: ", peak_areas_path)
      NULL}
    
    migration_times <- if (file.exists(migration_times_path)) {
      read.csv(migration_times_path, stringsAsFactors = FALSE, check.names = FALSE)} else {
      warning("File not found: ", migration_times_path)
      NULL}
  
    list(peak_areas = peak_areas, migration_times = migration_times)
  }
  
  #Define reactive variables for the dataframes
  peak_areas_data <- reactiveVal(NULL)
  migration_times_data <- reactiveVal(NULL)

  # A different saved analysis must be loaded before downstream transformations.
  observeEvent(input$results_folder2, {
    peak_areas_data(NULL)
    migration_times_data(NULL)
  }, ignoreNULL = FALSE, priority = 100)
  
  #Observer for the action button
  observeEvent(input$load_migration_area, {
    req(input$results_folder2)
    files <- grab_metabolite_files(input$results_folder2)
    peak_areas_data(files$peak_areas)
    migration_times_data(files$migration_times)
    
    showNotification("Data loaded", type = "message")
  })
  
  #Set results folder to work out of
  main_folder2 <- reactive({
      req(input$results_folder2)
      input$results_folder2
    })
  
  #Observe functions for selecting files and results folders to work with
  observe({
    results_folders_reactive()  # Refresh when completed results are discovered.
    files <- uploadedmz5()$FileName
    filename <- input$file_selector2
    if (length(filename) != 1L || !filename %in% files) filename <- head(files, 1L)
    folders <- list_results_folders(filename)
    current <- isolate(input$results_folder2)
    updateSelectInput(session, "results_folder2", choices = c("Select a saved analysis" = "",
                      stats::setNames(folders, folders)),
      selected = if (length(current) == 1L && current %in% folders) current else "")
  })
  observe({
    files <- uploadedmz5()$FileName
    current <- isolate(input$file_selector2)
    updateSelectInput(session, "file_selector2", choices = files,
      selected = if (length(current) == 1L && current %in% files) current else head(files, 1L))
  })
  
  
  #####5.2 Plot for m/z vs migration time #####
  #Function to extract m/z and names from the name vector
  extract_mass <- function(name_vec) {
    mz <- as.numeric(sub("_.*", "", name_vec))
    names <- sub("^[^_]*_", "", name_vec)
    data.frame(mz = mz, name = names)
  }
  
  #Reactive expression to create a data frame for plotting
  mzMT_data <- reactive({
    
    #Initialize variables 
    name_vec <- c(refMassListData()$name, massData()$name)
    req(migration_times_data(), name_vec, input$peak_number)
    mz_df <- extract_mass(name_vec)
    migration_times <- migration_times_data()
    base_file_name <- sub("\\.mz5$", "", input$file_selector2)
    
    #Filter data
    start_row <- which(migration_times$file.name == base_file_name)
    next_file_row <- which(migration_times$file.name != "" & migration_times$file.name != base_file_name & seq_along(migration_times$file.name) > start_row)[1]
    end_row <- ifelse(is.na(next_file_row), nrow(migration_times), next_file_row - 1)
    
    #Filter the data for the selected base file name block
    filtered_migration_time <- migration_times[start_row:end_row, ]
    migration_times_filtered <- filtered_migration_time[filtered_migration_time$peak.number == input$peak_number, ]
    migration_times_filtered <- migration_times_filtered[, !colnames(migration_times_filtered) %in% c("file.name", "peak.number")]
    
    #Extract m/z values from the column names in migration_times_filtered and match them
    mz_values <- as.numeric(sub("^([0-9.]+)_.*", "\\1", colnames(migration_times_filtered)))
    matched_columns <- colnames(migration_times_filtered)[match(mz_df$mz, mz_values)]
    
    #Extract names from metadata to colour code metabolites and internal standards
    ref_names <- sub("^[0-9.]+_", "", refMassListData()$name)
    
    #Create a data frame for plotting
    plot_df <- data.frame(
      mz = mz_df$mz,
      migration_time = unlist(migration_times_filtered[, matched_columns]),
      name = mz_df$name,
      Type = ifelse(mz_df$name %in% ref_names, "Internal Standard", "Metabolite")
    )
  })
  
  #Make plotly plot of m/z vs. migration time. 
  output$mz_vs_migration_plot <- renderPlotly({
    plot_df <- mzMT_data()
    
    #Make the plot
    p <- ggplot(plot_df, aes(x = migration_time, y = mz, text = name, fill = Type)) +
      geom_point(shape = 21, size = 3, alpha = 0.7, color = "black") + 
      scale_fill_manual(values = c("Metabolite" = "green", "Internal Standard" = "darkorange")) +
      labs(title = "m/z versus Migration Time (min) for all analytes", x = "Migration Time (min)", y = "m/z") +
      scale_x_continuous(limits = c(0, 40), breaks = seq(0, 50, by = 2)) +
      theme_classic() +  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, size = 16))
    
    ggplotly(p, tooltip = "text") %>% layout(legend = list(orientation = "h", x = 0.4, y = -0.2))
  })
  
  #Render the peak number selector
  output$peak_number_selector <- renderUI({
    num_peaks <- parametersData()$number.of.injections
    selectInput("peak_number", "Select Peak Number:", choices = 1:num_peaks)
  })
  
  
  #####5.3 Normalize to Creatine, F-Phe, or Cl-Tyr#####
  #Show modal for normalization options
  observeEvent(input$normalize, {
    showModal(modalDialog(
      title = "Normalization Options",
      selectInput("normalize_by", "Normalize By", choices = c("114.0667_Creatinine", "184.0774_F-Phe", "216.0427_Cl-Tyr", "182.0617_F-Phe","198.0566_F-Tyr", "214.0271_Cl-Tyr","207.0115_NMS", "212.9412_Bromohypoxanthine")),
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_normalize", "Normalize")
      ),
      easyClose = TRUE
    ))
  })
  
  #Event handler for normalizing the data.
  observeEvent(input$confirm_normalize, {
    removeModal()
    
    #Load data and initialize the environment youre working in.
    name_vec <- c(refMassListData()$name, massData()$name)
    results_folder2 <- main_folder2()
    area <- peak_areas_data()
    
    #Check if data is NULL or empty to prevent crashing 
    if (is.null(area) || nrow(area) == 0) {
      showNotification("No data uploaded. Please upload data before proceeding.", type = "error")
      return()
    }
    
    #Function for normalizing peak areas
    normalize_peak_areas <- function(data, normalize_by) {
      if (!(normalize_by %in% name_vec)) {
        showNotification("Selected normalization metabolite not found in data columns.", type = "error")
        return(data)}
      
      #Convert non-numeric values to NA and ensure all columns are now numeric
      data <- data %>%
        mutate(across(all_of(name_vec), ~ suppressWarnings(ifelse(grepl("^[0-9.]+$", .), as.numeric(.), NA_real_))))
      if (all(is.na(data[[normalize_by]]))) {
        showNotification("Normalization column does not contain any numeric values.", type = "error")
        return(data)}
      
      #Normalize each row by the corresponding area for the metabolite you are using to normalize.
      data <- data %>%
        rowwise() %>%
        mutate(across(all_of(name_vec), ~ ifelse(!is.na(.data[[normalize_by]]) & !is.na(.), . / .data[[normalize_by]], .)))
      return(data)
    }
    
    #Apply function to data
    normalized_data <- normalize_peak_areas(area, input$normalize_by)
    
    #Update reactive data
    peak_areas_data(normalized_data)
     
    #Save normalized data to a new CSV file so original datafile is preserved
    write.csv(normalized_data, file.path(results_folder2, "Normalized Metabolite Peak Areas.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    showNotification("Data successfully normalized", type = "message")
  })
  
  
  #####5.4 Replacing missing values in data#####
  observeEvent(input$missing_data, {
    
    #Intialize environment, define variables
    results_folder2 <- main_folder2()
    data <- peak_areas_data()
    name_vec <- c(refMassListData()$name, massData()$name)
    
    #Check if data is NULL or empty to prevent crashing 
    if (is.null(data) || nrow(data) == 0) {
      showNotification("No data uploaded. Please upload data before proceeding.", type = "error")
      return()
    }
    
    #Convert non-numeric values to NA
    data <- data %>%
      mutate(across(all_of(name_vec), ~ suppressWarnings(ifelse(grepl("^[0-9.]+$", .), as.numeric(.), NA_real_))))
    
    #Replace missing values with the minimum/five
    if (input$missing_data_method == "Minimum Values/5") {
      for (col in name_vec) {
        min_value <- min(data[[col]], na.rm = TRUE) / 5
        data[[col]][is.na(data[[col]])] <- min_value}
    } 
    
    #Replace missing values values with 0
    else if (input$missing_data_method == "Missing values = 0") {
      data[is.na(data)] <- 0}
    
    #Update reactive data
    peak_areas_data(data)
    
    #Save to a .CSV file
    write.csv(data, file.path(results_folder2, "Corrected Missing Values Metabolite Peak Areas.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    showNotification("Data successfully corrected", type = "message")
  })
  
  
  #####5.5 Connect Peak Position with Metadata#####
  
  #Observe action button for connecting metadata
  observeEvent(input$connect_metadata, {
    results_path <- main_folder2()
    data_folder <- file.path(results_path, "Data")
    
    #Look for path to metadata file
    metadata_path <- file.path(data_folder, "metadata.csv")
    
    if (!file.exists(metadata_path)) {
      showNotification("metadata.csv not found in Data subfolder. Are you working out of the right results folder?", type = "warning")
      return(NULL)
    }
    
    #Load data
    metadata <- read.csv(metadata_path, check.names = FALSE)
    area <- peak_areas_data()
    
    #Check if metadata is NULL or empty 
    if (is.null(metadata) || nrow(metadata) == 0) {
      showNotification("No metadata uploaded. Please upload data before proceeding.", type = "error")
      return()
    }
    
    #Check if data is NULL or empty 
    if (is.null(area) || nrow(area) == 0) {
      showNotification("No data uploaded. Please upload data before proceeding.", type = "error")
      return()
    }
    
    #identify peak columns in metadata file. Designed to by dynamic with changing peak numbers
    peak_columns <- grep("^[0-9]+$", names(metadata), value = TRUE)
    
    #Remove empty columns (just in case) and transform data to long format
    metadata <- metadata[, colSums(is.na(metadata)) < nrow(metadata)]
    metadata_long <- metadata %>%
      pivot_longer(cols = all_of(peak_columns), names_to = "peak.number", values_to = "PBM Sample ID") %>%
      select(`Data file`, peak.number, `PBM Sample ID`)
    
    #Filter metadata file to connect filenames in actual results with the samples IDs for each file & peak
    filter_data_by_file_name <- function(data, file_name) {
      start_row <- which(data$file.name == file_name)
      next_file_row <- which(data$file.name != "" & data$file.name != file_name & seq_along(data$file.name) > start_row)[1]
      end_row <- ifelse(is.na(next_file_row), nrow(data), next_file_row - 1)
      data[start_row:end_row, ]}
    
    merged_data <- do.call(rbind, lapply(unique(area$file.name[area$file.name != ""]), function(file_name) {
      filtered_area <- filter_data_by_file_name(area, file_name)
      filtered_metadata <- metadata_long[metadata_long$`Data file` == file_name, ]
      merged <- merge(filtered_area, filtered_metadata, by.x = "peak.number", by.y = "peak.number", all.x = TRUE)
      merged <- merged %>% select(-`Data file`) %>% relocate(`PBM Sample ID`, .after = file.name)
      return(merged)
    }))
    
    #Trim any leading or trailing whitespace (found in some of the metadata files for unknown resason)
    merged_data$`PBM Sample ID` <- trimws(merged_data$`PBM Sample ID`, whitespace = "[\\h\\v]")
    
    #Add merged data to reactive data file and save it to a .CSV
    peak_areas_data(merged_data)
    
    write.csv(merged_data, file.path(results_path, "merged_peak_area.csv"), row.names = FALSE, fileEncoding = "UTF-8")
    showNotification("Metadata succesffuly connected to peak areas", type = "message")
  })
  
  
  #####5.6 Control Charts#####
  #Update the selectInput choices based on the reactive value
  observe({
    ref_data <- refMassListData()
    if (!is.null(ref_data) && nrow(ref_data) > 0) {
      updateSelectInput(session, "selected_metabolite", choices = ref_data$name)
    }
  })
  
  ######5.6.1 Control Charts for Non-normalized data######
  observeEvent(input$create_control_chart, {
    req(input$selected_metabolite)
    
    #Initialize environment, load data, acquire names of selected metabolite for control chart
    selected_metabolite <- input$selected_metabolite
    results_path <- main_folder2()
    files <- grab_metabolite_files(results_path)
    
    #Load basic raw matrix. This data wont have any correction for missing values or anything normalized
    raw_peak_areas <- files$peak_areas
    
    #Extract the data for the selected metabolite
    raw_metabolite_areas <- raw_peak_areas[[selected_metabolite]]
    
    #Generate control chart
    output$control_chart_raw <- renderPlotly({
      plot_ly() %>%
        add_trace(x = seq_along(raw_metabolite_areas), y = raw_metabolite_areas, type = 'scatter', mode = 'lines+markers', name = 'Metabolite Peak Area') %>%
        add_lines(x = seq_along(raw_metabolite_areas), y = mean(raw_metabolite_areas), line = list(color = 'blue'), name = 'Mean') %>%
        add_lines(x = seq_along(raw_metabolite_areas), y = mean(raw_metabolite_areas) + 3 * sd(raw_metabolite_areas), line = list(color = 'black', dash = 'dash'), name = 'Upper Control Limit') %>%
        add_lines(x = seq_along(raw_metabolite_areas), y = mean(raw_metabolite_areas) - 3 * sd(raw_metabolite_areas), line = list(color = 'black', dash = 'dash'), name = 'Lower Control Limit') %>%
        layout(title = paste("Control Chart for", selected_metabolite, "(Raw Data)"),
               xaxis = list(title = "Sample", showgrid =FALSE, dtick = 5),
               yaxis = list(title = "Peak Area", showgrid =FALSE),
               legend = list(orientation = 'h', x = 0.5, xanchor = 'center', y = -0.2), 
               font = list(family = "sans-serif"))
    })
    showNotification("Control Chart Generated", type = "message")
  })
  
  
  ######5.6.2 Control Charts for normalized data######
  observeEvent(input$create_control_chart, {
    req(input$selected_metabolite)
    
    #Initialize environment, load data, acquire names of selected metabolite for control chart
    selected_metabolite <- input$selected_metabolite
    peak_areas <- peak_areas_data()
    metabolite_data <- peak_areas[[selected_metabolite]]
    
    #Generate control chart
    output$control_chart_processed <- renderPlotly({
      plot_ly() %>%
        add_trace(x = seq_along(metabolite_data), y = metabolite_data, type = 'scatter', mode = 'lines+markers', name = 'Metabolite Peak Area') %>%
        add_lines(x = seq_along(metabolite_data), y = mean(metabolite_data), line = list(color = 'blue'), name = 'Mean') %>%
        add_lines(x = seq_along(metabolite_data), y = mean(metabolite_data) + 3 * sd(metabolite_data), line = list(color = 'black', dash = 'dash'), name = 'Upper Control Limit') %>%
        add_lines(x = seq_along(metabolite_data), y = mean(metabolite_data) - 3 * sd(metabolite_data), line = list(color = 'black', dash = 'dash'), name = 'Lower Control Limit') %>%
        layout(title = paste("Control Chart for", selected_metabolite, "(Processed Data)"),
               xaxis = list(title = "Sample", showgrid =FALSE, dtick = 5),
               yaxis = list(title = "Peak Area", showgrid =FALSE),
               legend = list(orientation = 'h', x = 0.5, xanchor = 'center', y = -0.2),
               font = list(family = "sans-serif"))
    })
  })
  
  
  #####5.7 CV and PCA determination for QC variance checks#####
  ######5.7.1 Define functions to calculate CV######
  
  #CV calculation function 
  calculate_cv <- function(data){
    cv <- apply(data, 2, function(x) sd(x) / mean(x) * 100)
    return(cv)
  }
  
  #Reactive expression to calculate CVs
  cv_data <- reactive({
    data <- peak_areas_data()
    metabolite_data <- data[, 4:ncol(data)]
    cv <- calculate_cv(metabolite_data)
    return(cv)
  })
  
  ######5.7.2 Define functions for creating PCA plots######
  #Function to perform Pareto scaling
  pareto_scale <- function(data) {
    scaled_data <- scale(data, center = TRUE, scale = sqrt(apply(data, 2, sd)))
    return(scaled_data)
  }
  
  #Function for removing columns with zero variance
  remove_zero_variance <- function(data) {
    data_matrix <- as.matrix(data)
    zero_variance_cols <- apply(data_matrix, 2, function(x) var(x) == 0)
    data_matrix <- data_matrix[, !zero_variance_cols]
    return(data_matrix)
  }
  
  #Function to perform PCA with additional checks
  perform_pca <- function(data) {
    data <- remove_zero_variance(data)
    
    #Check for missing or infinite values
    if (any(is.na(data)) || any(is.infinite(data))) {
      showNotification("PCA Data contains missing or infinite values.", type = "error")
    }
    #Identify columns with all missing or infinite values
    cols_to_exclude <- which(apply(data, 2, function(x) all(is.na(x)) || all(is.infinite(x))))
    if (length(cols_to_exclude) > 0) {
      excluded_cols <- colnames(data)[cols_to_exclude]
      showNotification(paste("Excluding columns with all missing or infinite values:", paste(excluded_cols, collapse = ", ")), type = "warning")
      data <- data[, -cols_to_exclude]}
    
    #Check for any remaining missing or infinite values
    if (any(is.na(data)) || any(is.infinite(data))) {
      showNotification("Data still contains missing or infinite values. Filter failed. Please manually edit data.", type = "error")
      return(NULL)
    } else {
      showNotification("Filtering successful.", type = "message")
    }

    scaled_data <- pareto_scale(data)
    pca_result <- PCA(scaled_data, graph = FALSE)
    return(pca_result)
  }
  
  #PCA expression
  pca_data <- reactive({
    data <- peak_areas_data()
    metabolite_data <- data[, 4:ncol(data)]
    print(head(metabolite_data))
    pca_result <- perform_pca(metabolite_data)
    return(pca_result)
  })
  
  ######5.7.3 Function for extracting QC and non-QC catagories from user input######
  #Separate QCs and other samples based on user input
  average_cv <- reactive({
    data <- peak_areas_data()
    
    #Ensure column name exists
    if (!"PBM Sample ID" %in% colnames(data)) {
      showNotification("Column 'PBM Sample ID' not found in data. Please connect metadata to results", type = "error")
      return()
    }
    
    #Extract user specified dynamic categories
    user_categories <- strsplit(input$category_input, ",")[[1]]
    user_categories <- sapply(user_categories, function(x) {
      kv <- strsplit(x, "=")[[1]]
      setNames(trimws(gsub("'", "", kv[2])), trimws(kv[1]))
    }, USE.NAMES = FALSE)
    
    #Exclude specific sample IDs from CV analysis.
    exclude_ids <- strsplit(input$exclude_ids, ",")[[1]]
    exclude_ids <- trimws(exclude_ids)
    data <- data %>% filter(!`PBM Sample ID` %in% exclude_ids)
    
    qc_samples <- data %>% filter(grepl(user_categories[["QC"]], `PBM Sample ID`, ignore.case = TRUE))
    non_qc_samples <- data %>% filter(!grepl(user_categories[["QC"]], `PBM Sample ID`, ignore.case = TRUE))
    #Deal with cases where columns contain all missing or all infinite samples. Occurs when there is ZERO detection of a metabolite in ANY sample which makes the data either all NAs or ALL infinite. Occurs rarely but will crash app if not handled correctly.
    #Identify columns with all missing or infinite values in QC and non-QC samples
    qc_exclusions <- which(apply(qc_samples[, 4:ncol(qc_samples)], 2, function(x) all(is.na(x)) || all(is.infinite(x))))
    if (length(qc_exclusions) > 0) {
      excluded_qc_cols <- colnames(qc_samples)[qc_exclusions + 3]
      showNotification(paste("Excluding columns with all missing or infinite values in QC samples:", paste(excluded_qc_cols, collapse = ", ")), type = "warning", duration = 6)
      qc_samples <- qc_samples[, -qc_exclusions - 3] 
    }
    
    non_qc_exclusions <- which(apply(non_qc_samples[, 4:ncol(non_qc_samples)], 2, function(x) all(is.na(x)) || all(is.infinite(x))))
    if (length(non_qc_exclusions) > 0) {
      excluded_non_qc_cols <- colnames(non_qc_samples)[non_qc_exclusions + 3] 
      showNotification(paste("Excluding columns with all missing or infinite values in non-QC samples:", paste(excluded_non_qc_cols, collapse = ", ")), type = "warning", duration = 6)
      non_qc_samples <- non_qc_samples[, -non_qc_exclusions - 3] 
    }
    
    #Ensure the remaining columns have the same number of rows
    common_cols <- intersect(colnames(qc_samples), colnames(non_qc_samples))
    qc_samples <- qc_samples[, common_cols]
    non_qc_samples <- non_qc_samples[, common_cols]
    
    #Check for any remaining missing or infinite values in samples
    if (any(sapply(qc_samples[, 4:ncol(qc_samples)], is.na)) || any(sapply(qc_samples[, 4:ncol(qc_samples)], is.infinite))) {
      showNotification("QC samples still contain missing or infinite values. Filter failed. Please manually edit data.", type = "error")
      return(NULL)
    }
    
    if (any(sapply(non_qc_samples[, 4:ncol(non_qc_samples)], is.na)) || any(sapply(non_qc_samples[, 4:ncol(non_qc_samples)], is.infinite))) {
      showNotification("Non-QC samples still contain missing or infinite values. Filter failed. Please manually edit data.", type = "error")
      return(NULL)
    }
    
    showNotification("Filtering successful.", type = "message")
    
    #Calculate CVs for QC and non-QC samples
    qc_cv <- calculate_cv(qc_samples[, 4:ncol(qc_samples)])
    non_qc_cv <- calculate_cv(non_qc_samples[, 4:ncol(non_qc_samples)])
    
    #Combine QC and non-QC samples for PCA and store values
    combined_data <- rbind(qc_samples, non_qc_samples)
    pca_result <- perform_pca(combined_data[, 4:ncol(combined_data)])
    list(qc_cv = mean(qc_cv), non_qc_cv = mean(non_qc_cv), non_qc_cv_values = non_qc_cv, pca_result = pca_result, qc_samples = qc_samples, non_qc_samples = non_qc_samples, qc_cv_values = qc_cv)
  })
  
  
  ######5.7.4 Output for CV processing######
  observeEvent(input$compute_cv, {
    #######5.7.4.1 CV plot#######
    
    #Check if peak_areas_data() is NULL
    if (is.null(peak_areas_data())) {
      showNotification("No data uploaded. Please upload data before proceeding.", type = "error")
      return()}
    
    avg_cv <- average_cv()
    
    median_cv <- median(avg_cv$non_qc_cv_values, na.rm = TRUE)
    median_qc_cv <- median(avg_cv$qc_cv_values, na.rm = TRUE)
    
    #Create a data frame for the average CV values
    cv_table <- data.frame(
      Type = c("Technical CV (%)", "Biological CV (%)"),
      `Average CV` = c(avg_cv$qc_cv, avg_cv$non_qc_cv))
    
    #Render the table
    output$average_cv_table <- renderTable({
      cv_table})
    
    #Extract m/z values and names and create a dataframe to plot m/z vs CV
    mz_names <- extract_mass(colnames(peak_areas_data())[4:ncol(peak_areas_data())])
    
    #Filter out excluded columns that contained NAs or Infs
    valid_indices <- which(!is.na(avg_cv$non_qc_cv_values) & !is.infinite(avg_cv$non_qc_cv_values))
    filtered_mz_names <- mz_names$mz[valid_indices]
    filtered_colnames <- colnames(peak_areas_data())[4:ncol(peak_areas_data())][valid_indices]
    filtered_cv_values <- avg_cv$non_qc_cv_values[valid_indices]
    
    #Finalize dataframe
    plot_data <- data.frame(mz = filtered_mz_names, cv = filtered_cv_values, name = filtered_colnames)
    
    #Sort the data frame by m/z values
    plot_data <- plot_data[order(plot_data$mz), ]
    
    #Create the plot
    output$cv_plot <- renderPlotly({
      plot_ly(plot_data, x = ~mz, y = ~cv, type = 'scatter', mode = 'lines+markers', text = ~name, hoverinfo = 'text+x+y', name = 'Non-QC CV') %>%
        add_trace(x = ~mz, y = rep(median_cv, length(plot_data$mz)), type = 'scatter', mode = 'lines', line = list(dash = 'dash', color = 'black'), name = 'Median Non-QC CV') %>%
        add_trace(x = ~mz, y = avg_cv$qc_cv_values, type = 'scatter', mode = 'lines+markers', marker = list(color = '#780116'), line = list(color = '#780116'), name = 'QC CV') %>%
        add_trace(x = ~mz, y = rep(median_qc_cv, length(plot_data$mz)), type = 'scatter', mode = 'lines', line = list(dash = 'dash', color = 'green'), name = 'Median QC CV') %>%
        layout(title = "CV per Metabolite (QC and non-QC samples)",
               xaxis = list(title = "m/z"),
               yaxis = list(title = "CV (%)"),
               legend = list(orientation = 'h', x = 0.5, xanchor = 'center', y = -0.2))
    })
    
    #######5.7.4.2 Create PCA plot#######
    output$pca_plot <- renderPlot({
      
      #Intialize environment 
      qc_samples <- avg_cv$qc_samples
      non_qc_samples <- avg_cv$non_qc_samples
      pca_data <- data.frame(avg_cv$pca_result$ind$coord)
      pca_data$SampleType <- rep(c("QC", "Non-QC"), c(nrow(qc_samples), nrow(non_qc_samples)))
      variance_explained <- avg_cv$pca_result$eig[1:2, 2]
      results_path <- main_folder2()
      data_folder <- file.path(results_path, "Data")
      
      #Generate plot 
      pca_plot <- ggplot(pca_data, aes(x = Dim.1, y = Dim.2, color = SampleType)) +
        geom_point() +
        stat_ellipse() +
        scale_color_manual(values = c("QC" = "#38A700", "Non-QC" = "#780116")) +
        theme_classic() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
        ) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
        labs(
          title = "Quality Control PCA Plot for QCs and Non-QC Samples (PC1 & PC2)",
          x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
          y = paste0("PC2 (", round(variance_explained[2], 2), "%)"), colour = "Sample Type"
        )
      
      #Save the plot
      ggsave(file.path(data_folder, "QC_PCA_plot.png"), plot = pca_plot, width = 9, height = 6, dpi = 300)
      
      #Display the plot
      pca_plot
      
    })
  })
  
  
  #####5.8 Display Data#####
  #Display data
  output$peak_area_table <- renderDT({
    datatable(peak_areas_data(), extensions = 'FixedColumns', options = list(
      scrollX = TRUE,
      scrollY = "400px",
      scrollCollapse = TRUE,
      paging = FALSE,
      pageLength = 20,
      lengthMenu = c(20, 50, 100, "All"),
      fixedColumns = list(leftColumns = 4)  
    ))
  })
  
  #####5.9 Integrated Batch Correction#####
  ######5.9.1 Initialize Environment######
  observeEvent(input$load_batchcorr_packages,{
    #Note: I am aware that installing packages when the ShinyApp is already running is not ideal. However, for some reason, to install the package RcppArmadillo, which is a dependency of ChemometricsWithR and BatchCorrMetabolomics which we need to run the Batch Correction, R requires a C++ compiler of C++14. The only way I was able to achieve this was by finding the file Makevars.win located where R is stored on your computer and adding the line: CXX_STD=CXX14. Since the app crashes without this, I didnt want to include this as a requirement on setup... Also, BatchCorr wont be done with every run so why load it before?
    
    showNotification("Don't forget to type Yes/No or Y/N in R console to install packages", type = "warning")
    
    #Install packages
    if (!requireNamespace("RcppArmadillo", quietly = TRUE)) {
      install.packages("RcppArmadillo")}
    library(RcppArmadillo)
    
    #Install ChemometricWithR. If you just try to have R download this when installing BatchCorrMetabolomics, it fails. Hence the separation
    if (!requireNamespace("ChemometricsWithR", quietly = TRUE)) {
      remotes::install_github("rwehrens/ChemometricsWithR")}
    library(ChemometricsWithR)
    
    #Install BatchCorrMetabolomics after installing troublesome and required dependencies 
    if (!requireNamespace("BatchCorrMetabolomics", quietly = TRUE)) {
      remotes::install_github("rwehrens/BatchCorrMetabolomics")}
    library(BatchCorrMetabolomics)
    
    showNotification("Environemnt Initialized", type = "message")
  })
  
  #Load Batch Corr Order
  batch_order <- reactive({
    req(input$batch_corr_order)  #Ensure file is uploaded
    read.csv(input$batch_corr_order$datapath, stringsAsFactors = FALSE)
  })
  
  
  ######5.9.2 Perform Batch Correction######
  observeEvent(input$batch_corr,{
  
    #Load data.
    Rawdata <- peak_areas_data()
    batch_metadata <- batch_order()
  
    #Perform log10 transform on raw data before correction
    data_BatchCorr <- Rawdata
    
    #Define which columns are your metabolites
    metCols <- 4:ncol(data_BatchCorr)
    
    #Apply log10 in place
    data_BatchCorr[ , metCols] <- log10(data_BatchCorr[ , metCols])
    
    #Create variables necessary for batch correction with initial conditions
    minBatchOccurrence.Ave <- 2
    minBatchOccurrence.Line <- 2
    conditions <- c("")           #Replace missing variables with not available
    experiments <- c(t(outer(c("Q", "S"), conditions, paste, sep = "")))
  
    methods <- rep("lm", length(experiments))
    methods[grep("c", experiments)] <- "tobit"
    LOD <- as.numeric(as.character(min(data_BatchCorr[ , metCols][!is.na(data_BatchCorr[ , metCols])])))
    imputeValues  <- rep(NA, length(experiments))
    
    refSamples <- list(
      "Q" = which(batch_metadata$SCode == "ref"),
      "S" = which(batch_metadata$SCode != "ref")
    )
    strategies <- rep(c("Q", "S"), each = length(conditions))
    
    #Perform batch correction ONLY on numeric metabolite columns
    print("Test 1")
    BatchCorrResults <- lapply(seq(along = experiments), function(ii)
      apply(data_BatchCorr[ , metCols], 2, doBC,
            ref.idx   = refSamples[[ strategies[[ii]] ]],
            batch.idx = batch_metadata$Batch,
            minBsamp  = minBatchOccurrence.Line,
            seq.idx   = batch_metadata$Seq,
            method    = methods[ii],
            imputeVal = imputeValues[ii])
    )
    print("Test 2")
    
    #Naming list items with the type of correction
    names(BatchCorrResults) <- experiments
    
    #Anti-log of corrected results
    BatchCorrResultsAntilog <- BatchCorrResults
    BatchCorrResultsAntilog$Q <- 10^(BatchCorrResultsAntilog$Q)
    BatchCorrResultsAntilog$S <- 10^(BatchCorrResultsAntilog$S)
    
    #Write corrected results to CSV files
    write.csv(BatchCorrResults,        file = "BatchCorrResults.csv")
    write.csv(BatchCorrResultsAntilog, file = "BatchCorrResultsAntilog.csv")
    
  })
  
  ####6. Reporting####
  #Create reactive variable for storing the Matrix to display it.
  MetaboloMatrix <- reactiveVal(NULL)
  
  #####6.1 Generating Data Matrix for MetaboAnalyst#####
  observeEvent(input$generate_matrix, {
    
    #Initialize environment
    data <- peak_areas_data()
    if (is.null(data) || nrow(data) == 0) {
      showNotification("No data available to generate matrix", type = "error")
      return(NULL)}
    results_path <- main_folder2()
    
    #Remove unnecessary columns, transpose the data, and acquire sampleIDs 
    data <- data[, !names(data) %in% c("peak.number", "file.name")]
    transposed_data <- t(data)
    transposed_df <- as.data.frame(transposed_data)
    sampleIDs <- transposed_df["PBM Sample ID", ]
    
    #Parse user-specified categories
    user_categories <- strsplit(input$category_input, ",")[[1]]
    user_categories <- sapply(user_categories, function(x) {
      kv <- strsplit(x, "=")[[1]]
      setNames(trimws(gsub("'", "", kv[2])), trimws(kv[1]))
    }, USE.NAMES = FALSE)
    
    #Create CLASS row based on dynamic criteria specified by text input. Allows for dynamic categorizers
    class_row <- sapply(1:ncol(transposed_df), function(i) {
      entry <- sampleIDs[i]
      
      #Check if entry is numeric and if a numeric category is specified
      if ("Numeric" %in% names(user_categories) && suppressWarnings(!is.na(as.numeric(entry)))) {
        return(user_categories[["Numeric"]])
      }
      
      #Check user-specified categories from text input
      for (pattern in names(user_categories)) {
        if (pattern != "Numeric" && grepl(pattern, entry, ignore.case = TRUE)) {
          return(user_categories[[pattern]])
        }
      }
      
      #If a class assignment is missing, assign a blank value to it.
      return("")
    })
    
    #Add class row to transposed data.
    transposed_df <- rbind(class_row, transposed_df)
    rownames(transposed_df)[1] <- "CLASS"
    
    #fix order so CLASS is second row
    transposed_df <- transposed_df[c("PBM Sample ID", "CLASS", rownames(transposed_df)[!rownames(transposed_df) %in% c("PBM Sample ID", "CLASS")]),]
    
    #Save Data
    MetaboloMatrix(transposed_df)
    tryCatch({
      write.table(transposed_df, file.path(results_path, "Metaboanalyst_matrix.csv"), row.names = TRUE, col.names = FALSE, sep = ",", fileEncoding = "UTF-8")
      showNotification("Matrix Made", type = "message")
    }, error = function(e) {
      showNotification("Failed to save matrix: " + e$message, type = "error")
    })
  })
  
  
  #Render the Matrix 
  output$matrix_table <- renderDataTable({
    datatable(MetaboloMatrix(), options = list(
      scrollX = TRUE,
      scrollY = "400px",
      paging = FALSE
    ))
  })
  
  
  
  ####7. Drug Analysis Tab####
  #for now i think I will start from the adherence data structure -> will build in the adherence code later?
  #7.1 Load data and define reactive values
  
  #Reactive values
  adherence_patient_data <- reactiveVal(NULL)
  drugbank_db <- reactiveVal(NULL)
  twosides_db <- reactiveVal(NULL)
  
  
  #Load adherence data
  observeEvent(input$adherence_results,
               {
                req(input$adherence_results)
                 #load data from uploaded file
                 adherence_patient_data(readxl::read_excel(input$adherence_results$datapath))
               })
  

  #Load drug data 
  observeEvent(input$drug_databases, {
                 #Load TWOSIDES database
    twosides_db(parseTWOSIDES("Drug_Data/TWOSIDES.csv.gz"))
    print("Loaded TWOSIDES")
    
    print("Checking for Drug Databank")
    
    #Check for .rds DrugBank file
    #If present, load into reactive value
    if (file.exists("Drug_Data/drugbank_db.rds")){
      drugbank_db(readRDS("Drug_Data/drugbank_db.rds"))
      
      print(".RDS found. Data loaded")
    } else{
      print("Data not found, parse?")
      #Show modal for asking if users want to parse the database
      showModal(modalDialog(
        title = "DrugBank Database Not Found", "Parse now?",
        tags$br(), tags$br(),
        tags$b("Estimated time: 30-40 min."),
        footer = tagList(
          modalButton("No (skip)"),
          actionButton("confirm_parse_drugbank", "Yes, parse now",class = "btn-primary")
        )
      ))
    }
  })
  
  observeEvent(input$confirm_parse_drugbank,{
    removeModal()
    
    #Initialize waiter loading screen
    initial_message <- sample(loading_messages, 1)
    waiter_show(
      html = tagList(
        spin_flower(),
        h4(id = "loading-message", initial_message, style = "color:black;")
      ),
      color = "rgba(210,210,210,0.7)"
    )
    
    #Java script for switching between different loading messages
    js_script <- sprintf("
    setTimeout(function() {
      var messages = %s;
      function updateMessage() {
        var elem = document.getElementById('loading-message');
        if (elem) {
          var newMessage = messages[Math.floor(Math.random() * messages.length)];
          elem.innerText = newMessage;
        }
      }
      setInterval(updateMessage, 4000); // Update every 5 seconds
    }, 2000);
  ", jsonlite::toJSON(loading_messages, auto_unbox = TRUE))
    runjs(js_script)
    
    #Exclude drugabnk catagories we dont need
    exclude <- c("manufacturers", "prices", "atc_codes", "international_brands", 
                 "mixtures", "packagers", "calculated_properties", "pdb_entries","patents","external_links", "atc_codes","sequences","external_identifiers","pdb_entries","ahfs_codes","snp_adverse_reactions", "snp_effects")
    
    #parse drugbank
    drugbank <- parseDrugBank(
      db_path = "Drug_Data/full database.xml",
      drug_options = setdiff(drug_node_options(), exclude),
      parse_salts = FALSE,
      parse_products = FALSE,
      references_options = NULL,
      cett_options = c("targets", "enzymes", "transporters")
    )
    
    #Save file
    saveRDS(drugbank, "Drug_Data/drugbank_db.rds")
    drugbank_db(drugbank)
    
    #hide waiter and confirm
    waiter_hide()
    showNotification("DrugBank loaded and saved successfully!", type = "message")
  })
  
    
}#Closing bracket
  
