# Existing-analysis cleanup: Engine keeps the run picker and uses sidebar navigation;
# Visualization hides its redundant file picker for saved runs. New analyses and
# the retained pseudo-file workflow keep file selection. The folder reminder is removed.
# Clear/reselect fix: run and analysis dropdowns start with an explicit empty
# option so Shiny reports cleared selections and allows choosing the same run again.
# Startup/display fix: Engine uses a normal layout container because it has no
# sub-tabs. Real tab sets elsewhere retain their existing panels and controls.
# Open existing analysis selects a run folder under Results; no raw mz5 selection
# or upload is needed. Visualization selects a saved analysis for that run.

# Analysis modes added: New analysis uploads real mz5 files; Open existing analysis
# uses the filenames recorded inside the selected saved run folder.
# The legacy pseudo-file controls are retained and styled green for later removal.
# Use this ui.R together with server.R from the same analysis-mode branch.

#PeakMeister2.0_UI was written by Erick Helmeczi and Liam Surry 
#Initialized 2024-11-12 @McMaster University
#Department of Chemistry and Chemical Biology 

#### 1. Setup ####
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pacman::p_load("shiny","shinyBS", "tools","waiter","shinyFiles","RColorBrewer","shinycssloaders", "shinydashboard", "DT", "shinyalert","ggpubr", "plotly","hash","pracma", "tidyverse", "stats", "DescTools", "xcms", "rlang","markdown", "openxlsx","readxl","writexl","fontawesome", "MSnbase", "mzR","shinyjs","FactoMineR", "factoextra","jsonlite","remotes","dbparser","R.utils","RSQLite","DBI","shinyWidgets",
               install = TRUE)

#Title of app
title <- "PeakMeister v2.0"

#define UI
ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(
    title = title
  ),
  #### 2. Sidebar Main Tabs ####
  #Set up sidebar and tabs
  dashboardSidebar(
    sidebarMenu(id = "main_navigation",
      #About tab -> for README and instructions
      menuItem("About",
        tabName = "about",
        icon = icon("info-circle")),
      #Paramaters tab -> user supplied mass list and instrument settings
      menuItem("User Supplied Parameters",
        tabName = "userparameters",
        icon = icon("screwdriver-wrench")),
      #Engine tab -> where the analysis actually occurs and where users input their datafiles 
      menuItem("Engine",
               tabName = "engine",
               icon = icon("heartbeat")),
      #Visualization tab for checking electropherograms
      menuItem("Visualization",
        tabName = "processing",
        icon = icon("area-chart")),
      #Manual adjustment tab for processing data
      menuItem("Downstream Processing",
        tabName = "datamanipulation",
        icon = icon("searchengin")),
      #Reporting tab for exporting results/outputting true detections
      menuItem("Reporting",
        tabName = "reporting",
        icon = icon("file-export")),
      #Drug analysis tab for assesing adherence and outputting results of potential drug interactions
      menuItem("Drug Analysis",
               tabName = "drugs",
               icon = icon("pills"))
    )
  ),
  #### 3. About Main Tab Content####
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper {
          padding-bottom: 100px;
        }
        .selectize-dropdown-content {
        max-height: 200px;
        overflow-y: auto;
        }
      .flower-spinner .dots-container .bigger-dot {
        background-color: black;
      }
      .flower-spinner .dots-container .bigger-dot .smaller-dot {
        background-color: black;
      }
      .flower-spinner .dots-container::before,
      .flower-spinner .dots-container::after {
        background-color: black;
      }
      .box.box-primary {
        border-top: 3px solid #1999CC;
        border-radius: 4px;
        margin-bottom: 8px;
      }
      
      .box-header {
        cursor: pointer;
        padding: 10px 15px;
        background-color: #f8f9fa;
        border-bottom: 1px solid #dee2e6;
      }
      
      .box-title {
        font-size: 15px;
        font-weight: normal;
        font-family: inherit;
        color: #FFFFFF;
      }
      
      .box {
        margin-bottom: 6px !important;
      }
      
      .box-tools .btn {
        color: #1999CC;
      }
      "))
    ),
    tags$style(HTML("
      #pseudo_file_workflow,
      #pseudo_file_workflow h4,
      #pseudo_file_workflow label,
      #pseudo_file_workflow .form-control,
      #pseudo_file_workflow .btn,
      #pseudo_file_workflow table th,
      #pseudo_file_workflow table td,
      #pseudo_file_workflow input,
      #pseudo_file_workflow select { color: #188038 !important; }
    ")),
    useShinyjs(),
    use_waiter(),

    tabItems(
      #About tab content
      tabItem(
        tabName = "about",
        fluidPage(
          column(width = 12, 
                 h3("Documentation"),
                 box(title = "Disclaimer", status = "primary", solidHeader = TRUE,
                     width = NULL, collapsible = TRUE, collapsed = TRUE,
                     uiOutput("disclaimerContent")),
                 box(title = "README", status = "primary", solidHeader = TRUE,
                     width = NULL, collapsible = TRUE, collapsed = FALSE,
                     uiOutput("readmeContent")),
                 box(title = "Updates", status = "primary", solidHeader = TRUE,
                     width = NULL, collapsible = TRUE, collapsed = TRUE,
                     uiOutput("updatesContent")),
                 box(title = "License", status = "primary", solidHeader = TRUE,
                     width = NULL, collapsible = TRUE, collapsed = TRUE,
                     uiOutput("licenseContent"))
              )
            )
        ),
      
      #### 4. User Supplied Inputs Main Tab Content####
      tabItem(
        tabName = "userparameters",
        tabsetPanel(
          id = "userParameters",
          #make subtabs for Mass List, Refererence List, & Paramaters. 
          #Option to either upload and excel file or for you to manually create your own through the app 
          
          #####4.1. Mass List Subtab##### 
          tabPanel("Mass List",
                   fileInput("massList", "Upload Mass List and Parameters (Excel)", accept = c(".xlsx")),
                   div(style = "overflow-x: auto;",
                   DT::dataTableOutput("massListData")),
                   fluidRow(
                     #Create input columns for required inputs from excel file
                     column(3, textInput("name", "Name")),
                     column(3, numericInput("mz", "mz", value = 0)),
                     column(3, numericInput("extraction.window.ppm", "Extraction Window (ppm)", value = 0)),
                     column(3, textInput("interference", "Interference")),
                     column(3, numericInput("interference.comigration.threshold.seconds", "Interference Comigration Threshold (seconds)", value = 0)),
                     column(3, numericInput("minimim.peak.width.seconds", "Minimum Peak Width (seconds)", value = 0)),
                     column(3, numericInput("migration.window.seconds", "Migration Window (seconds)", value = 0)),
                     column(3, numericInput("peak.space.tolerance.percent", "Peak Space Tolerance (%)", value = 0)),
                     column(3, numericInput("snr.threshold", "SNR Threshold", value = 0)),
                     column(3, textInput("smoothing.kernel", "Smoothing Kernel")),
                     column(3, numericInput("smoothing.strength", "Smoothing Strength", value = 0))
                   ),
                   #For organizational purposes, separate peak migration time inputs from other paramaters
                   tags$strong(h4("Peak Migration Time Inputs")),
                   #Function for creating a column for each peak. 
                   #May include option to change the number of peaks in the future. For now, its hardcoded at 13
                   fluidRow(
                     lapply(1:13, function(i) {
                       column(3, numericInput(paste0("peak.", i, ".mt"), paste("Peak", i, "MT"), value = 0))
                     })
                   ),
                   #Button for adding a row
                   actionButton("addMassRow", "Add Mass"),
                   #Button for Deleting a selected row
                   actionButton("deletemassrow", "Delete Selected Row"), 
                   #Button for clearing entire table 
                   actionButton("clearMassList", "Clear All")),
          
          #####4.2 Reference Mass List Subtab#####
          tabPanel("Reference Mass List",
                   fileInput("massList", "Upload Mass List and Parameters (Excel)", accept = c(".xlsx")),
                   div(style = "overflow-x: auto;",
                   DT::dataTableOutput("refMassListData")),
                   fluidRow(
                     #generate input columns for reference mass list
                     column(3, textInput("name", "Name of Reference Mass")),
                     column(3, numericInput("mz", "mz (reference mass)", value = 0)),
                     column(3, numericInput("extraction.window.ppm", "Extraction Window (ppm)", value = 0)),
                     column(3, textInput("class", "Class")),
                     column(3, numericInput("min.mt.min", "Min MT Min", value = 0)),
                     column(3, numericInput("max.mt.min", "Max MT Min", value = 0)),
                     column(3, numericInput("peak.fwhm.tolerance.multiplier", "Peak FWHM Tolerance Multiplier", value = 0)),
                     column(3, numericInput("peak.space.tolerance.percent", "Peak Space Tolerance (%)", value = 0)),
                     column(3, textInput("smoothing.kernel", "Smoothing Kernel")),
                     column(3, numericInput("smoothing.strength", "Smoothing Strength", value = 0))
                   ),
                   #For organizational purposes, separate peak migration time inputs from other paramaters
                   tags$strong(h4("Peak Migration Time Inputs")),
                   #Function for creating a column for each peak. 
                   #May include option to change the number of peaks in the future. For now, its hardcoded at 13
                   fluidRow(
                     lapply(1:13, function(i) {
                       column(3, numericInput(paste0("ref.peak.", i, ".mt"), paste("Peak", i, "MT"), value = 0))
                     })
                   ),
                   #Button for adding a row 
                   actionButton("addRefMassRow", "Add Reference Mass"),
                   #Button for deleting reference mass row 
                   actionButton("deleterefmassrow", "Delete Selected Row"),
                   #Button for clearing whole table 
                   actionButton("clearRefMassList", "Clear All")),
                   
          #####4.3 Parameters list subtab#####
          tabPanel("Parameters",
                   fileInput("massList", "Upload Mass List and Parameters (Excel)", accept = c(".xlsx")),
                   div(style = "overflow-x: auto;",
                   DT::dataTableOutput("parametersData")),
                   fluidRow(
                     #Generate input columns for parameters datatable
                     column(3, numericInput("number.of.injections", "number.of.injections", value = 1)),
                     column(3, numericInput("ref.mass.one", "ref.mass.one", value = 0)),
                     column(3, numericInput("ref.mass.two", "ref.mass.two", value = 0)),
                     column(3, numericInput("ref.mass.window.ppm", "ref.mass.window.ppm", value = 0)),
                     column(3, numericInput("ref.mass.minimum.counts", "ref.mass.minimum.counts", value = 0)),
                     column(3, selectInput("apply.mass.calibration", "apply.mass.calibration", choices = c("Yes", "No"))),
                     column(3, selectInput("apply.smoothing", "apply.smoothing", choices = c("Yes", "No"))),
                     column(3, selectInput("plot.format", "plot.format", choices = c("Sample", "Metabolite"))),
                     column(3,selectInput("Manual.Indexes","Manual.Indexes", choices = c("Yes", "No")))
                   ),
                   #Action button for adding rows
                   actionButton("addParamRow", "Add Parameter"),
                   #Button for deleting row
                   actionButton("deleteparamrow", "Delete Selected Row"),
                   #Button for clearing table 
                   actionButton("clearparam", "Clear All")),
          
          #####4.4 Project Information Subt ab#####
          tabPanel("Project Information",
                   fluidPage(
                     #Set HTML formating fonts and line spacing and allow for text wrapping
                     tags$style(HTML("pre {
               font-family: Arial, sans-serif;
               line-height: 2.0;
               white-space: pre-wrap; /* Allows wrapping */}")),
                     #Generate input columns 
                     column(width = 12,
                            h3("Project Information"),
                            fluidRow(
                              column(3, textInput("projectName", "Project Name (no spaces -> REQUIRED)")),
                              column(3, textInput("projectDescription", "Project Description")),
                              column(3, textInput("projectSupervisor", "Project supervisor (no spaces)")),
                              column(3, textInput("projectContact", "Project Contact (email)"))
                                    ),
                            #Button for saving project info to reactive value -> will be used to generate a metadatafile later on
                            actionButton("saveProjectInfo", "Save Project Info"),
                            #Button for clearing project information
                            actionButton("clearproject", "Clear All"),
                            #text output for project info 
                            tags$hr(),
                            h4("Current Project Info"),
                            verbatimTextOutput("currentProjectInfo")
                          )
                        )
                      ),
          
          #####4.5: Manual Migration Index Input#####
          tabPanel("User Supplied Migration Indexes",
                   fileInput("userMTI", "User Supplied Migration Indexes (Excel)", accept = c(".xlsx")),
                   div(style = "overflow-x: auto;",
                       DT::dataTableOutput("userMTIdata")),
                   fluidRow(
                     #generate input columns for reference mass list
                     column(3, textInput("name", "Name of Metabolite")),
                     column(3, textInput("left_is", "Left internal standard")),
                     column(3, textInput("right_is", "Right internal sdtandard")),
                     column(3, textInput("description", "Description")),
                   ),
                   #For organizational purposes, separate peak migration time inputs from other paramaters
                   tags$strong(h4("Peak Migration Time Inputs")),
                   #Function for creating a column for each peak. 
                   #May include option to change the number of peaks in the future. For now, its hardcoded at 13
                   fluidRow(
                     lapply(1:13, function(i) {
                       column(3, numericInput(paste0("ref.peak.", i, ".mt"), paste("Peak", i, "MT"), value = 0))
                     })
                   ),
                   #Button for adding a row 
                   actionButton("addMTIrow", "Add manual MTI row"),
                   #Button for deleting reference mass row 
                   actionButton("deleteMTIrow", "Delete Selected MTI row"),
                   #Button for clearing whole table 
                   actionButton("clearMTIlist", "Clear All")
                   )
          
                  ),
        #Action button for saving manually input data
        downloadButton("downloadData", "Save Table (.xlsx)")
        ),
      
      ####5. Engine####
      tabItem(
        tabName = "engine",
        div(
          id = "Setup",
          fluidRow(
            column(6,
                   radioButtons("analysis_mode", "Choose analysis mode",
                                choices = c("New analysis" = "new", "Open existing analysis" = "existing"),
                                selected = "new"),
                   helpText("Switching modes clears the selected file list. Saved results are kept."),
                   conditionalPanel("input.analysis_mode === 'new'",
                     tags$h4(style = "text-decoration: underline;", "File Uploads for Initializing Data Preprocessing"),
                     fileInput("mz5Files", "Upload .mz5 Files", accept = c(".mz5"), multiple = TRUE),
                     checkboxInput("autorun", "Start run automatically after upload", value = FALSE)
                   ),
                   conditionalPanel("input.analysis_mode === 'existing'",
                     tags$h4("Open existing analysis"),
                     selectInput("existing_run", "Select run folder in Results:",
                                 choices = c("Select a run" = "")),
                     helpText("Choose a saved run. The original .mz5 file is not needed.")
                   ),
                   actionButton("clearFiles", "Clear Files"),
            ),
            column(1,),
            column(5,
                   div(id = "pseudo_file_workflow",
                   tags$h4(style = "text-decoration: underline;", "Generate Pseudo .mz5 Files to Access Stored Datafiles"),
                   fileInput("excelFile", "Upload Excel File of File Names (.xlsx)", accept = c(".xlsx")),
                   DTOutput("mz5fileTable"),
                   actionButton("addRow", "Add Empty Row"),
                   actionButton("generateFiles", "Generate Pseudo .mz5 Files")
                   )
            )
          )
        ),
        conditionalPanel("input.analysis_mode === 'new'",
        fluidRow(
          column(6,
                 #Big red button for initializing the run.This is where most of the calculations are done
                 actionButton("initialize", "Initialize Run", 
                              style = "background-color:#1999CC;
                      color:#000000;
                      border-color:#000000;
                      border-style:double;
                      border-width:10px;
                      border-radius:0%;
                      font-size:40px; width: 100%;")),
              )
        ),
       DT::dataTableOutput("fileTable"),
      #Create progress bar for monitoring computation 
       uiOutput("progressBar")
      ),
      
      
      ####6. Visualization tab####
      tabItem(tabName = "processing",
              fluidRow(
                column(4,
                       #Option for uploading a saved ,RData file to allow users to edit their data at another time generated by the app -> requires users to also upload their files so file_list is populated 
                       fileInput("RDatainput", "Uploaded .RData file", accept = c(".RData")),
                       #Button for reseting uploaded plot data
                       actionButton("resetplotlydata","Clear uploaded data"),
                       #Button for selecting files -> will subset/filter data to choose plots
                       conditionalPanel("output.show_visualization_file_selector === 'true'",
                         selectInput("file_selector", "Select File:", choices = NULL)),
                       
                       #Select results folder to work out of 
                       selectInput("results_folder", "Select saved analysis for this file:",
                                   choices = c("Select a saved analysis" = "")),
                       textOutput("results_folder_status"),
                       #Select plots, display plots
                       tags$h4(style = "text-decoration: underline;", "Loaded Plots"),
                       #Datatable for displaying modified plots
                       DTOutput("plot_table"), 
                       tags$h4(style = "text-decoration: underline;", "Edited plots to be regenerated"),DTOutput("modified_peak_plots_table")),
                
                column(8,
                       tabsetPanel(
                         tabPanel("Plotting",
                       #Display selected plot 
                       plotlyOutput("selected_plot", height = "600px"),
                       #toggle integration mode to allow for zooming and drawing on the graph
                       shinyWidgets::materialSwitch("draw_mode", "Draw Mode", status = "primary"),
                       #Line positions
                       textOutput("red_line_position"),
                       textOutput("blue_line_position"),
                       #HTML tags for setting button colour for aesthetics 
                       tags$style(HTML("
                        #adjust_all_baselines {background-color: #71A9F7; color: black;}
                        #delete_peak {background-color: #2E86AB; color: black;}
                        #undo {background-color: #243B4A; color: white;}
                        #other_options {background-color: #0E3B43; color: white;}
                        #manual_integrate {background-color: #00A8E8; color: black;}
                        #adjust_indiv_baseline {background-color: #28536B; color: white;}
                        ")),
                       #Action button for deleting peaks
                       actionButton("delete_peak", "Delete Selected Peaks"),
                       #Action button for undoing peak deletion
                       actionButton("undo", "Undo Peak Deletion"),
                       #Action button for manual integration 
                       actionButton("manual_integrate", "Manually Adjust Integration"),
                       #Action button for adjusting baseline (across individual peaks)
                       actionButton("adjust_indiv_baseline", "Adjust Individual Baseline"),
                       #Action button for adjusting baseline (across all peaks)
                       actionButton("adjust_all_baselines", "Adjust All Baselines"),
                       #Action button for regenerating plots
                       # actionButton("regenerateplots", "Regenerate Changed Plots"),
                       actionButton("other_options", "Other"),
                       #action button for interferences
                       actionButton("interferences","Mark Peak Interference"),
                       #Action button for deleting all but the first and last peaks.
                       actionButton("deletemiddlepeaks","Delete Middle Peaks (Keeps first and last)"),
                       #Action button for deleting all peaks
                       actionButton("deleteallpeaks", "Delete all peaks"),
                       #Data table for peak positions
                       tags$h4(style = "text-decoration: underline;", "Peak Position Information"),
                       DTOutput("peak_info_table"),
                         ),
                       
                       tabPanel("Dual EIE View",
                                tags$h4(style = "text-decoration: underline;", "Combined Electropherograms"),
                                #Output for combined electropherograms
                                plotlyOutput("combined_plot", height = "600px"),
                                tags$h4(style = "text-decoration: underline;", "Select Plots to Compare"),
                                #Datatables for selecting plots for the combvined electropherograms
                                column(4,
                                DTOutput("plottable1")),
                                column(4,
                                DTOutput("plottable2")),
                                column(4,
                                DTOutput("plottable3")),
                    )
                  )
                )
              )
            ),
      
      #### 7. Downstream processing tab####
      tabItem(
        tabName = "datamanipulation",
          fluidRow(
            column(3,
                   tags$h4(style = "text-decoration: underline;", "Set working environment and initialize data"),
                   #Select results folder to work out of 
                   selectInput("results_folder2", "Select Results Folder:",
                               choices = c("Select a saved analysis" = "")),
                   #Button for loading migration data and peak area data to reactive value
                   actionButton("load_migration_area", "Load Migration and Peak Area Data"),
                   #m/z vs MT controls
                   tags$h4(style = "text-decoration: underline;", "m/z vs MT controls"),
                   #Button for selecting files -> will subset/filter data to choose plots
                   selectInput("file_selector2", "Select File:", choices = NULL),
                   #Button for selecting peak for m/z versus migration time plot
                   uiOutput("peak_number_selector"),
                   #Manipulate Data
                   tags$h4(style = "text-decoration: underline;", "Editing your data matrix"),
                   #Option to select a missing data correction method
                   selectInput("missing_data_method", "Select Method for Replacing Missing Data", choices = c("Minimum Values/5", "Missing values = 0")),
                   #Action button for replacing missing data with selected button
                   actionButton("missing_data", "Replace Missing Data"),
                   #button for normalizing data to a selected metabolite
                   actionButton("normalize", "Normalize Data"),
                   #Control charts
                   tags$h4(style = "text-decoration: underline;", "Generate Control Charts"),
                   selectInput("selected_metabolite", "Select Metabolite", choices = NULL),
                   actionButton("create_control_chart", "Create Control Chart"),
                   #Metadata
                   tags$h4(style = "text-decoration: underline;", "Connecting metadata"),
                   #Button for connecting metadata to results
                   actionButton("connect_metadata","Connect Metadata to Results")
                   
            ),
            column(9,
                   tabsetPanel(
                     #Plotting m/z versus migration time
                      tabPanel("m/z vs Migration Time",
                            plotlyOutput("mz_vs_migration_plot", height = "600px")),
                      #Generating and plotting control charts
                      tabPanel("Control Charts",
                               
                              tags$h4(style = "text-decoration: underline;", "Control chart for non-normalized data"),
                              plotlyOutput("control_chart_raw"),
                              tags$h4(style = "text-decoration: underline;", "Control chart for normalized data"),
                              plotlyOutput("control_chart_processed")),
                      #Analyzing variance (CVs and PCAs) 
                      tabPanel("Data Quality & Variability",
                               
                              fluidRow(
                                column(3,
                                       
                                      textInput("category_input", "Specify QC identifier (e.g., QC=NIST)", value = "QC="),
                                      HTML("<p><strong>Note:</strong> Identifier must match label in metadata.<p>"),
                                      textInput("exclude_ids", "Exclude PBM Sample IDs (e.g., SP1,SP2,Blank)", ""),
                                      HTML("<p><strong>Note:</strong> This allows you to exclude samples from your CV and PCA assessment. Exclusions are case sensitive and comma-separated.<p>"),
                                      actionButton("compute_cv", "Calculate CVs"),
                            
                                      tags$h4(style = "text-decoration: underline;", "Average CVs (%)"),
                                      tableOutput("average_cv_table")),
                                
                                column(9,
                                      plotOutput("pca_plot")
                                      )
                                    ),
                              fluidRow(
                                column(12,
                                      tags$h4(style = "text-decoration: underline;", "CV based on m/z of metabolite"),
                                      plotlyOutput("cv_plot")
                        )
                      )
                    ),
                      tabPanel("Metadata",
                               tags$h4(style = "text-decoration: underline;", "Peak Area Datatable"),
                             DTOutput("peak_area_table")
                    ),
                    
                    tabPanel("Batch Correction",
                           #initialize environment
                           actionButton("load_batchcorr_packages", "Initialize Packages"),
                           fileInput("batch_corr_order", "Upload Batch Correction Order", accept = c(".csv")),
                           actionButton("batch_corr", "Perform Batch Correction")
                    )
                  )
                )
              )
            ),
      
      ####8. Reporting Tab####
      tabItem(
        tabName = "reporting",
        fluidRow(
          column(3,
                 tags$h4(style = "text-decoration: underline;", "Generate Matrix for MetaboAnalyst"),
                 actionButton("generate_matrix","Generate Data Matrix for MetaboAnalyst"),
                 textAreaInput("category_input", "Specify classes for MetaboAnalyst matrix", rows = 5),
                 HTML("<p><strong>Note:</strong> Class assignments for your matrix can be both specific and/or generic. Specific class assignments include assignments such as 'SP1 = Standard1', 'SP2 = Standard2', and 'NIST = SRM30061'. Generic assignments would include assignments such as assigning all numeric values to Positive_Samples such as 'Numeric = Positive_Samples'.<br><br>An example of a set of class assignments would be: 'SP1=Standard1,SP2 =Standard2, NIST = SRM361, Numeric=Pos_Sample' <br><br>Unassigned sample classes are left blank. Class names are case sensitive, and multiple class specifications are separated by a comma. Class names must be one word (no spaces)</p>")
          ),
          column(9,
                 tags$h4(style = "text-decoration: underline;", "MetaboAnalyst Matrix"),
                 dataTableOutput("matrix_table")
                )
              )
            ),
      
      ####9. Drug Analysis Tab####
      tabItem(
        tabName = "drugs",
        fluidRow(
          column(3,
                 #load adhernece results if not already present
                 fileInput("adherence_results", "Upload Adherence Results (Excel)", accept = c(".xlsx")),
                 #load drugbank data and drug interaction data
                 HTML("<p><strong>Note:</strong> Users must acquire drug databases independently.<p>"),
                 actionButton("drug_databases", "Load Drug Databases")
                 )
        )
      )
      
          ),
    #Display logo at bottom right of page 
    tags$div(
      style = "position: absolute; bottom: 10px; right: 10px;",
      tags$img(src = "PeakMeisterLogo.png", style = "width: auto; height: auto;")
    )
  )
)#Closing Bracket for UI dashboiard
