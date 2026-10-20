---
title: "PeakMeister2.0"
author: "Erick Helmeczi & Liam Surry"
bibliography: bibliography.bib
output: html_document
---


## Software for Pre-Processing MSI-CE-MS Datasets

### **Description**
***

Multisegment injection-capillary electrophoresis-mass spectrometry (MSI-CE-MS) is a robust analytical technique capable of rapidly acquiring metabolomic data (<4min/sample) via serial injections of samples and quality controls within every analytical run. The high-throughput nature of this platform makes it an attractive technique for large scale metabolomic experiments where keeping costs down and meeting project deadlines are of utmost importance. However, until now, data collected by MSI-CE-MS had to be manually pre-processed by an experienced analyst due to the high migration time variability of CE-MS and the complexity of multiplexed data sets, which reduce the compatibility of this technique with other pre-processing tools. Unfortunately, manually pre-processing data is slow, expensive, and tedious, ultimately decreasing the merits of this technique. Thus, we are introducing PeakMeister, an open-source software written in the R statistical environment for the automated pre-processing of targeted full-scan MSI-CE-MS data.

### **Current Features**
***

The key differentiating feature of PeakMeister from other currently available software tools for pre-processing metabolomic datasets is its use of migration indexes to predict the elution time of analytes in MSI-CE-MS datasets. Thus, to achieve correspondence between analytical runs, PeakMeister does not perform migration time alignments or converting the time dimension to electrophoretic mobility. Instead, PeakMeister computes migration indexes for each analyte and uses the migration time of internal standards, or any reliable signals with sufficient signal-to-noise ratios, to compute the migration times of analytes which can then be used for peak annotation and integration. Additionally, as PeakMeister was designed for multiplexed datasets, it also uses the spaces between analytical peaks to confirm and adjust peak annotation, as the gaps between peaks are typically consistent in MSI-CE-MS experiments. Results produced by PeakMeister are saved and include:

1. A table containing the migration indexes or relative migration times used to annotating analytes
2. A copy of the parameters used to process the data
3. Plots of each extracted ion electropherogram which can be used to check for proper peak annotation and integration
4. A table containing the migration times of peaks or expected peak positions
5. A table containing the peak areas
6. The raw metadata used by the program which can be used further in data analysis steps

Migration time indexes are determined with the following equations


$$
MTI = \frac{{MT_{Metabolite} - MT_{IS_1}}}{{MT_{IS_2} - MT_{IS_1}}}
$$ 


Where:

1. $MTI$ = The Migration time index
2. $MT_{Metabolite}$ = Metabolite migration time 
3. $MT_{IS_1}$ = Internal Standard 1 migration time
4. $MT_{IS_1}$ = Internal Standard 2 migration time

MTI's are then used to select peaks based on their predicted migration times (MT's) within a defined time window using the following equation 

$$
MT_{Metabolite} = MTI \left( (MT_{IS_2}) - (MT_{IS_1}) \right) + (MT_{IS_1})
$$

Where 

1. $MT_{Metabolite}$ is the predicted migration time of the metabolite

This app would not be possible without @helmeczi_software_2024, @deutsch_file_2012, and @wilhelm_mz5_2011 
### **How to use**
***
PeakMeister2.0-UI includes the use of a graphical user interface (GUI) run by shiny intended to be used with zero knowledge of code. This requires the user to have [R](https://www.r-project.org/) and [R studio](https://posit.co/download/rstudio-desktop/) installed. PeakMeister2.0 UI requires [R version 4.2](https://mirror.csclub.uwaterloo.ca/CRAN/) or [R version 4.3](https://mirror.csclub.uwaterloo.ca/CRAN/). Due to outdated packages, PeakMeister2.0_UI will not work on R version 4.4.

#### **Setting up**

1.	Create a new project in R studio. This can be done by clicking the dropdown menu in the top right corner of R studio or clicking File -> New project

2.	In “New Project,” select version control and click “Git – Clone a project from a Git repository”

3.	Select where on your computer you want this project to be stored and paste the following URL: https://github.com/liamthepineapple/PeakMeister2.0_UI.git

4.	This app is designed in a way where you never need to go into the complex code. Open the file titled app.R. This is how you launch the app.

5.	Select all on this page and hit “control-enter” to run all selected lines of code. Alternatively, you can click the “Run app” button in the top right corner of your script. This will initialize the app and install the required packages.

6.	Occasionally, the packages “xcms” and “ggpubr” will not install correctly. The user may have to install these manually. To install these packages: 

#### xcms
```r 
BiocManager::install(“xcms”)
library(xcms)
```
#### ggpubr
```r
install.packages(“ggpubr”)
library(ggpubr)
```

Note: The user may require [Rtools](https://cran.r-project.org/bin/windows/Rtools/) to be installed which must be the same as the version of R installed.

PeakMeister has only two requirements to get you up and running and all required template files and folders are included in the latest releases:
  * Convert all your data files to open-source mz5 files using [ProteoWizard](https://proteowizard.sourceforge.io/)
  * Provide a targeted mass list and the corresponding parameters using the provided "Mass List and Parameters.xlsx" template

Use the project file to open the R script titled "app.R" and execute the script to begin pre-processing your data

### **Detailed Usage**
***

Although PeakMeister is a R-based software tool, users will typically require little to no knowledge of coding to pre-process their data as an excel sheet containing all pre-processing parameters is provided. In this detailed usage overview, the purpose of each parameters and how they can be manipulated to accurately annotate MSI-CE-MS data will be explained.

#### **Mass List and Paramaters File**

##### Sheet 1: *Mass List*

1. name - Provide a name for each analyte in your study. All names must be unique, so if you have multiple unknown compounds, use names such as Unknown-1, Unknown-2, Unknown-3, etc.
2. mz - Provide the mass-to-charge to be extracted for each analyte in your study.
3. extraction.window.ppm - This is the extraction window used to extract each mz value provided with units of ppm. We found that a minimum mass window of ~30-35 ppm is required depending on the mass of the analyte, however this will likely be dependent on the mass spectrometer used duing data acquisition.
4. interference - Here you can designate analyte interferences. Write the name of the analyte or internal standard that is a significant interference and PeakMeister will annotate any instances of overlap.
5. interference.comigration.threshold - This parameter only needs to be set when an interference is also provided. This is the window used to determine if two peaks are overlapping enough to be considered an interference. The widths of your peaks will determine how high this value needs to be set.
6. minimum.peak.width.seconds - This parameter is used to filter out noise during peak detection. Only peaks with a width equal to or greater than this value will be considered during annotation.
7. migration.window.seconds - After migration time prediction, PeakMeister will look for analyte peaks at the expected migration time +/- the time set for this parameter.
8. peak.space.tolerance.percent - Peaks are expected to be spaced approximaely equally apart. This parameter designates how differently they can deviate from the median spacing before being subject to reanalysis by PeakMeister. 
9. snr.threshold - After identifiying which signals correspond to the analyte peaks, PeakMeister will compute their S/N to determine if they should be recorded as "<LOD".
10. smoothing.kernal - Smoothing is performed using the "ksmooth" function from the [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/ksmooth.html) package. This parameter sets the kernal to be used.
11. smoothing.strength - Smoothing is performed using the "ksmooth" function from the [stats](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/ksmooth.html) package. This parameter sets the bandwidth to be used.
12. The last columns of of this sheet are where peak migration times for each analyte are to be set. All of these values need to be taken from the same run, including those that will be set for the internal standards in the next sheet.

##### Sheet 2: *Reference Mass List*

For this sheet, I will focus on the sections not discussed above. Features used in this section can be used as landmark peaks that can help predict the migration times of analytes in the previous sheet. These landmark peaks use a separte algorithm for annotation which requires their peak areas to be reliable and always the most intense within the EIE. Thus, internal standards are typically used here, however analytes that are very concentrated in your sample matrix may also be used. 

1. class - Use "Reference" for any compound you want to be used to calculate migration indexes for analytes on the previous sheet. Typically this should be every compound listed here, however the option is availble to not use them. Simply type "Not Reference" instead if this is the case.
2. min.mt.min - Minimum migration time cutoff. Peaks that elute before this threshold will not be considered as possible reference standard peaks.
3. max.mt.min - Maximum migration time cutoff. Peaks that elute after this threshold will not be considered as possible reference standard peaks.

##### Sheet 3: *Parameters*

1. number.of.injections - Number of injections used during data acquisition.
2. ref.mass.one - Lower lock-mass used for accurate mass correction.
3. ref.mass.two - Upper lock-mass used for accurate mass correction.
4. ref.mass.window.ppm - mass window to search for lock-mass peak.
5. ref.mass.counts - Minimum peak hight requirement of the lock-mass for correction to be applied.
6. apply.mass.correction - Should the lock mass correction be applied. Correction will only be applied if "Yes" is set for this parameter.
7. apply.smoothing - Should smoothing be applied to the EIEs. Smoothing will only be applied if "Yes" is set for this parameter.

#### **App Tabs**
***

Upon Launching the app, you're taken to a homepage with 6 tabs displayed on the left. The tabs are as following:

##### 1. **About**
***

This tabs contains information about the app including the detailed README file for how the app works. Tabs are dropdown menus that display information once clicked. The following tabs are present

1. Disclaimer -> provides a disclaimer that this app is still under development
2. README -> README on how the app works
3. Updates -> location for any recent updates
4. License -> Includes MIT license

##### 2. **User Supplied Parameters**
***

User supplied parameters is the location where an individual will specify the paramaters used in the PeakMeister computations. Users can upload excel files with the required information (supplied in the Sample Files folder) which will automatically populate the data table with the contained information. If a table is not available, users can input their data using the app interface in addition to saving this data to a prelabeled parameters datafile. Uploaded data can also be edited direcly in the table by double clicking column. The tabs are as following

1. Mass List 

Required for computation to run. Will auto-populate when user uploads Mass List and Parameters Excel File. Contains information on the metabolites to be analyzed. 

2. Reference Mass List

Required for computation to run. Will auto-populate when user uploads Mass List and Parameters Excel File. Contains information on refernece internal standrads/masses to use for migration time index calcuations. 

3. Parameters

Required for computation to run. Will auto-populate when user uploads Mass List and Parameters Excel File. Contains information about specific parameters for running the calculation including reference masses, and number of injections. Manual.Indexes will determine if you need to use the User Supplied Migration indexes tab or not

4. Project Information (not currently functional)

Tab for adding basic project metadata. This includes project name, a description of it, a supervisor, and a contact for who is running/analyzing the samples. This will then be incldued into a file thats included with the output results folder

5. User Supplied Migration Indexes

Only required for when users select "Yes" in the Manual.Indexes option in the Parameters tab. Will auto-populate when users upload excel file titled User Supplied Migration Indexes. Provides option for users to supply their own migration time indexes. 

##### 3. **Engine**
***

This tab is where the major computation is performed. The metaphorical engine of our metaphorical car. This requires .mz5 files to run.Users will click on "Browse" to select a .mz5 file to upload. This file can be stored anywhere on your computer and the app will upload it. When the file is uploaded, shiny stores it until the app is closed. Uploaded files will be displayed in a rendered data table so the user can check the files that they have selected. Once files have been uploaded, users can click the big button labelled "Initialize Run" which will start the button. This button will not do anything until the uploaded .mz5 files are uploaded and additionally requires the Mass List and Parameters to be uploaded. An option exists for a user to check a checkbox that will automatically start the run after the very last data file has been uploaded. This allows the user to upload a large amount of files and walk away from their computer knowing the app will automatically start once the last file is uploaded. The computation will then be performed on your samples. The run will give the following files/folder:

1. Mass List and Parameters (.xlsx)

A copy of the Mass List and Parameters used for the computation

2. Migration Index Summary (.xlsx)

A table containing the migration indexes or relative migration times used to annotating analytes

3. Plots (folder)

Plots of each extracted ion electropherogram which can be used to check for proper peak annotation and integration

4. Metabolite Migration Times (.xlsx) 

A table containing the migration times of peaks or expected peak positions

5. Metabolite Peak Areas (.xlsx)

A table containing the peak areas

6. plotly_objects (.RData)

A Rdata file containing all generated interactive plotly plots. ALlows for visualization after the app has been closed

7. Data (folder)

Stores the raw metadata for each datafile. This includes all information required for generating the plots and calculating peak areas.

##### *Opening an existing analysis*

In **Engine**, choose **Open existing analysis**, then select a run folder under `Results`. The original `.mz5` file is not needed. For example, a run created from `Run01.mz5` is stored as `Results/Run01/Results YYYY-MM-DD/`.

Selecting a run automatically selects and loads its most recent saved analysis. The date in the folder name determines the newest analysis; reruns on the same date are ordered numerically (`Results YYYY-MM-DD 10` is newer than `Results YYYY-MM-DD 2`). The **Select saved analysis for this file** dropdown in **Visualization** allows switching to an older version. A manually selected older version stays selected while the run remains open. If a saved analysis has no SQLite database, the existing database migration prompt is shown.

**Clear Files** clears the active selection without deleting the original files or saved results. The same run can be selected again afterward. No empty `.mz5` files are generated or required.

##### 4. **Visualization**
***

Use **Visualization** in the left sidebar to view the analysis selected in Engine. For an existing analysis, its filename is selected internally and the **Select File** control is hidden. For new analyses with multiple uploaded files, **Select File** remains available to switch between files. Select a plot from **Loaded Plots** to display its electropherogram and peak information.

The saved-analysis dropdown is populated and selected automatically when opening an existing run. It remains available for switching between dated analysis versions. Stored data and plots are loaded from the selected analysis folder.

There are two subtabs to display the plots in which are titled "Plotting" or "Dual EIE View." The default is the "Plotting" tab.

##### *Plotting subtab*

This is the default tab for displaying the plots. The user selects a plot from the "Loaded Plots" table which then causes the selected plot to render. Additionally, the data associated with all the peaks will display underneath the plot with information on the peak number, the apex of the peak in seconds, and the maximum intensity. A number of buttons are available for editing the plots which are described below:

1. *Delete Selected Peaks*

This button is for deleting peaks such as any instances where PeakMeister has misannoatated a peak. Users can select indivdual peaks or multiple peaks at once by selecting the corresponding row in the displayed table of the plot data. Peaks are then deleted by clicking "Delete Selected Peaks" after row selection. This will change the metadata associated with the plot. Deleted peaks are replaced with the label "NPD" which corresponds to "No peak detected." This does **not** directly update the stored peak area information, but it stores the metadata information with the deleted peaks for use in regenerating the peak area excel file at a later step. Upon deleting a peak, the peak name will be marked as "edited" and appear in the table of edited plots to be regenerated. 

2. *Undo Peak Deletion*

This button is for undoing peak deleting **PRIOR** to regenerating the plots. This is used in case the user accidentally deletes the wrong. It will return the metadata to whatever it was prior to the last peak deletion event. Note: this can only undo the previous peak deletion and this will reset if the user switches files.  

3. *Manually Adjust Integration* 

This button is for adjusting the integration of peaks. Occasionally, PeakMeister will either integrate the wrong peak or perform a poor integration of the peak area. The user can specify the left and right boundaries to integrate between by using the red and blue lines that are rendered with each plot. Each of these lines can be dragged and indicates a value on the x axis to integrate between.The red line indicates the left boundary and the blue line indicates the right boundary. These line positions will be saved as the new start and end time points for the peak being integrated. The app then looks into the electropherogram data and integrates the data found between these two plots. The peak apex can be redefined in this step as the program simply looks for the highest point in between these two boundaries and chooses that as the peak apex. The peak being integrated can either be selected from the Peak Position Information table displayed below the plot or the user can enter the peak number into a popup that appears if the user has not selected anything in the table. Only one peak integration can be adjusted at a time, however, the peak integration can be done as many times as a user wants for a plot before clicking "regenerate changed plots." Additionally, users could adjust integration on as many peaks as they wanted for as many plots as they wanted as long as the file isn't changed. Once "Manually Adjust Integration" is clicked, the program will update the associated metadata and saved Peak Area excel file with the new data.     

4. *Adjust Individual Baseline*

This button is for adjusting baselines on a individual peak-by-peak basis. Once clicked, baseline adjustment mode is activated which means that the plot can now be clicked which will record a "x" and "y" position. The "click" will define the y value of the baseline to be used for the peak. The same peak selection logic used for manually adjusting the integration is applied here as well. The peak where you are adjusting the baseline for can either be indicated by selecting it in the datatable below the plot or a popup will appear where the user can manually enter the peak number being adjusted. This function will then reintegrate the peak with the new baseline and update all associated metadata and peak area information. 

5. *Adjust All Baselines*

This button is similar in function to the "adjust individual baseline" button, only this will adjust the baseline on **ALL** peaks. The same logic from the individual baseline adjustment function is applied for choosing the new baseline. This will also reintegrate the associated metadata and update the medata date and peak area excel file.  

6. *Regenerate Changed Plots*

All modifications done to plots mark the plots as "modified" which will display in a table titled "Edited plots to be regenerated." The "regenerate changed plots" button essentially remakes the plots with the new edited data. The default is just updating the plotly plots but if a user has the "Mass List and Parameters" excel sheet upload onto the app, it will also update the saved .png plots in the "plots" subfolder of each results folder. The stored plotly objects are compressed into a smaller file size which takes a few more seconds but decreases file size by an order of magnitude in cases with a large number of plots. This process does take time so it is advised that plots are regenerated after all edits for a specific data file have been performed. 

##### *Dual EIE View subtab*

This tab allows the user to display two electropherograms for use in identifying comigrating ions. Upon selecting this tab, the user will see two tables appear. These tables indicate to the program which plots you would like to compare. The table on the left will display the plot on top while the table on the right will display the bottom plot. Users will select the plots they wish to display by selecting the plot from the table. This will render a combined plot which can be zoomed into. Note: this can only be used to compare two plots from the **same** data file.


##### 5. **Downstream Processing**
***

This tab is where users can make large-scale modifications to their data after it has been edited or adjusted in the visualization tab. This requires you to have the 'Metabolite Peak Areas.csv' and the 'Metabolite Migration Times.csv' files. The user also must upload the 'Mass List and Parameters' file in the 'User Supplied Parameters' tab. Choose the existing run in Engine to identify its saved data, or upload real files for a new analysis. Visualization's selected analysis is also selected in Downstream Processing. The downstream results-folder control remains available to change the analysis used there. Once a results folder is set, the user can click the button 'Load Migration and Peak Area Data' which will grab the required .csv files from the results folder you selected and upload their information to the app. The peak area data can be viewed at any time in the Metadata subtab.


The app automatically generates a plot of the migration time versus the m/z and groups the data as internal standards and metabolites. The user can select the file they wish to look at in addition to selecting the peak being displayed.


As mentioned earlier, this tab is designed to allow the user to make large scale edits to their data. The app gives the user two different methods for handling missing data: replacing the missing data by the minimum peak area divided by five or replacing the missing data with zero. The user can select their method and then click the 'Replace Missing Data' button and all missing data in the peak area data will be replaced using the method you selected. This includes the <LOD, interfered, or NPD peak areas. This new data will then be saved to a new .csv file titled 'Corrected Missing Values Metabolite Peak Areas.csv' in the selected results folder. This will **not** modify the original 'Metabolite Peak Areas.csv', ensuring transparency along the data processing pipeline. 


The user can then normalize their data to one of the metabolites by clicking the 'Normalize Data' button. Upon clicking, a popup will appear will the user can select the metabolite they wish to normalize by. The current options for normalizing are using 144.0667_Creatinine, 184.0774_F-Phe, or 216.0427_Cl-Tyr. If the user wishes to normalize using a metabolite not listed here, they can simply add it to the server.R script in section 5.3 by modifying the following line:

```r
SelectInput("normalize_by", "Normalize By", choices = c("114.0667_Creatinine", "184.0774_F-Phe", "216.0427_Cl-Tyr"))
```

Once normalized, the new data will again be saved as a new .csv file titled 'Normalized Metabolite Peak Areas.csv' and the original peak area data will remain untouched.


Once the missing values have been replaced and the data has been normalized (normalization is not required, but encouraged), the user can connect existing metadata to their peak area results. This assigns a sample ID to the peak position. The name of the column for the sample ID is currently hard-coded to 'PBM Sample ID'. This can be edited in the server.R code as you see fit. This requires the user to create and place a .csv file into the 'Data' subfolder with the name 'metadata.csv'. **The name of this file is case sensitive and must be named metadata.** Once clicked, a new .csv file will be generated called 'merged_peak_areas.csv' The data stored on the app itself will reflect all these changes, but the original Metabolite Peak Areas file will remain unchanged. The tab 'Reporting' and the subtab 'Data Quality and Variability' require the metadata connected to the peak area results.    

##### *Control Charts*

This tab is for generating control charts for your data. This tab requires you to have loaded your migration and peak area data into the app in addition to having uploaded the 'Mass List and Parameters' file. Users can select a control chart they wish to display from the dropdown menu. Two control charts are generated when the user clicks the 'Create Control Chart' button. The first control chart shows the raw data from the 'Metabolite Peak Areas.csv' file. The second chart displays the same data after it has been normalized (if the user has normalized their data to a metabolite). 

##### *Data Quality and Variability*

This tab is for assessing the quality of the data runs through the use of CVs and a PCA. This is limited to comparing **two groups**: QC and non-QC. This requires the user to have connected the metadata to their peak areas. This tab allows the user to dynamically specify which samples are QC samples, which samples are not QC samples, and which samples to exclude from this analysis. To specify which samples are QC samples, the user can type 'QC=XYZ' in the text box titled 'Specify QC identifier' where 'XYZ' is PBM Sample ID for the quality control samples in your study. This entry cannot have any whitespace (no spaces), must be a single word; however, it is not case sensitive. PeakMeister2.0_UI additionally allows you to exclude samples from your variance assessments. Such samples would include things such as calibration solutions. These exclusions can be specified by typing the specific PBM sample IDs to exclude in the provided text box with each entry separated by a comma. For example, if the user wanted to exclude the sample IDs SP1, SP2, the blanks, and some calibration solutions, the user would type into the box 'SP1,SP2,Blank,C1,C2,C3,C4,C5,C6'. These exclusions **are** case sensitive.


Once these catagories have been specified, the user can click the 'Calculate CVs' button. This generates the following items. The first item generated is a PCA plot showing PC1 and PC2 for your QCs and non QCs. This PCA is also saved as a .png file with a quality of 300DPI in the 'data' subfolder located in the results folder you are working out of. The button also generates a table that lists the technical CV and biological CV as a percentage. The final item generated is a plot of the CV based on the m/z of the metabolite for both QC and non-QC samples, showing the median CV for each group.      

##### *Metadata*

This tab displays the peak area data you are working with after it has been uploaded onto the app.

##### 6. **Reporting**
***

This tab is for generating a datamatrix in the required format for [MetaboAnalyst](https://www.metaboanalyst.ca/). This requires the users to have connected the metadata to their peak area results. To process your data using the Generic format statistical analysis, MetaboAnalyst requires the data to contain a "CLASS" row, which is a way of categorizing your data into groups. This tab is designed to allow the user to dynamically assign classes to their data. In addition, it transforms and transposes your data into the required format to run your statistical analysis. Similar to the 'Data Quality and Variability' subtab, the syntax is 'PBM Sample ID = CLASS specifier' such as 'SP1 = Standard'. Class assignments can be both specific and generic. A specific assignment includes things such as 'SP1=Standard', 'SP2 = Standard', 'SRM = SRM3061'and '303 = Positive_Sample'. Generic class assignments would include something such as assigning all numeric values (typically your actual sample IDs) as 'Positive_Samples'. For example, to use a generic assignment the user would type 'Numeric=Positive_Sample'. Unassigned classes will be left blank. Class names are case sensitive (the value you are assigning a class to must be named what it is in the PBM Sample ID column) and multiple class specifications are separated by a comma. Class names must be one word. Whitespace is allowed when specifying multiple classes. An example of class assignments would be 'SP1 = Standard, SP2 = Standard, NIST = SRM3061, Numeric = Pos_Sample'. Upon clicking the 'Generate Data Matrix for MetaboAnalyst', the data will be transformed, a matrix will be made with your specified class assignments, and the matrix will be saved as a .csv file titled ' Metaboanalyst_matrix.csv' in the results folder you specified in the 'Downstream Processing' tab.   

### Copyright
***

PeakMeister is licensed under the [MIT](https://choosealicense.com/licenses/mit/) license

### References
***

