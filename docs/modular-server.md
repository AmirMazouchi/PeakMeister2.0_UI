# Modular server

Branch: `feature/modular-server`, based on `feature/new-or-open-analysis`.
All application functions remain in `server.R`. The existing `ui.R` is compatible.

Visualization's unused `.RData` upload box and **Clear uploaded data** button
have been removed. To reopen another data file, select its run folder in
**Engine → Open existing analysis**, then choose the saved analysis version
in **Visualization**. This sequence can be repeated to switch runs and versions.

The Initialize button validates the selected uploads and calls
`run_new_analysis()`. That function prepares the parameters and calls
`process_mz5_file()` for each uploaded file. The per-file function calls the
processing tasks in order and records their existing timing stages.

The UI page/trigger column identifies where the user starts each task. Helpers
are generally called indirectly through server functions, rather than directly
by `ui.R`. Engine tasks also run when **Start run automatically after upload**
is enabled. Their parameter tables come from **User Supplied Parameters**;
processing starts on **Engine**.

| Task | Function(s) | UI page / trigger |
| --- | --- | --- |
| Catch action/update errors and keep the session open | `observe_action()`, `observe_updates()`, `with_action_error()`, `report_action_error()`, `action_label()` | All pages → existing button, upload, selection and background-update handlers |
| Coordinate uploaded files and process each file | `run_new_analysis()`, `process_mz5_file()` | Engine → New analysis → Initialize Run |
| Validate parameters and calculate/override migration indexes | `prepare_analysis_parameters()` | Engine → New analysis → Initialize Run |
| Read the temporary raw file | `read_mz5_data()` | Engine → New analysis → Initialize Run |
| Optional mass calibration | `calibrate_masses()` | Engine → New analysis → Initialize Run |
| Prepare extraction windows and sort spectra | `prepare_extraction_spectra()` | Engine → New analysis → Initialize Run |
| Extract traces and construct `eie_df` | `extract_electropherograms()`, `build_eie_dataframe()` | Engine → New analysis → Initialize Run |
| Optional smoothing | `smooth_electropherograms()` | Engine → New analysis → Initialize Run |
| Detect, integrate and filter internal standards | `process_internal_standard_peaks()` | Engine → New analysis → Initialize Run |
| Detect, integrate and align metabolites | `process_metabolite_peaks()` | Engine → New analysis → Initialize Run |
| Flag low/undetected and interfered peaks | `flag_peaks_below_lod()`, `flag_interfered_peaks()` | Engine → New analysis → Initialize Run |
| Prepare annotations, integration ribbons and axes | `build_plot_data()` | Engine → New analysis → Initialize Run |
| Save PNGs and build interactive plots | `save_static_plots()`, `build_interactive_plots()` | Engine → New analysis → Initialize Run |
| Prepare report tables | `build_result_tables()` | Engine → New analysis → Initialize Run |
| Save RData, SQLite, CSVs and Plotly objects | `save_analysis_snapshot()`, `save_analysis_database()`, `export_result_csvs()`, `save_interactive_plots()` | Engine → New analysis → Initialize Run |
| Stage timing and error/interrupt reports | `new_run_timer()` | Engine → New analysis → Initialize Run |
| Delete selected/middle/all peaks and undo | `delete_selected_peaks()`, `delete_middle_peaks()`, `delete_all_peaks()`, `undo_peak_deletion()` | Visualization → Plotting → Delete Selected Peaks / Delete Middle Peaks / Delete all peaks / Undo Peak Deletion, respectively |
| Regenerate changed plots | `regenerate_changed_plots()` | Visualization → Plotting → Other → Regenerate Changed Plots |
| Mark peak interferences | `mark_peak_interferences()` | Visualization → Plotting → Mark Peak Interference |
| Normalize and replace missing data | `normalize_peak_areas()`, `replace_missing_peak_areas()` | Downstream Processing → left controls → Normalize Data (then Normalize in the dialog) / Replace Missing Data, respectively |
| Locate shared metadata and join it to sample injections | `metadata_file_paths()`, `connect_peak_area_metadata()` | Downstream Processing → left controls → Connect Metadata to Results; merged table appears in Metadata |
| Control charts | `build_control_chart()` | Downstream Processing → left controls → Create Control Chart; plots appear in Control Charts |
| CVs and PCA | `compute_quality_metrics()`, `calculate_cv()`, `perform_pca()`, `build_quality_cv_plot()`, `build_quality_pca_plot()` | Downstream Processing → Data Quality & Variability → Calculate CVs |
| Batch correction | `load_batch_correction_packages()`, `apply_batch_correction()` | Downstream Processing → Batch Correction → Initialize Packages / Perform Batch Correction, respectively |
| MetaboAnalyst matrix | `build_metaboanalyst_matrix()` | Reporting → Generate Data Matrix for MetaboAnalyst |
| Drug databases | `load_drug_databases()`, `parse_and_save_drugbank()` | Drug Analysis → Load Drug Databases; parsing starts with Yes, parse now in the DrugBank Database Not Found dialog |

Data-processing helpers accept tables/objects explicitly and return their results.
GUI action helpers use their own Shiny session's reactive state. The observers
remain responsible for reacting to buttons and updating the GUI. Existing saved
analysis discovery, SQLite helpers, manual integration and baseline functions
remain in place.

Metadata is shared between dated analyses at `Results/<run>/metadata.csv`.
`metadata_file_paths()` derives that path from the selected analysis, so switching
analyses within a run reuses the same metadata while switching runs uses a
different file. The shared file takes priority over the older
`<analysis>/Data/metadata.csv` location, which remains a fallback. Older flat
results layouts keep their existing analysis-local lookup. Missing-file messages
show the expected location. No metadata file is moved or overwritten automatically;
merged output remains `<selected analysis>/merged_peak_area.csv`.

The extraction/peak-processing formulas, filters, exported column names, saved-object names,
folder layout and GUI controls retain their previous behavior. A redundant
assignment of the last plot's already-built data was removed. Temporary raw
copies and pipeline database writer connections now also close on errors.

The Data Quality & Variability calculation now validates metadata sample IDs,
the QC identifier, group sizes (at least two QC and two non-QC samples), and
numeric peak areas. Entirely unusable metabolite columns are excluded with a
message; remaining missing/nonnumeric/infinite values must be corrected before
calculation. Undefined zero-mean CVs are reported as missing rather than infinite.
CV plots retain the correct metabolite/value association after exclusions and
sorting by m/z.

PCA now uses base R's `stats::prcomp()` SVD instead of the failing
FactoMineR/irlba path. The previous centering and population-standard-deviation
scaling after Pareto scaling are preserved. The returned scores and eigenvalue
table match the fields used by the GUI. Too few varying metabolites, insufficient
rank or a PCA error produce a message while retaining the CV results. Confidence
ellipses are drawn only for groups with at least four points and full rank in
the two displayed dimensions. Invalid input clears previous results; the button
handler catches errors without ending the session.
PCA panel checks explicitly call `shiny::validate()` and `shiny::need()` so
other packages' validators cannot mask Shiny's display validation.

All existing event handlers now register through `observe_action()` and plain
background observers through `observe_updates()`. These call a shared error
guard, which stops a failed action, displays a persistent dismissible notification
with the action name and error message, and logs the condition/call to the R
console. Loading waiters are released after errors; a failed database migration
also dismisses its non-closeable progress modal. Missing-input `req()` conditions
and Shiny validation retain their normal silent/inline behavior. Handler options
such as priority, ignoreInit, ignoreNULL and once are forwarded to Shiny.

The parameter-workbook download content is guarded too: a failed write removes
the partial temporary download file. The MetaboAnalyst save-error message now
formats the original I/O error correctly. Plot/table outputs continue to use
Shiny's existing inline error handling.

These guards cover ordinary R errors in actions/updates; they do not undo
reactive-state changes or file/database writes completed before an error. Existing
pipeline error timing and resource cleanup still run before the error reaches
the UI guard. Startup failures and termination of the R process are outside this
session-level protection.

Validation:

- `Rscript --vanilla tests/gui-startup.R`
- `Rscript --vanilla tests/saved-run-selection.R`
- `Rscript --vanilla tests/modular-server.R`
- `Rscript --vanilla tests/quality-metrics.R`
- `Rscript --vanilla tests/action-errors.R`
- `Rscript --vanilla tests/metadata-location.R`
- `Rscript --vanilla tests/modular-equivalence.R /path/to/previous/server.R`

The processing tests use synthetic spectra/traces and a trapezoid-integration
fixture. They exercise the actual peak algorithms, CSV/RData/SQLite exports,
Plotly traces/layouts, both PNG naming layouts, multiple-file processing, timing
and failure cleanup. PNG writing and raw MSnbase decoding use test adapters.
The equivalence check compares those outputs with the pre-refactor server
(commit `3604a346c92089342519746d578ecb19f0bf4488`) and compares the moved GUI
action bodies. GUI startup and saved-run tests use the actual handlers.

The quality-metrics tests use a 12-sample, 20-metabolite fixture to verify PCA
against an independent correlation-matrix eigendecomposition. They exercise
the actual Calculate CVs observer, including the reported NULL-object failure,
continued CV display, invalid input, and recovery on the next valid calculation.
The observer is also tested with `jsonlite` loaded after Shiny and an S3
validator that reproduces the reported `validate(NULL)` method error.
The native PCA scores and eigenvalues were also compared with FactoMineR's
full-SVD reference source, independently of irlba.

The action-error tests audit observer registrations and active button coverage,
then exercise real metadata, batch, workbook upload, integration, raw-run timing,
and download callback failures. They verify that the session stays open,
notifications identify the action, loading waiters are hidden, failed downloads
are removed, and corrected input can be retried. They also verify silent
validation, early returns, observer isolation and once behavior.

The metadata-location tests use the actual load/connect observers to check
shared-file reuse across versions, precedence over local files, switching runs,
legacy fallback, missing-file messages, and analysis-specific merged exports.

A complete analysis of a real `.mz5` file was not run in this environment;
MSnbase/mzR and the full scientific dependency set are unavailable here.
