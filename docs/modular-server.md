# Modular server

Branch: `feature/modular-server`, based on `feature/new-or-open-analysis`.
All application functions remain in `server.R`. The existing `ui.R` is compatible.

The Initialize button validates the selected uploads and calls
`run_new_analysis()`. That function prepares the parameters and calls
`process_mz5_file()` for each uploaded file. The per-file function calls the
processing tasks in order and records their existing timing stages.

| Task | Function(s) |
| --- | --- |
| Validate parameters and calculate/override migration indexes | `prepare_analysis_parameters()` |
| Read the temporary raw file | `read_mz5_data()` |
| Optional mass calibration | `calibrate_masses()` |
| Prepare extraction windows and sort spectra | `prepare_extraction_spectra()` |
| Extract traces and construct `eie_df` | `extract_electropherograms()`, `build_eie_dataframe()` |
| Optional smoothing | `smooth_electropherograms()` |
| Detect, integrate and filter internal standards | `process_internal_standard_peaks()` |
| Detect, integrate and align metabolites | `process_metabolite_peaks()` |
| Flag low/undetected and interfered peaks | `flag_peaks_below_lod()`, `flag_interfered_peaks()` |
| Prepare annotations, integration ribbons and axes | `build_plot_data()` |
| Save PNGs and build interactive plots | `save_static_plots()`, `build_interactive_plots()` |
| Prepare report tables | `build_result_tables()` |
| Save RData, SQLite, CSVs and Plotly objects | `save_analysis_snapshot()`, `save_analysis_database()`, `export_result_csvs()`, `save_interactive_plots()` |
| Stage timing and error/interrupt reports | `new_run_timer()`, `run_new_analysis()` |
| Delete selected/middle/all peaks and undo | `delete_selected_peaks()`, `delete_middle_peaks()`, `delete_all_peaks()`, `undo_peak_deletion()` |
| Regenerate changed plots and mark interferences | `regenerate_changed_plots()`, `mark_peak_interferences()` |
| Normalize and replace missing data | `normalize_peak_areas()`, `replace_missing_peak_areas()` |
| Join metadata to sample injections | `connect_peak_area_metadata()` |
| Control charts, CVs and PCA | `build_control_chart()`, `compute_quality_metrics()`, existing `calculate_cv()` and `perform_pca()` |
| Batch correction | `load_batch_correction_packages()`, `apply_batch_correction()` |
| MetaboAnalyst matrix | `build_metaboanalyst_matrix()` |
| Drug databases | `load_drug_databases()`, `parse_and_save_drugbank()` |

Data-processing helpers accept tables/objects explicitly and return their results.
GUI action helpers use their own Shiny session's reactive state. The observers
remain responsible for reacting to buttons and updating the GUI. Existing saved
analysis discovery, SQLite helpers, manual integration and baseline functions
remain in place.

The scientific formulas, filters, exported column names, saved-object names,
folder layout and GUI controls retain their previous behavior. A redundant
assignment of the last plot's already-built data was removed. Temporary raw
copies and pipeline database writer connections now also close on errors.

Validation:

- `Rscript --vanilla tests/gui-startup.R`
- `Rscript --vanilla tests/saved-run-selection.R`
- `Rscript --vanilla tests/modular-server.R`
- `Rscript --vanilla tests/modular-equivalence.R /path/to/previous/server.R`

The processing tests use synthetic spectra/traces and a trapezoid-integration
fixture. They exercise the actual peak algorithms, CSV/RData/SQLite exports,
Plotly traces/layouts, both PNG naming layouts, multiple-file processing, timing
and failure cleanup. PNG writing and raw MSnbase decoding use test adapters.
The equivalence check compares those outputs with the pre-refactor server
(commit `3604a346c92089342519746d578ecb19f0bf4488`) and compares the moved GUI
action bodies. GUI startup and saved-run tests use the actual handlers.

A complete analysis of a real `.mz5` file was not run in this environment;
MSnbase/mzR and the full scientific dependency set are unavailable here.
