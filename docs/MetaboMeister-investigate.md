# Investigate raw data

Branch: `MetaboMeister`, based on `feature/modular-server` at `052bb2b`.
Replace both `server.R` and `ui.R` with the matching files from this branch.
All runtime functions remain in `server.R`; no additional package is required
beyond the existing mzR, Plotly and Shiny dependencies.

Investigate appears between Engine and Visualization. It reads the original
MS1 spectra independently of the analysis workflow. No calibration, smoothing,
peak detection or peak fitting is applied, and no result files are written.

## Open raw data

1. Upload one or more real `.mz5` files in Engine, then select the desired upload
   in Investigate. Alternatively, choose **Select another .mz5 file** and upload
   the raw file directly in Investigate without starting a new analysis.
2. Click **Load raw file**. This reads the scan metadata and checks up to eight
   evenly spaced MS1 spectra to initialize the m/z and time ranges. Missing or
   inconsistent header mass ranges are corrected from that sample; the status
   identifies the range source. It does not automatically draw the entire file.
3. Set the ranges and click **Draw / apply 3D ranges** or
   **Plot / apply cross-section**.

Opening a saved analysis in Engine does not load its original spectra. The full
m/z–time–intensity data are not available in the saved extracted traces; supply
the original `.mz5` in Investigate to explore that analysis's raw measurements.
Investigate's loaded-file status identifies the raw file being shown.

## Interactive 3D view

- Axes: m/z, acquisition time in minutes, and measured intensity.
- Drag to rotate. The Plotly toolbar provides pan, zoom, camera reset and PNG
  export; scroll zoom is enabled. **Reset 3D view** restores the camera and the
  last applied axis ranges without rereading the file.
- Set independent m/z, time and intensity display ranges. Intensity defaults
  to an automatic range; turn that off to supply explicit limits.
- Optional log-intensity display makes small signals easier to see alongside
  large signals. Negative intensities are shown at zero height in log view;
  linear view retains their values.
- Choose 100/200/400 m/z bins and 80/160/320 time bins. The preview contains at
  most 128,000 cells, even for a very large raw file. It reports the number of
  source measurements and MS1 scans represented.
- Each cell shows the maximum measured intensity in that m/z/time bin;
  empty bins are zero. The grid is a display summary rather than a peak-area
  measurement. Narrowing ranges increases the detail available per bin.
- Click the surface to populate **Centre m/z**. The orange band shows the
  currently selected m/z window. Click **Plot / apply cross-section** to extract
  the corresponding trace. A click uses the preview bin centre, not an exact
  measured mass. Enter an exact centre m/z for a narrow ppm window, or narrow
  the 3D display range first to obtain finer bins.

## Exact cross-section

The window is a centre m/z plus a **± half-width**, expressed in ppm or Da:

`half-width in Da = centre m/z × ppm / 1,000,000`

For example, 100 m/z ± 10 ppm includes values from 99.999 to 100.001 m/z.
Both endpoints are included.

The trace is calculated directly from the original MS1 scans, independently of
the 3D grid. Mean intensity is the default and follows PeakMeister's existing
window-aggregation convention. Sum and maximum intensity are also available.
An empty window contributes zero for that scan. Invalid/nonfinite data points
are excluded. This feature does not locate or integrate peaks.

The trace has its own time and intensity ranges. Drag to zoom, use the toolbar,
or select **Reset trace view** to restore the last applied ranges. Changing only
intensity limits reuses the cached data. Changing the window, extraction time
range or aggregation requires another extraction.

## Resource handling and errors

Spectra are read in batches of eight, rather than loading the full raw run into
an MSnbase object. Only one preview and one exact trace are cached per session.
File handles close when the loaded file is replaced, the investigation is
cleared, or the session ends. Opening a file that cannot be read and invalid
range requests use the existing guarded error handling, allowing retry.

## Code sections

Existing numbered section comments retain their original numbering and order.
New server code sits after 3.6 and before 3.7, under **Matabo3.6**, with subsections
for raw reading, extraction, GUI actions and plot outputs. The sidebar addition
uses **Matabo2.1** in `ui.R`; the page uses **Matabo5.1** after Engine and before
Visualization. The original analysis functions are unchanged.

`tests/investigate.R` checks exact extraction, MS1 selection, ppm/Da boundaries,
batch/grid limits, plots, caches, error recovery, handle cleanup and the actual
GUI controls. Existing action-guard, startup and scientific equivalence checks
also cover this branch.

The mzR integration check uses ProteoWizard's genuine HDF5/mz5 reference file,
`example_data/tiny.pwiz.1.1.mz5`, from commit
`c9a8cadc630568b519e1fee00ce2aec71c28d695`. The reader and exact aggregates are
checked independently of the preview. Chromium checks cover actual WebGL
rotation/reset, linked selection, ranges, tab switching and raw-file upload.
Large laboratory files were not available for a performance benchmark.
