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
2. Click **Load raw file**. This reads scan metadata and makes one streamed pass
   over all valid MS1 spectra to measure native sampling, check whether m/z
   coordinates are shared, and initialize the observed ranges. The progress bar
   identifies this pass; it adds reading time on load. These measurements are
   cached until the file is cleared/replaced. Plots are drawn only on request.
3. Set the ranges and click **Draw / apply 3D ranges** or
   **Plot / apply cross-section**.

Opening a saved analysis in Engine does not load its original spectra. The full
m/z–time–intensity data are not available in the saved extracted traces; supply
the original `.mz5` in Investigate to explore that analysis's raw measurements.
Investigate's loaded-file status identifies the raw file being shown.

## Native sampling information

The loaded-file summary reports:

- MS1 scan/time-slice count and number of distinct acquisition times.
- Minimum positive time difference between adjacent distinct MS1 times, in
  minutes and seconds. Duplicate timestamps do not give a zero spacing.
- Minimum positive difference between adjacent distinct m/z positions **within
  an individual spectrum**, measured across all valid MS1 spectra. Differences
  between slightly shifted masses in different scans are not used.
- The number of shared native m/z positions when nonempty spectra have identical
  coordinates. Empty spectra still count as time slices and contribute zeros.
- If m/z coordinates vary, the minimum–maximum number of distinct m/z positions
  per spectrum, including empty spectra. Such data have no single shared native
  m/z axis, so a universal native m/z-slice count would be misleading.

These are measured **sampling spacings**, not an estimate of the instrument's
mass resolving power or peak widths. A spacing is reported as unavailable when
there are not two distinct readings to compare. Invalid/nonfinite measurements
are excluded. The bounds come from the actual valid MS1 data, including files
with missing or inconsistent header ranges.

When coordinates differ between scans, an adjustable **Shared m/z alignment
width (Da)** appears. Each scan's readings in a grid cell are averaged first;
missing cells contribute zero. Its initial width is at least the measured
minimum spacing and gives at most about 200,000 initial positions. This is a
chosen alignment grid, labelled separately from measured native spacing; it is
not interpolation or a claim about instrument resolution. Its last cell can
be narrower at the requested upper bound. It starts at the selected m/z lower
limit, so changing that limit can change alignment. Narrow the m/z range or
increase the width if a request exceeds 2,000,000 alignment positions.

Shared native coordinates are used directly, including nonuniform spacing;
the alignment control is hidden and its value has no effect in that case.

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
- Set positive whole **m/z grouping** and **time grouping** factors. Factor
  **1** keeps each position/scan; **10** combines each consecutive group of 10.
  Time groups follow scan acquisition order rather than equal-duration bins.
  For varying m/z coordinates, mass grouping follows the aligned positions.
  The last group uses its actual size if it contains fewer than the factor.
- Loading initializes factors to give roughly **50 × 40 display positions**
  (or fewer for a small file), rather than fixed factors that could create a
  huge plot. Both factors remain editable. The live summary shows source and
  resulting position counts, while the drawn plot reports its applied factors.
  A preview is limited to 1,000,000 cells; increase factors or narrow ranges
  for larger requests. A single resulting row/column uses 3D lines/markers
  because a surface requires at least two positions on each axis.
- The surface uses a light blue–cyan–yellow–orange–red palette with brighter
  lighting, avoiding the previous dark purple low-intensity background.
- Each cell shows **mean intensity** across its consecutive mass positions
  and scans, including zeros for missing positions. If several raw readings
  share one native/aligned mass position, they are averaged per scan first.
  Group centres are means of the represented mass coordinates/scan times;
  groups need not have equal physical width when native spacing is irregular.
  This is a display summary, not a peak-area measurement.
- Click the surface to populate **Centre m/z** or **Slice time**, depending on
  the active slice axis. The orange marker shows the selected m/z or time
  window. Click **Plot / apply cross-section** to extract it. A click uses a
  group centre, which can differ from a recorded mass or scan time. Enter an
  exact value for precise windows, or reduce grouping factors for finer detail.

## m/z-window cross-section

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

## Time-window cross-sections

Choose **Time → mass spectrum**, enter a centre time in minutes and a positive
± window in seconds. The lower graph then plots intensity against m/z, with
separate spectrum m/z-range controls and the existing intensity-range/reset
controls. Only MS1 scans inside the inclusive time window are used.

Two displays are available:

1. **Combine spectra (mean):** align scans to an explicit m/z grid (default
   bin width 0.01 Da). Compute each scan's mean intensity within each mass bin,
   then average those bin intensities across all selected scans. A missing bin
   in a scan contributes zero, so a scan with more raw points does not receive
   extra weight. The status reports the number of scans and bin width. The last
   mass bin can be narrower where it meets the requested m/z upper limit.
2. **Overlay spectra:** choose the number of spectra, default 5 (4 is also
   supported). Select up to that many closest scans inside the time window,
   ordered by absolute distance from its centre; ties use the earlier time.
   Each curve keeps its original m/z points and has a different colour. Legend
   entries give actual acquisition times and time distances, nearest first.
   Click a legend entry to hide/show a curve. Empty spectra in the chosen m/z
   range are omitted and the status reports the selected and plotted counts.

Combined spectra are averages on their own stated mass grid; they do not use
the coarse 3D preview. Overlaid spectra are not mass-binned. Both are calculated
from the original data with no smoothing, calibration or peak fitting. A
combined spectrum is limited to 1,000,000 mass bins; narrow the m/z range or
increase its bin width for a larger request.

## Median m/z spectrum

Click **Calculate median spectrum** for a separate m/z-versus-intensity graph.
It uses every MS1 scan inside the **current 3D time range controls**, including
both endpoints, and the current 3D m/z limits. At each shared native m/z position
(or aligned position if axes differ), there is one intensity per scan; missing
signals and empty scans contribute zero. Multiple readings at the same
position are averaged within that scan before taking the median across scans.
No 3D grouping is applied to this graph, and changing only grouping factors
does not invalidate its cache. A range/alignment change displays a reminder to
recalculate; the status always identifies the plotted time range and scan count.

The median is exact, including the mean of the two middle intensities for an
even number of scans. Short-lived peaks can have **zero median** when absent
from most scans; the status explains an all-zero result. This is expected for
many transient signals and does not mean the original data have disappeared.
Use the time-window combine/overlay views to inspect signals in shorter windows.
The graph supports zoom, pan, PNG export and **Reset median view**.

The implementation avoids allocating a dense scans-by-m/z matrix. It counts
positive/negative observations with implicit zeros first, skipping value sorting
where those counts already prove a zero median. Nonzero values are cached for
small requests, and only the central ranks are selected rather than fully
sorting every position. For larger requests it reads the selected scans once
more and partitions only candidate nonzero values into temporary binary files.
Partitions are processed individually with a limited number of open files;
they are deleted on success or error. No new package or compiled extension is
required. Repeating the same median request reuses the result. Initial reading
and alignment are still proportional to the selected raw data size; the feature
does not promise instant processing of large files.

## Resource handling and errors

Spectra are read in batches of eight, rather than loading the full raw run into
an MSnbase object. Native sampling information, one preview, the current
cross-section and one median result are cached per session.
File handles close when the loaded file is replaced, the investigation is
cleared, or the session ends. Opening a file that cannot be read and invalid
range requests use the existing guarded error handling, allowing retry.
Raw loading and calculation actions explicitly call `shiny::withProgress()` and
`shiny::setProgress()` so other loaded packages cannot mask progress updates.

## Code sections

Existing numbered section comments retain their original numbering and order.
New server code sits after 3.6 and before 3.7, under **Matabo3.6**, with subsections
for raw reading, extraction, GUI actions and plot outputs. The sidebar addition
uses **Matabo2.1** in `ui.R`; the page uses **Matabo5.1** after Engine and before
Visualization. The original analysis functions are unchanged.

`tests/investigate.R` checks native sampling, fixed/varying mass axes, adjacent
grouped means, independently calculated exact medians (odd/even, signed, sparse,
RAM and disk paths), MS1/ppm/Da boundaries, batch/grid limits, plots, caches,
error recovery, temporary-file/handle cleanup and the actual GUI controls.
Existing action-guard, startup and scientific equivalence checks also cover
this branch.

The mzR integration check uses ProteoWizard's genuine HDF5/mz5 reference file,
`example_data/tiny.pwiz.1.1.mz5`, from commit
`c9a8cadc630568b519e1fee00ce2aec71c28d695`. The reader and exact aggregates are
checked independently of the preview. Chromium checks cover actual WebGL
rotation/reset, linked selection, grouping factors, independent time-range
medians, ranges, tab switching and raw-file upload.
Large laboratory files were not available for a performance benchmark.
