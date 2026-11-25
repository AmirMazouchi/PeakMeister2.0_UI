# Investigate shared Engine data

Use matching `server.R` and `ui.R` from branch **MetaboMeister**. Runtime helpers
remain self-contained in these files and use the existing application packages.

## Engine owns the upload

Upload `.mz5` files in Engine. Each upload is read once into memory. Investigate
automatically selects that file and populates its ranges, sampling counts,
controls and status before fine-cache calculation starts. The upload handler first
flushes the filename and a preparing status to the browser, then queues Engine's
read and sampling work. Preparation reports progress after `File Reading Complete`.
Cache preparation also runs after the UI update; its failure leaves the Engine
filename and spectra available. A reader failure displays the filename and an explicit error status.
There is no Investigate upload or raw-file-open button.
Selecting another already loaded Engine file reuses its memory snapshot.

Engine retains the original uncalibrated object. Investigate references its Spectrum
objects and materializes only bounded scan batches; it does not duplicate the full
run into numeric matrices. Source fingerprints are calculated per scan during
sampling, avoiding a second whole-run serialization. Analysis receives a deep copy
of its mutable assay-data environment, so mass calibration and extraction cannot
change Investigate's original spectra. Cross-sections, 3D previews and Median
read the in-memory snapshot, not the original file again.

**Start run automatically after upload** is enabled by default. Turn it off in
Engine if only inspection is wanted. The 3D time-start default is **1.2 minutes**,
also after loading a file; the recorded end populates its upper range. Grouping
defaults remain 40 mass readings and 20 scans. The 3D preview uses grouped means;
cross-sections use the original scans. Existing rotation, zoom and resets remain.

## Fine Maximum and Minimum caches

**Maximum** is the default m/z statistic. Both extrema are calculated together
from all valid MS1 scans within the selected calculation limits. The defaults
cover the recorded run, independently of the 3D grouping and Median controls.

All scans and slices use one global grid: equal **0.001 Da** cells anchored at
integer multiples of 0.001 Da, with intensities plotted at cell centers. This
grid is used even when scans share identical native mass coordinates. Multiple
readings in a scan's cell are averaged before aggregation. Only internal empty
cells between adjacent occupied nonzero cells are linearly interpolated. Measured
zeros, outer gaps and empty scans retain the original zero semantics.

New caches use compressed blocks and 10-second time slices, with partial slices
at the recorded start and end. Intervals include the start and exclude the end,
except at the recorded end, which is included. A boundary scan belongs to one
slice only, including duplicated acquisition times. Empty time slices have zero
scan count and no summaries; they never add false zero observations. A window
containing no scans displays a gap and reports zero scans.

Combining slices takes the element-wise **maximum of slice maxima** or
**minimum of slice minima**, never an average. These summaries are separate from
Engine's scientific analysis, which retains its existing calculations.

The paired caches are stored alongside dated analysis folders:

```text
Results/<mz5 basename>/
  max_mz_spec_dmz0.001.sqlite
  min_mz_spec_dmz0.001.sqlite
  Results <date>/
```

To recover existing completed files, use Investigate → **Saved spectra**, choose
its run, then click **Load / plot saved spectra**. No mz5 upload or recalculation
is required. Select Maximum or Minimum to switch the curve. Keep the two SQLite
files in their existing `Results/<mz5 basename>/` directory. Validation failures
are reported explicitly and do not replace the files.

Each database stores indexed `(slice, mass-block)` records with binary double
intensities in blocks of 4,096 cells. Queries fetch at most 64 blocks at a time
and compute each block's offsets directly, avoiding both whole-cache allocation
and repeated full-axis scans. Initial plotting shows progress. A selected region such as 162–165 m/z reads
only its intersecting blocks and selected nonempty slices, plus line-continuity
neighbors. It does not combine the entire mass range while zoomed in.

Metadata records the original spectra/time fingerprint, filename, common grid,
slice boundaries/counts, interpolation and alignment rules, schema version,
completion flag and generation. Both databases must match. Incomplete,
outdated, structurally damaged or mismatched pairs are rejected. A matching
Engine upload reuses completed caches without traversing its scans again.
Different raw data with the same filename invalidates the old pair.

Reopening a saved run in Engine automatically loads its matching completed
caches without an mz5 read. It validates the stored pair/settings; a subsequent
raw upload additionally validates its spectra/time fingerprint. Raw-only views
(3D, cross-sections and Median) still require the original upload through Engine;
the fine extrema caches do not contain every original spectrum.

## Window controls and plotting

For Maximum and Minimum, enter **Calculation start (min)** and
**Calculation end (min)** and select **Apply time limits**. **Full** uses this
entire interval. The other widths are 20 sec, 30 sec, 1 min, 2 min, 3 min,
4 min, 5 min and 10 min. The position slider moves the window within the
selected interval. If the selected interval is shorter than the width, the
plot uses the entire interval. The display shows the actual endpoints.

With the original mz5 loaded through Engine, the app reads any scans needed
from partially covered cache slices to honor the entered endpoints. With only
saved spectra loaded, choose saved time boundaries. Older caches remain
readable with their original boundaries. Reopen the original mz5 through
Engine to generate the finer cache and use the short widths.

Dragging sends positions live, including before release. Region summaries update
only the SVG curve and hover label. The plot widget, m/z zoom and intensity-axis
limits remain unchanged. Zooming into a new region queries its cached mass blocks.
Reset restores the full mass/intensity view. Broad views retain per-interval
minima and maxima at their original coordinates, preserving narrow peaks within
the bounded overview; detailed zoom restores the cached points.

## Median

Median remains available with independent time-start/end, every-Nth-scan and
alignment-width controls. Click **Calculate Median** after editing them. The
existing exact sparse/partitioned calculation, interpolation, native-axis behavior
and scan-stride semantics remain in place. It uses Engine's uncalibrated snapshot.
The fine 0.001 Da slice workflow does not apply to Median.

## Validation

- `tests/investigate.R`: existing numerical, display, 3D, cross-section and Median
  checks plus the new shared/cache tests and actual UI defaults.
- `tests/investigate-cache.R`: independent extrema oracle, global grid, duplicate
  boundary times, final/empty slices, region reads, persistence, source/settings
  validation, incomplete-cache rejection, saved-run reuse and actual GUI handlers
  retaining both axes during live curve updates.
- `tests/shared-engine-mz5.R <fixture.mz5>`: real MSnbase data, production analysis
  reuse without a second read, and calibration isolation from original spectra.
- `tests/shared-analysis-equivalence.R <pre-change server.R>`: two-file synthetic
  comparison of CSVs, saved objects, SQLite tables, plots and pipeline stages.
- `tests/investigate-mz5.R <fixture.mz5>`: genuine mz5 reading and independently
  checked trace/spectrum extraction. The fixture used is ProteoWizard's
  `example_data/tiny.pwiz.1.1.mz5`.

The optional `tests/investigate-large-cache.R` scale check creates two approximately
637 MiB caches, verifies independently calculated full-run Maximum/Minimum and a
narrow m/z region, then reopens the completed files without raw spectra. The cloud
check used about 345 MiB peak resident memory; full queries took approximately
9–13 seconds, and the narrow query about 0.07 seconds. Times depend on hardware.

A Chromium check also recovers completed caches in a fresh session with no mz5
upload. It uploads through Engine and drags the actual slider while zoomed
to 162–165 m/z and fixed intensity limits. It verifies updates before release,
unchanged axes and the same SVG widget. Large laboratory files, full user-data
analysis and RStudio's embedded viewer were not available for validation.

## Engine mass calibration workbook

Each new Engine analysis writes `Mass Calibration - <dated analysis folder>.xlsx`
in `Results/<mz5 basename>/` after the calibration step. The filename includes
the dated folder name and any rerun suffix, so earlier workbooks are preserved.
Opening a saved analysis does not regenerate its workbook.

- **Mass List** and **Reference Mass List** include every workbook target. The
  input m/z is the value used for extraction. The modeled corrected m/z columns
  show the minimum, mean and maximum of the scan-specific correction evaluated
  at that input m/z across calibrated scans. These are model estimates.
- **Calibration References** lists both expected lock masses from the Parameters
  sheet and their measured m/z before calibration in scan processing order.
  A missing lock pair skips that scan. A lock mass measured on a skipped scan
  can still have a before value.
- **Summary** identifies the analysis and records whether calibration was applied,
  disabled, or enabled with no usable scan pairs. Unavailable measurements and
  unapplied corrections appear as blank cells.

Run the workbook check with `Rscript --vanilla tests/mass-calibration-report.R`.
