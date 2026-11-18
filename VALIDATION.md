# Release validation: 0.1.0

Checked on 18 November 2026.

- 69 numerical, import, spectral reconstruction, and CLI tests passed. These include independent synthetic
  lifetime recovery, sparse/dense equivalence, zero photon counts, all seven
  contributed multi-curve exports, original time-label retention, and explicit
  source-curve/wavelength handling.
- Headless Chromium checks editing actual manual-table cells, applying values
  and downloading the same normalized intensities and manual-source metadata.
  Manual entry completed without JavaScript errors.
- Manual-table AppTest checks verify fixed wavelengths/rows, blank initial cells,
  missing/negative/above-one/nonfinite/all-zero rejection, exact preservation
  of entered normalized scale, the independent reconstruction formula, unit
  TRANES wavelength area, export values/provenance and preview-rerun persistence.
- Spectral workflow checks use independently generated and globally fitted decay
  curves to verify that wavelengths, component amplitudes and lifetimes connect
  directly to a separate normalized steady-state file. Shuffled wavelength rows
  retain the correct association; a normalized subset whose maximum is below
  one is used unchanged. Missing, extra, duplicate and ambiguous rows are
  rejected in the exact matching path.
- Full-spectrum tests verify that division by the source maximum and min-max
  normalization both occur before sampling at the fitted wavelengths. Source
  extrema outside the fitted wavelength subset are included. Signed raw
  baselines are accepted only with explicitly selected min-max normalization.
- PCHIP and B-spline curves pass through the calculated samples. PCHIP preserves
  nonnegative bounds; a known B-spline negative overshoot is retained as a
  display artifact. Changing connections leaves all TRES/TRANES matrix values
  and their sampled-grid spectral areas unchanged. A CSV round-trip check
  verifies that normalized and matched steady-state exports can be loaded again
  and reproduce the same TRES/TRANES values.
- Single-decay AppTest checks exercise the shared kinetics preview, handles
  above the curve, bidirectional channel/time snapping, whole-record and
  selected-window views, and disabling fits with fewer than four channels.
  Selecting 5–30 ns in the bundled single decay retains all 1,601 input bins
  while the objective uses channels 200–1,200 inclusive (1,001 bins).
- Tail-mode checks verify a separate post-peak starting interval and an uploaded
  decay with a different grid, fitted without an IRF. All 601 imported bins are
  retained and the selected objective uses 321. Global preview/trace selection
  checks retain the full 600-by-7 matrix with 451 selected bins for 2–20 ns, and
  the Continue button still opens the saved global result in TRES/TRANES.
- The supplied reconstruction notebook was read as reference code. Its widget
  callbacks and file-writing operations were not executed. Independent scalar
  calculations using the supplied three-component table reproduce the TRES
  engine within 5.6e-17 absolute intensity, and TRANES within 6.1e-18. Numerical
  time integration recovers the supplied steady-state values; trapezoidal
  integration of each TRANES spectrum over its actual, uneven wavelength grid
  gives one to numerical precision. See the
  [notebook review](docs/notebook_review/NOTEBOOK_REVIEW.md) and its adjacent
  reference matrices for the calculation, interpolation checks, and assumptions.
- Reconstruction checks cover dense steady-state interpolation, coverage
  rejection, normalization, wavelength/column alignment, tail-reference
  conversion, zero-area rejection, per-wavelength amplitude-scale invariance,
  and exported normalization metadata. The supplied lifetimes are reference
  inputs, not independently established ground truth.
- Streamlit 1.67.0 AppTest exercised the bundled single-decay fit and a full
  51-trace global fit, including residuals, autocorrelation, amplitude spectra,
  import metadata and result download rendering.
- Spectral-workflow AppTest checks exercise an actual two-component, seven-wavelength global
  fit, the Continue button, automatic parameter handoff, exact normalized
  steady-state upload, full-spectrum peak and min-max paths, signed raw input,
  display-style invariance, spectral ZIP export and clearing of invalid results.
  A subsequent single-decay fit preserves the earlier saved global result.
  Synthetic steady-state inputs are used only for UI tests; no steady-state
  spectrum is assigned to the separate measured demonstration dataset.
- Linked fit-range controls were checked in both directions, including channel
  snapping, data/method reset, invalid short windows, individual-trace preview,
  and zoom. For a selected 6–30 ns window, all 1,999 imported bins and 51 traces
  remain available to reconvolution; the objective uses exactly channels
  240–1,200 inclusive (961 channels).
- Headless Chromium checks the single-decay preview, numerical fields, native
  range-handle updates, whole/selected views and an actual 5–30 ns fit. It verifies
  the displayed 1,001-channel fit interval and completes without JavaScript
  errors. Screenshots show the settled full-record and selected-window views.
  The global/reconstruction route is covered by the interface checks described above.
  Runtime logs contain no deprecated width-argument notices.
- Official Streamlit 1.51.0 API signatures were checked for the width arguments
  used by figures, dataframes, buttons and downloads. Runtime tests used 1.67.0.
- Release text/code, importer metadata, CLI help and example output were scanned
  for removed product branding.
- The experimental user guide covers existing-environment startup, first
  installation, single/global fitting, range selection, result interpretation,
  connected spectral reconstruction, exports and troubleshooting.
- Beginner guide commands, environment names, file paths and launchers were
  cross-checked against the release. Windows command files are supplied but were
  not executed on Windows in this Linux validation environment. Whether a local
  application-control policy accepts installed libraries must be checked there.

Reproduced example lifetimes (ns):

| Example | Result |
|---|---|
| Single decay, one component | 2.549156 |
| Synthetic global, truth 0.7 and 3.3 | 0.704709 and 3.300522 |
| Contributed measured global, two components | 2.333760 and 3.882271 |

The measured global reduced Poisson deviance is approximately 1.381. This indicates
remaining mismatch under the selected shared-IRF exponential model. The result is
an illustrative model fit, not ground truth; local uncertainties exclude model
and calibration systematics. See SCIENTIFIC_NOTES.md.

Run `python -m unittest discover -s tests -v` for the suite,
`python examples/run_examples.py` to reproduce the packaged fit results, and
`python examples/run_spectra_example.py` to reproduce the contributed
three-component spectral reconstruction. The latter creates numerical matrices,
normalization metadata, and a TRES/TRANES figure under
`example_results/spectral_reconstruction`.

Numerical test environment:

```json
{
  "python": "3.12.14",
  "numpy": "2.3.5",
  "scipy": "1.17.0",
  "pandas": "2.2.3",
  "matplotlib": "3.10.8",
  "threadpoolctl": "3.6.0"
}
```

## GitHub repository packaging checks

The first public release includes the application, numerical tests and repository support files. The 69-test suite passes in the repository copy. The project builds a Python wheel locally; this is a library/CLI packaging check, not a standalone GUI installer. CITATION.cff passes the official CFF 1.2.0 schema. Workflow and issue-template YAML, ignore rules, source-data line-ending rules and local documentation links were checked. GitHub-hosted Windows/Linux jobs have not run yet; their results become available after uploading the repository and running Actions.
