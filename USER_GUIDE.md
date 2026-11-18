# TCSPC Workbench 0.1.0 — experimental user's guide

This guide follows a measurement from exported photon counts to a saved lifetime fit, then, when appropriate, to global analysis and reconstructed TRES/TRANES. You can use the controls without writing or changing Python code.

The workbench runs on your computer. A Python process does the calculation, and a browser page displays the controls and plots. After the packages are installed, ordinary local analysis does not require an internet connection or cloud account.

## 1. Choose where to start

| Your situation | Start here |
|---|---|
| You already have a compatible Python environment | Keep its working environment and follow Section 2. |
| You have Anaconda Navigator but have never run the workbench | Follow Section 3 once, then Section 2 for later sessions. |
| The app is already open and you want to learn the fitting controls | Try the supplied single-decay example in Section 5. |
| You are ready to use your own decay and IRF | Prepare the files in Section 4, then follow Section 6. |
| You measured decays at several wavelengths | Follow global analysis in Section 9, then spectral reconstruction in Section 10. |

Keep the original acquisitions and a note of your time calibration, excitation repetition rate, detection settings, and any experimental corrections. The workbench needs the exported numbers and the correct units; a filename alone does not establish how a measurement was acquired.

## 2. Open an existing working installation

### Extract the complete version

1. Download the release ZIP to your computer.
2. In Windows File Explorer, right-click the ZIP and choose **Extract All**. Do not run the application while browsing inside the compressed ZIP.
3. Open the extracted folder. It may contain another folder named `TCSPC_Workbench`; open that too.
4. Stop when you can see `app.py`, `requirements.txt`, and a folder named `tcspc` together. Keep the whole folder, including its `examples` and documentation.

A `.py` file contains Python instructions. Double-clicking it may open a text editor or a terminal that closes immediately. Start the workbench with the command below instead. The `.html` guides can be opened directly in a browser without starting Python.

### Activate your working environment and launch

If another app server is still running, first press **Ctrl+C** in its terminal. A browser tab left open on an older server can keep showing the older version.

Open **Anaconda Prompt** from the Windows Start menu. In Navigator, the equivalent is **Environments → triangle beside your working environment → Open Terminal**.

Activate the environment that already runs the workbench. If you used the name from the installation guide, type:

```bat
conda activate tcspc_workbench
```

If your working environment has another name, replace `tcspc_workbench` with that name. A compatible existing environment can be reused. You can list installed environment names with:

```bat
conda env list
```

Copy the extracted program folder's location from File Explorer's address bar. Change into that folder, using your actual location in place of this example:

```bat
cd /d "C:\Users\YourName\Downloads\TCSPC_Workbench"
```

Check that you reached the correct folder:

```bat
dir app.py requirements.txt
```

Both files should be listed. Then start the app:

```bat
python -m streamlit run app.py
```

Your browser should open. Otherwise, open the **Local URL** printed in the terminal, usually `http://localhost:8501`. If it prints a different port, use that address. `localhost` means this computer; you do not need to connect to a remote workspace.

Check that the page says **Version 0.1.0**. Keep the terminal open while using the app. You may minimize it. Closing only the browser tab does not stop the Python process; press **Ctrl+C** in the terminal when finished. Download your results before stopping or refreshing the session.

For the next session, repeat activation, `cd /d`, and the launch command. Installation is normally a one-time step. After activating the environment, `START_WORKBENCH_WINDOWS.bat` is an optional shortcut from that same terminal; it does not select your environment for you.

## 3. First installation with Anaconda

Skip this section if your existing environment already has the required dependencies.

If Anaconda or Miniconda is not installed, use [Anaconda's official download page](https://www.anaconda.com/download), install the version matching your computer, and reopen the Start menu afterward. Navigator manages environments; commands are entered in **Anaconda Prompt** or an environment's terminal.

For 64-bit Windows on Intel or AMD, open Anaconda Prompt and paste this as one complete command:

```bat
conda create -n tcspc_workbench --override-channels -c defaults python=3.12 numpy scipy matplotlib pandas threadpoolctl "streamlit>=1.51,<2" "blas=*=mkl"
```

Review the package list and enter `y` when asked to proceed. Wait for the downloads and installation to finish. The quotes are part of the command. This installs the numerical and interface packages in a separate environment named `tcspc_workbench`.

Activate it and check the packages:

```bat
conda activate tcspc_workbench
python -c "import numpy, scipy.optimize, scipy.signal, matplotlib.pyplot, pandas, streamlit, threadpoolctl; print('Dependencies OK')"
```

Continue when `Dependencies OK` appears. Then follow Section 2 to change to the extracted folder and launch. Windows on ARM requires packages appropriate to that architecture; the MKL-specific command above is for Intel/AMD. Alternative Python, macOS, and Linux installation instructions are in [INSTALLATION.md](INSTALLATION.md).

If a package fails to import, resolve that environment error before fitting data. The error is occurring before the measurement calculation starts. Section 12 explains the common cases.

## 4. Prepare experimental files

### A decay or an IRF: two columns

Export time in the first column and raw photon counts in the second. Use a separate file for a shared measured IRF. A short format example is:

```csv
time_ns,counts
0.000,0
0.025,2
0.050,1
0.075,0
```

The small table illustrates the layout; it is not a complete fluorescence decay. Keep all intended time channels, including zero-count channels. For Poisson fitting, use raw nonnegative counts before normalization, smoothing, or background subtraction; background is represented in the fit model.

Choose the calibrated input unit in the app: `ns`, `ps`, `us`, `ms`, or `s`. A plain header such as `time_ns` is a useful label but does not automatically override the unit selector. The workbench reports fitted times and lifetimes in ns.

The decay and shared IRF need compatible time calibration and matching bin widths. Do not move the IRF by changing its time column just to make the display line up; the shift setting describes a timing offset within the model.

### Several wavelength decays: one common time column

For global analysis, put time in column 1 and a decay at each known wavelength in the subsequent columns:

```csv
time_ns,500,520,540
0.000,1,2,0
0.025,0,1,3
0.050,2,0,1
0.075,1,1,0
```

Each numeric header is its column's wavelength in nm. Rows are time bins. A file with 51 wavelength decays therefore has 52 columns including time. Do not add the shared IRF as another wavelength column in this plain matrix. Upload it separately.

If wavelength labels are absent from a plain matrix, provide one wavelength for each decay column in the app, in exactly that column order. Generic curve indices such as `0, 1, 2` are not wavelength calibration.

### Export from Origin or a similar worksheet

Use these settings for the simplest import:

| Export setting | Recommended value |
|---|---|
| Text encoding | UTF-8; plain ASCII characters are also valid UTF-8. |
| Column separator | Comma for CSV. |
| Decimal separator | Decimal point, such as `0.025`. |
| Numeric precision | Full precision. |
| Header | One row: time label followed by the real wavelengths, or time/count labels for a two-column file. |
| Extra row indices, comments, units and other worksheet-header rows | Exclude them from this simple table. |

In Origin, put the intended headers in the Long Name row and export that row. Turn off extra Short Name, Units, Comments, and row-index rows for this layout. Check the exported file in a plain-text editor: wavelength labels must be present in the file, not only visible in the worksheet.

Consistent tabs, semicolons, or spaces are also supported by the text reader, but comma CSV is the recommended starting point. Renaming a binary measurement file to `.csv` does not convert it. Export numeric text from the acquisition software.

Small rounding differences in printed TCSPC times are handled automatically. Counts remain in their original channels. You normally do not need to rebuild an experimentally uniform time axis. Duplicate, missing, or genuinely irregular channels need investigation. Detailed bundle and source-curve formats are in [DATA_FORMATS.md](DATA_FORMATS.md).

## 5. Learn with the supplied single-decay example

This is a useful first session because both decay and IRF are supplied.

1. Choose **Analysis → Single decay** in the left sidebar.
2. Choose **Source → Bundled decay + IRF**.
3. Under **Method**, choose **IRF reconvolution**. Use **Exponential components = 1** and **Poisson likelihood (raw counts)**.
4. Open **Background, shift, and excitation**. Keep **Fit a background for each decay** selected. Set **Subtract IRF background (counts/bin)** to `1.4`.
5. Select **Limit IRF pulse range**. Set its first and last included pulse times to `3.8` and `8` ns. Select **Fit one shared IRF shift**; the default shift bounds can remain for this example. Leave repetition-period correction off.
6. Inspect **Kinetics preview and fit interval**. Start with the whole imported time record. Try moving the range handles and using the numerical time boxes, as explained in Section 7.
7. For the reference worked result, return the fitting interval to the whole imported record, then press **Run fit**.
8. Inspect **Most recent fit**, its fitted curve, residuals, autocorrelation, and diagnostics. The worked lifetime is near **2.549 ns** with these particular settings.
9. Select **Download fit results (JSON + CSV ZIP)** and save the ZIP in a folder for this analysis.

The supplied text example assumes ns calibration. The baseline, pulse range, and approximate lifetime above apply to this example only. They are not instrument constants or known ground truth for your own measurements. Selecting a different fit window can change the result; compare the settings as well as the number.

## 6. Fit your own normal TCSPC decay

1. Select **Analysis → Single decay**, then **Source → Upload ASCII data**.
2. Choose the correct **Time unit in uploaded files**.
3. Upload your fluorescence file under **Decay or TRES ASCII file**. For a single decay it contains time and counts only.
4. Upload the corresponding measured **IRF ASCII file (optional for tail fit)** when using reconvolution.
5. Open **Inspect input** if you want to check the first imported rows. Check the loaded time range and bin-width readout against your acquisition.
6. Choose the fitting method using the table below. Start with one exponential component if that is a reasonable model for the experiment.
7. Set the fitting interval in the kinetics preview. Enter the background and pulse information appropriate to your measurement; the example's values are not defaults for an unknown IRF.
8. Press **Run fit**, inspect the result and residuals, then download the result ZIP.

| Method | When to use it | Range selection |
|---|---|---|
| **IRF reconvolution** | You have a measured IRF and want the model to account for the excitation/detector response. | Choose the relevant observed decay region. The model can include the response-affected rise; retain the full imported source record for convolution. |
| **Exponential tail fit** | You deliberately fit a later region where the IRF-affected rise is negligible. An IRF upload is unnecessary. | Start after that rise and fit a useful falling region with enough signal to distinguish decay from background. The automatic initial start is only a starting suggestion. |

For reconvolution, do not automatically discard the rising portion merely because a tail fit would discard it. For a tail fit, do not start in the response-affected rise. If changing the range changes a lifetime substantially, investigate the model, response, background and signal quality rather than treating one selected window as proof.

Use additional components when supported by the measurement and residual behavior. A more flexible curve can look better while its individual lifetimes and amplitudes remain poorly determined. Comparing different component models is most useful with the same observations, fit range and fit statistic.

## 7. Select and inspect the fitting interval

The same **Kinetics preview and fit interval** panel is available for normal single-decay TCSPC and global fitting.

1. Set **Kinetics preview view → Whole record** to locate the pulse, decay and late background.
2. Move the two ends of **Fit interval handles (time bins)** immediately above the preview curve. The numerical start/end boxes remain in the left settings column. The left handle sets the first included channel; the right sets the last. The handle numbers are channel indices, while the linked time boxes and readout show ns.
3. For a precise limit, enter **First fitted time (ns)** or **Last fitted time (ns)**. Values snap to the nearest available channel. The handles, numerical fields, shaded region and vertical guides update together.
4. Select **Kinetics preview view → Selected fit window** to enlarge that chosen region. Switch back to **Whole record** to see its position in the acquisition.
5. Check the displayed interval and included-channel count, then press **Run fit**. At least four channels are required; useful experimental windows usually need many more.

The shaded region marks the channels used by the next fit. A zoom changes the displayed x-axis; it does not independently choose a different fitting interval. The existing **Most recent fit** is updated only after another fit is run. Its caption records the interval actually used, so you can compare a new preview selection with the previous result.

For a global dataset, **Kinetics preview trace** chooses the sum or one wavelength for display. The selected fit interval still applies to every selected decay, not just the preview trace. This trace selector is unnecessary for a single decay.

The IRF overlay is scaled to the decay height for display only. Determine its background from the original IRF counts in an independently selected baseline region, not from the height of that scaled overlay.

Choosing a window does not crop, rebin, normalize, or move the imported counts. Reconvolution uses the full imported record before applying the selected fitting interval to the comparison with the data.

Compare the [whole-record interval example](docs/interface_single_range.png) with the [selected-window zoom](docs/interface_single_zoom.png).

## 8. Inspect the fitting result

| Result | What it tells you |
|---|---|
| **Convergence** | Whether the numerical optimizer met its stopping conditions. A converged fit can still use an unsuitable model. |
| **Lifetime (ns)** | The fitted exponential time constants under the selected model and settings. |
| **Approximate standard error** | Local parameter uncertainty conditional on the chosen model and measured IRF; it does not include calibration or model errors. |
| Data/fit overlay | Whether the predicted counts describe the observed decay over the included range. A log display makes late features easier to see. |
| Residual plot | Differences under the chosen statistic. Look for a pattern at the rise, a long trend, oscillation or a late offset rather than only a smooth fit line. |
| Residual autocorrelation | Whether neighboring residuals retain a pattern. Strong structure can indicate mismatch; random-looking residuals alone do not identify a physical mechanism. |
| Reduced Pearson χ² and reduced Poisson deviance | Different summaries of agreement. Their meanings differ and a value near 1 alone is not proof of the model. |
| Amplitude spectra in a global result | The fitted pre-exponential factors at each wavelength. These are inputs to reconstruction, not the reconstructed TRES themselves. |

Keep a note of the model size, method, interval, time calibration, IRF preparation and background assumptions with the result. If the optimizer reports nonconvergence or a parameter is at a bound, investigate before relying on that parameter. [SCIENTIFIC_NOTES.md](SCIENTIFIC_NOTES.md) explains the mathematical conventions and limits.

## 9. Fit several wavelength decays globally

Global analysis uses shared lifetime components for all selected wavelength decays while fitting separate amplitudes and backgrounds at each wavelength.

1. Choose **Analysis → Global TRES** and **Source → Upload ASCII data**.
2. Upload a time-by-wavelength matrix from Section 4 and the matching shared IRF for reconvolution.
3. Choose the time unit. With a numeric wavelength header, leave the wavelength text field blank. Otherwise enter the known wavelengths in decay-column order.
4. Choose the component count and fitting settings. Inspect the heatmap, then a wavelength's kinetics using **Kinetics preview trace**.
5. Use the shared fitting-window handles, numerical fields and zoom controls from Section 7.
6. Press **Run fit**. Inspect several wavelengths in **Inspect fitted trace**, their residuals and amplitude spectra.
7. Save **Download fit results (JSON + CSV ZIP)**. Use **Continue to TRES / TRANES** for the next calculation.

To try a measured example, upload `examples/measured_TRES.csv` and `examples/measured_IRF.csv`, select ns, and use the header wavelengths. The worked settings in [INSTALLATION.md](INSTALLATION.md) reproduce the included two-component example. That measured dataset has no supplied matching steady-state spectrum; it cannot supply the missing fluorescence intensities by itself.

For multi-curve ASCII bundles, select the actual IRF and decay curve indices in the import controls. A generic curve index does not tell the workbench its wavelength or polarization role. [DATA_FORMATS.md](DATA_FORMATS.md) gives the supplied examples' assignments.

## 10. Continue to TRES and TRANES

The completed global fit already contains wavelengths, amplitudes and paired lifetimes. Click **Continue to TRES / TRANES** and keep **Fitted-decay source → Latest global fit**. You do not need to create and upload another amplitude table or retype its lifetimes.

The global result remains available if you subsequently fit a separate single decay in this running session. A restarted or reloaded session is not a saved analysis archive; download the outputs.

### Choose the matching steady-state input

The file-upload paths use a two-column file with wavelength in nm and fluorescence intensity:

```csv
wavelength_nm,intensity
500,0.25
520,0.80
540,0.60
```

These are layout examples, not intensities for a packaged experimental dataset.

| Steady-state input | Use when | Calculation |
|---|---|---|
| **Normalized values at fitted wavelengths** | You already have normalized intensities at exactly the wavelengths of the global fit. | Rows are matched by wavelength and used directly, without interpolation or another rescaling. Include one unique row per fitted wavelength; row order may differ. |
| **Full steady-state spectrum** | You have an emission spectrum with its own wavelength sampling. | Normalize the full spectrum, then match or linearly interpolate its values at the global-fit wavelengths. It must cover all those wavelengths. |
| **Enter normalized intensities manually** | You have normalized intensities to type or paste. | Fills and locks the wavelength column; uses the completed intensity column directly after Apply intensities. |

In the first path, intensities must be between 0 and 1, with some positive signal. Their maximum can be below 1 if your fitted wavelengths cover only part of an already normalized spectrum. Missing, duplicate or extra wavelengths do not belong to that exact-matching path; use the full-spectrum path for a denser spectrum.

In the full-spectrum path, select **Normalization method**:

- **Divide by maximum (I / max)** divides all intensities by the full uploaded spectrum's maximum. It requires nonnegative signal and a positive maximum.
- **Scale minimum to 0 and maximum to 1** uses `(I − minimum) / (maximum − minimum)`. This subtracts the minimum before scaling, so it changes relative spectral shape when that minimum is not zero. It can accommodate a finite negative raw offset but cannot normalize a constant spectrum this way.

Use the spectrum for the same sample and compatible measurement conditions. Record blank, spectral-response and other experimental corrections. Min–max scaling is a chosen normalization, not a measured blank correction. Values outside measured wavelength coverage are not invented by extrapolation.

Inspect **Steady-state intensities matched to the global fit**. Open **Inspect full steady-state spectrum and normalization** to compare the original values, normalized spectrum and matched intensities.

### Enter normalized intensities directly

Choose **Steady-state input → Enter normalized intensities manually**. The table already contains the fitted wavelengths; that column is locked. Type or paste normalized fluorescence intensities into the intensity column and click **Apply intensities**. Fill every row with a finite value between 0 and 1; zeros are allowed, but the whole spectrum cannot be zero. Blank cells stay blank and must be completed.

Values are used directly, without interpolation or renormalization. Their maximum can be below 1 if the fitted wavelengths exclude the peak of a previously normalized spectrum. After editing a value, click **Apply intensities** again to update the spectra. The spectral ZIP records the entered values in `sampled_steady_state.csv` and identifies manual input in its metadata. No steady-state file upload is needed for this path.

See the [manual-entry table example](docs/interface_manual_table.png).

### Choose delays and plot connections

Enter nonnegative delays after excitation in **Reconstruction delays (ns)**. Move **Preview delay after excitation (ns)** to inspect one delay. **Maximum preview delay (ns)** sets that slider's extent; it does not change the fitted lifetimes. The initial delays for a linked fit scale with its longest lifetime. Times outside the measured kinetics depend on model extrapolation.

Choose **Spectral curve style**:

| Style | Display |
|---|---|
| **Smooth curve (PCHIP)** | A smooth connection that preserves the shape between nonnegative calculated samples. |
| **Straight lines** | Direct line segments between calculated samples. |
| **Cubic B-spline** | A smooth interpolating spline. It may overshoot between samples, including below zero; with fewer than four wavelengths a lower degree is used. |

The markers and exported matrices retain the calculated values at the fitted wavelengths. Changing the connection does not refit the kinetics or create new experimental spectral resolution. A spline dip below zero can be a display artifact rather than negative fluorescence.

TRES is the intrinsic model's relative emission spectrum at the selected delay, scaled using the steady-state intensities and each wavelength's fitted decay area. TRANES divides each TRES by its own wavelength area. Its sampled trapezoidal area is one over the fitted wavelength interval; it is not peak normalization. A smooth drawn connection can have a different area from that sampled-grid convention.

Tail-fit amplitudes refer to the first fitted tail channel. If reconstructing from a tail fit, the app asks for that channel's independently known elapsed delay after excitation before back-extrapolation. An absolute detector-time label alone does not determine that delay.

Save **Download TRES + TRANES matrices** after inspecting the spectra. The optional **Example or external fitted table** route is available for a separate worked example or results fitted elsewhere. It is unnecessary when continuing your own current global fit.

## 11. Save and organize the analysis

The workbench does not overwrite the input measurement files. Save the results deliberately before closing the app, refreshing the page or stopping Python.

| Download | Main contents |
|---|---|
| **Download fit results (JSON + CSV ZIP)** | `result.json` with fitted parameters/settings/diagnostics; `curves.csv` with observations, predictions, residuals and inclusion flags; `amplitudes.csv`; import and wavelength mappings when available. |
| **Download TRES + TRANES matrices** | `tres.csv`, `tranes.csv`, matched and normalized steady-state tables, fitted amplitudes, normalization/matching metadata and reconstruction-source information. |

For a global spectral analysis, save both ZIP files. A practical folder can contain the original decay matrix, IRF, steady-state spectrum, fit-result ZIP, reconstruction ZIP and an experimental-settings note. Use distinct filenames when comparing component counts or fitting windows so that one downloaded result does not hide another.

Raw fitting input uses **time rows and wavelength columns**. Reconstructed spectral output uses **wavelength rows and delay columns**. Import the latter into Origin as wavelength x versus the chosen delay columns. `sampled_steady_state.csv` holds intensities at fitted wavelengths; `normalized_steady_state.csv` holds the normalized complete source spectrum. Both are two-column tables that can be read again as spectrum inputs.

## 12. Troubleshooting

| What you see | What to do |
|---|---|
| `conda` is not recognized | Open Anaconda Prompt or the working environment's terminal in Navigator. |
| Environment name not found | Run `conda env list` and activate the environment you actually installed. Do not guess that an update requires a new one. |
| `app.py` cannot be found | Open the inner extracted program folder; `dir app.py requirements.txt` should show both files. |
| `No module named streamlit`, NumPy or SciPy | Check that the terminal uses the intended environment. Run the dependency import check in Section 3 there. |
| NumPy C-extension/DLL error with Application Control blocking a library | This is an installation/library-policy problem before fitting begins. Record the exact blocked filename and error. Use a trusted compatible package build; the official Conda/MKL route is an alternative build, not a guarantee that a local policy accepts it. Disabling security controls is not part of this guide. |
| A browser tab still shows an older version | Stop the old Python server with Ctrl+C and launch from the new folder. Use the new process's Local URL and check the displayed version. |
| The command is running but no browser opened | Open the Local URL printed in that terminal manually. Keep the process running. |
| Wavelength labels are requested | Supply a numeric header or one wavelength per decay column. Check the exported text, not only the worksheet's column properties. |
| A CSV import fails | Check encoding, separator, decimal point, extra headers, missing cells and whether the first column really contains time or wavelength for that file type. See Section 4. |
| The time axis is rejected | Check selected units, duplicate/missing channels and acquisition calibration. Tiny printed rounding is handled automatically; do not rebin a calibrated uniform histogram just to alter displayed decimals. |
| Run fit is disabled | Expand the selected interval to at least four time channels. For experimental analysis, choose enough of the decay to support the model. |
| The preview changed but the reported lifetime did not | Press Run fit again. Most recent fit remains the previous completed result until another fit runs. |
| IRF background is zero | Supply the independently measured IRF baseline if required. It is separate from the decay background. |
| Reconvolution asks for an IRF | Upload the measured response, or deliberately select Exponential tail fit and a suitable later region. |
| Normalized steady-state wavelengths do not match | Include exactly one row per fitted wavelength, or choose Full steady-state spectrum for a denser grid. |
| Steady-state coverage is too short | Supply a spectrum spanning every fitted wavelength. The workbench does not extrapolate beyond the measured spectrum. |
| A B-spline goes below zero between points | Select PCHIP or straight lines to inspect the same nonnegative calculated samples without that spline overshoot. |
| A fit does not converge or changes strongly with starting values | Inspect the range, backgrounds, response, component count, bounds and residuals. More optimizer attempts can help explore solutions; they do not establish the physical model. |

## 13. Short glossary

| Term | Meaning here |
|---|---|
| TCSPC time bin/channel | One interval of measured photon delay after an excitation reference. Counts belong to that channel. |
| IRF | Instrument response function: measured timing response to a short excitation event. |
| Reconvolution | Predicting the measured shape by combining an intrinsic exponential model with the measured IRF. |
| Lifetime | A fitted exponential time constant under the selected model. |
| Pre-exponential amplitude | Weight multiplying a lifetime component at a wavelength. |
| Global fit | Several decays fitted together with shared lifetime components. |
| Fit window | The selected time channels entering the comparison of model and data. |
| Preview zoom | The portion of the x-axis shown on screen; it is separate from selecting channels for fitting. |
| Residual | The remaining difference between model and observation under the chosen statistic. |
| Steady-state spectrum | Fluorescence intensity as a function of wavelength under steady measurement conditions. |
| TRES | Reconstructed time-resolved emission spectrum at a chosen delay. |
| TRANES | A reconstructed spectrum normalized to unit area over wavelength. |

## 14. Optional fitting-control reference

Use this table when you need to understand a setting beyond the basic workflow.

| Control | Meaning and practical use |
|---|---|
| **Exponential components** | Number of exponential decay components. Start with a justified small model; extra components require evidence and adequate time resolution/signal. |
| **Fit statistic** | Poisson likelihood suits raw photon counts, including zero-count bins. The other choices fit with observed-count weights or without weights. Keep the statistic consistent when comparing fits. |
| **Initial lifetimes (ns; optional)** | Starting guesses, such as `0.5, 3.0`; the optimizer may change them. Leave blank to use the workbench's starts. Supply one value per component when entering them. |
| **Fixed lifetimes (ns; optional)** | Values held constant. These override initial guesses. Use them only when fixing the components is part of the intended analysis. |
| **Fit a background for each decay** | Fits a constant background separately for each selected trace. Otherwise the entered decay-background value is fixed. |
| **Decay background (counts/bin)** | Fluorescence histogram baseline in counts per channel. It belongs to the decay model, not to the IRF correction. |
| **Subtract IRF background (counts/bin)** | Independently measured IRF baseline to subtract before normalizing the IRF pulse. Estimate it from an appropriate baseline region of that measurement. |
| **Limit IRF pulse range** | Selects the measured pulse and relevant response tail while excluding unused baseline. Do not cut away a real response tail merely to make the pulse look narrow. |
| **IRF shift / Fit one shared IRF shift** | Fixed or fitted timing offset between the measured response and decay. Positive values delay the IRF. Bounds should be consistent with plausible acquisition timing. |
| **Use repetition-period correction** | An advanced model for repeated excitation on a complete period grid. Enter the measured laser repetition period only when the imported grid represents one complete, nonduplicated period. Leave off when that condition is not met. |
| **Starting-point attempts** | Tries several optimizer starts. More attempts can help explore alternative solutions but increase run time. |
| **Evaluation limit per attempt** | Maximum optimizer work per start. Increasing it does not correct an unsuitable physical model. |

A zero IRF-background warning means no independently supplied baseline has been subtracted from the response. It does not mean the fluorescence decay background is necessarily zero. Determine the IRF baseline from its measurement instead of copying the decay background into that field.


For file layouts and curve assignments, read [DATA_FORMATS.md](DATA_FORMATS.md). For equations, assumptions and interpretation limits, read [SCIENTIFIC_NOTES.md](SCIENTIFIC_NOTES.md). Installation alternatives are in [INSTALLATION.md](INSTALLATION.md), and the checks performed for this release are in [VALIDATION.md](VALIDATION.md).
