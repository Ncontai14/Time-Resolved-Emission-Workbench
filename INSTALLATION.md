# TCSPC Workbench 0.1.0

TCSPC Workbench fits fluorescence photon-counting decays with a measured instrument response function (IRF). It also fits wavelength-resolved and repeated-decay series globally, with shared lifetimes and separate amplitudes and backgrounds for each trace. With an accompanying steady-state fluorescence spectrum, it reconstructs steady-state-scaled time-resolved emission spectra (TRES) and their area-normalized forms (TRANES).

The application runs on your computer and displays its controls in a web browser. The browser is the interface; Python performs the calculations. The package includes a command-line interface and a Python API for more advanced use.

For a step-by-step guide written for experimental users, open [USER_GUIDE.md](USER_GUIDE.md) or double-click `USER_GUIDE.html`. It covers the first installation, daily startup, fitting controls, fit-window selection, global analysis, TRES/TRANES and saving results without editing Python code. This README retains installation alternatives and reference workflows. Each guide has an offline HTML edition that opens directly in your browser.

## Using an existing compatible environment

If your Conda or Python environment already has the required dependencies, you can use it. Stop any running app with **Ctrl+C**, activate that environment, and launch `app.py` from this extracted folder. The environment-creation instructions below are for a first-time installation.

Check that the page displays **Version 0.1.0**. Opening a browser tab left connected to an older running server can show the old version.

## Before starting

You need the **whole extracted package**, a compatible Python environment, and an internet connection for the first dependency installation. Ordinary local fitting does not require a cloud account.

1. Download the ZIP and choose **Extract All** on Windows, or extract it with your usual archive tool.
2. Open the extracted folder. If it contains another folder, open that too until you see `app.py`, `requirements.txt`, and the `tcspc` folder together.
3. Keep those files together. The `.py` extension means Python source code. Opening `app.py` in a text editor shows the program; the launch command below starts its interface.
4. Choose **one installation route** below. If you already use Anaconda Navigator, use the recommended Anaconda route.

A terminal is a window in which you type commands and press Enter. On Windows, **Anaconda Prompt** is available from the Start menu after installing Anaconda. Commands shown in code blocks are copied into that window, one command at a time. Do not copy the surrounding explanation or any prompt such as `(base)`.

## Recommended Windows installation with Anaconda

This route creates a separate environment named `tcspc_workbench`. An environment is a folder containing a Python installation and the packages for this application. It avoids changing packages used by your other projects.

The command below uses packages from Anaconda's `defaults` channel and selects the MKL numerical-library build. This is intended for 64-bit Windows on Intel or AMD processors. Windows on ARM needs compatible packages for that architecture; the MKL selector is not appropriate there.

### 1. Open Anaconda Prompt

If Anaconda is not installed, obtain Anaconda Distribution or Miniconda from [the official Anaconda download page](https://www.anaconda.com/download). Choose the installer matching your Windows processor and follow its installation prompts. Reopen the Start menu after installation. Navigator is a convenient manager; the commands below run in its Python environment, not in a Navigator text editor.

Use **Start → Anaconda Prompt**. You can also use Anaconda Navigator: open **Environments**, click the triangle beside an environment, and choose **Open Terminal**.

### 2. Create the application environment

Paste this as **one complete command**:

```bat
conda create -n tcspc_workbench --override-channels -c defaults python=3.12 numpy scipy matplotlib pandas threadpoolctl "streamlit>=1.51,<2" "blas=*=mkl"
```

When asked whether to proceed, review the package list and enter `y`. Downloading and resolving packages may take several minutes. Wait until the command finishes and a prompt returns.

The quotes around package expressions are important in Windows terminals. This command installs all the application dependencies, including Streamlit. It replaces the pip-based dependency installation used in the alternative route below.

An equivalent environment specification is supplied as `environment-windows.yml`:

```bat
conda env create -f environment-windows.yml
```

Use the command-line package list **or** the environment file to create the environment, not both. The environment-file command must run from the extracted program folder or use that file's full path.

### 3. Activate the environment

```bat
conda activate tcspc_workbench
```

The prompt should now begin with `(tcspc_workbench)`. Activation tells the terminal which Python installation to use.

### 4. Check the installation

```bat
python -c "import numpy, scipy.optimize, scipy.signal, matplotlib.pyplot, pandas, streamlit, threadpoolctl; print('Dependencies OK')"
```

Continue when this prints `Dependencies OK`. An error at this stage concerns the Python environment or its libraries, before any measurement is fitted. See Troubleshooting below if it fails.

### 5. Go to the program folder

In File Explorer, open the folder containing `app.py` and copy its location from the address bar. In Anaconda Prompt, type `cd /d` followed by that path in quotation marks. For example:

```bat
cd /d "C:\Users\YourName\Downloads\TCSPC_Workbench"
```

**Replace the example path with your actual extracted folder.** The `/d` option also changes drive if the package is on another drive. Check that you are in the correct folder:

```bat
dir app.py requirements.txt
```

Both files should be listed.

### 6. Start the application

```bat
python -m streamlit run app.py
```

Your browser should open automatically. Otherwise, open the **Local URL** printed in the terminal, usually:

```text
http://localhost:8501
```

Use the printed port if it is different. Keep the terminal open while fitting data. Closing the browser tab does not stop the program; press **Ctrl+C** in the terminal to stop it.

After activating the environment and entering the program folder, you may alternatively run `START_WORKBENCH_WINDOWS.bat`. This launcher starts an already installed environment; it does not install dependencies or select a Conda environment for you.

## Opening the application next time

Open Anaconda Prompt, then run:

```bat
conda activate tcspc_workbench
cd /d "C:\Users\YourName\Downloads\TCSPC_Workbench"
python -m streamlit run app.py
```

Use your actual folder path. The environment creation and package installation are normally one-time steps.

If you already have a compatible environment, activate it and launch from the extracted program folder. Use your actual environment name:

```bat
conda activate YOUR_EXISTING_ENVIRONMENT
```

Replace `YOUR_EXISTING_ENVIRONMENT` with the installed environment name. A separate dependency update is needed only if the package check reports an unmet requirement. Then change to the new application folder and run `python -m streamlit run app.py`. Stop any older running copy first so that you do not keep viewing an old browser session.

## Alternative: standard Python without Anaconda

Use Python 3.12 for the most directly exercised installation. These commands create a local virtual environment called `.venv` inside the program folder. Use this route when you are not using the Conda setup above.

If Python is not installed, use [the official Python download page](https://www.python.org/downloads/) and choose a compatible installer. For a conventional Windows Python installation, enable **Add Python to PATH** when the installer offers it. Open a new terminal afterward so that it picks up the updated command paths.

On Windows, open a terminal in the folder containing `app.py`. In File Explorer, click that folder's address bar, type `powershell`, and press Enter. Then run:

```bat
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe -m streamlit run app.py
```

Using the environment's Python directly avoids PowerShell activation-policy problems. If the `py` launcher is unavailable but Python is installed, use `python -m venv .venv` for the first command.

On macOS or Linux:

```sh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

A general Conda specification is also supplied as `environment.yml`. Once a compatible environment is activated on macOS or Linux, `bash start_workbench.sh` launches the app from its folder.

## First fit: a single decay

1. In the sidebar, choose **Single decay** and **Bundled decay + IRF**.
2. Select IRF reconvolution, one exponential component, and Poisson likelihood.
3. In the background and excitation settings, enter an IRF background of **1.4 counts/bin**, enable the IRF pulse-range limit, and use **3.8–8 ns**.
4. Enable fitting of one shared IRF shift. Keep the decay background fitted.
5. In **Kinetics preview and fit interval**, inspect the whole record, move **Fit interval handles (time bins)** immediately above the curve or enter the first/last fitted times, and use **Kinetics preview view → Selected fit window** to enlarge the selection. The linked fields and shaded interval apply to a normal single decay as well as a global fit. Return to the whole imported fitting interval for this reference example.
6. Click **Run fit**. This worked example gives a lifetime near **2.549 ns** with these particular settings.
7. Inspect the data/fit overlay, residuals, autocorrelation, warnings, and fit statistics. Use the download button to save the result ZIP.

The IRF baseline and pulse limits above are illustrative choices for this sample. They are not instrument calibration constants for other measurements. The contributed measured examples have no supplied ground-truth lifetimes.

The time boxes, range handles and preview remain linked. Zooming changes the displayed x-axis; press **Run fit** after choosing a different fitting interval to update the result. See the [single-decay interval example](docs/interface_single_range.png) and the full operating procedure in [USER_GUIDE.md](USER_GUIDE.md).

For your own plain decay, choose **Upload ASCII data**, upload a time/count table, and upload the corresponding IRF. Choose the correct input time unit. For a tail fit, use a post-rise range where the IRF contribution is negligible; an IRF is then unnecessary.

## First global TRES fit

Select **Global TRES** and **Upload ASCII data**, then use:

- Decay file: `examples/measured_TRES.csv`.
- IRF file: `examples/measured_IRF.csv`.
- Input time unit: **ns**.
- Wavelength field: leave blank; the CSV header contains the 51 wavelengths from 490 to 590 nm in 2 nm steps.

For the packaged worked example, use two exponential components, initial lifetimes `2.5, 4.0` ns, Poisson likelihood, and fitted backgrounds. Set the IRF background to **0.38 counts/bin**, pulse range to **3.2–6.5 ns**, and fit a shared shift within **−0.5 to 0.5 ns**. Leave repetition-period correction disabled for this cropped example grid.

Before fitting, inspect an individual wavelength's kinetics in the global preview. Move the two ends of **Fit interval handles (time bins)** or enter **First fitted time (ns)** and **Last fitted time (ns)** numerically; the controls stay linked. The shaded region and vertical guides show which channels enter the objective. Use **Kinetics preview trace** to inspect the sum or an individual wavelength, and **Kinetics preview view** to choose **Whole record** or **Selected fit window**. A display zoom does not independently change the fitting interval.

See the [fit-range interface example](docs/interface_range.png) for the linked controls and shaded preview.

The shared lifetimes are approximately **2.334 and 3.882 ns** under that model. Residual statistics show some model mismatch. These are worked fits, not reference values or proof of a particular physical mechanism. The amplitude spectra have no automatic detector-response or acquisition-duration correction.

Raw multi-curve examples are under `examples/curve_exports/`. When uploading `TRES_IRF.dat` from that directory, assign curve `0` as the IRF and use the inclusive decay range `1-51`. Generic curve labels do not encode wavelengths; supply known labels if you want physical wavelength axes.

## Reconstructing TRES and TRANES

There are four different quantities to keep separate:

| Quantity | Meaning |
|---|---|
| Measured time-by-wavelength counts | The imported TCSPC histograms, including the instrument response and background. |
| Decay-associated amplitudes | The fitted pre-exponential factors plotted against wavelength. These alone are not reconstructed TRES. |
| Steady-state-scaled TRES | Intrinsic model intensity versus wavelength at selected delays, scaled to the measured steady-state fluorescence spectrum. |
| TRANES | Each reconstructed spectrum divided by its own area over wavelength, so its area is one. |

For wavelength `λ`, the workbench uses the fitted lifetimes and amplitudes to calculate

```text
I(λ,t) = SS(λ) × [Σ A_k(λ) exp(-t/τ_k)] / [Σ A_k(λ) τ_k]
TRANES(λ,t) = I(λ,t) / integral[I(λ,t) dλ]
```

The denominator `Σ A_k τ_k` is essential: it makes the time integral of the intrinsic decay equal to the supplied steady-state intensity. Multiplying unscaled pre-exponential factors by steady-state intensity alone would retain arbitrary per-wavelength decay scaling.

### Continue from a global fit

The global fit already provides the wavelengths, pre-exponential amplitudes, component lifetimes, and their pairing. A later single-decay fit does not replace the saved global result. Continue from that global result; there is no need to export and re-upload an amplitude table or enter its lifetimes again.

1. Complete the wavelength-resolved **Global TRES** fit and inspect its residuals. Physical wavelength labels are required for reconstruction.
2. Click **Continue to TRES / TRANES** beneath the result. On the spectral page, **Latest global fit** is the initial **Fitted-decay source**. It supplies the fitted inputs automatically from this application session.
3. Under **Steady-state input**, choose a file-upload path or enter normalized intensities directly in the manual table.
4. Inspect **Steady-state intensities matched to the global fit**. Open **Inspect full steady-state spectrum and normalization** to compare the supplied intensities, normalized spectrum, and values assigned to the global-fit wavelengths.
5. Enter up to 64 unique nonnegative **Reconstruction delays (ns)**, separated by commas or spaces. Move **Preview delay after excitation (ns)** for a single-time preview; its slider extent is controlled by **Maximum preview delay (ns)**.
6. Choose how plotted points are connected, inspect both spectra, then use **Download TRES + TRANES matrices**. The numerical matrices preserve the calculated values at the fitted wavelengths. Initial delay values for a linked fit scale with its longest lifetime; edit them for your scientific question.

| Steady-state input | What to upload | What the workbench does |
|---|---|---|
| **Normalized values at fitted wavelengths** | Wavelength in nm and normalized fluorescence intensity, one row per global-fit wavelength. | Matches the rows to the global result and uses those values directly, without further normalization or interpolation. Row order may differ. |
| **Full steady-state spectrum** | Wavelength in nm and measured fluorescence intensity; it can contain many more wavelength points than the lifetime measurements. | Normalizes the full spectrum using your chosen method, then samples it at the global-fit wavelengths by linear interpolation within its measured coverage. |
| **Enter normalized intensities manually** | No file; type or paste values into the table. | Keeps the fitted wavelengths fixed and uses entered intensities directly after Apply intensities. |

For the first path, upload the file under **Normalized steady-state intensities**. Intensities must be finite values from 0 to 1, with some positive signal. Every fitted wavelength needs one matching row, and unrelated extra rows are rejected. The largest supplied value need not be 1 when the global-fit wavelengths cover only part of a previously normalized spectrum.

For the full-spectrum path, upload **Steady-state fluorescence spectrum**. Under **Normalization method**, choose **Divide by maximum (I / max)** or **Scale minimum to 0 and maximum to 1**. These calculate `I / max(I)` and `(I - min(I)) / (max(I) - min(I))`, respectively. Normalization happens before sampling, so a maximum between two lifetime-measurement wavelengths is retained in the scaling. Min–max normalization subtracts the supplied spectrum's minimum; that changes spectral shape unless the minimum was already zero. It is a user-selected scaling convention, not an inferred physical blank correction. A constant spectrum cannot use min–max normalization. Values outside the measured wavelength coverage are not extrapolated.

The 51-wavelength measured example has no corresponding steady-state spectrum supplied. Its global fit can be inspected, but it needs the matching experimental fluorescence spectrum to reconstruct TRES. The three-component example below belongs to a different dataset and must not be used as its steady-state spectrum.

### Enter normalized intensities directly

Choose **Steady-state input → Enter normalized intensities manually**. The table already contains the fitted wavelengths; that column is locked. Type or paste normalized fluorescence intensities into the intensity column and click **Apply intensities**. Fill every row with a finite value between 0 and 1; zeros are allowed, but the whole spectrum cannot be zero. Blank cells stay blank and must be completed.

Values are used directly, without interpolation or renormalization. Their maximum can be below 1 if the fitted wavelengths exclude the peak of a previously normalized spectrum. After editing a value, click **Apply intensities** again to update the spectra. The spectral ZIP records the entered values in `sampled_steady_state.csv` and identifies manual input in its metadata. No steady-state file upload is needed for this path.

### Choosing the spectral curve style

The **Spectral curve style** selector offers **Smooth curve (PCHIP)**, **Straight lines**, and **Cubic B-spline** connections. PCHIP is the initial smooth connection and preserves the shape between nonnegative samples. A cubic B-spline can overshoot and draw values below zero between samples; those values are displayed rather than clipped. With fewer than four wavelengths, the B-spline uses a lower degree.

All styles pass through the same calculated wavelength samples. Style selection changes the connection drawn between them; it does not refit the kinetics, supply additional measured wavelengths, or change the exported TRES/TRANES matrices. TRANES normalization uses the actual sampled wavelength grid. The area under a smooth drawn curve may differ from that discrete trapezoidal area. The selected style is recorded with the reconstruction export.

### Optional external fitted results and the worked example

Under **Fitted-decay source**, **Example or external fitted table** remains available for results calculated outside this running workbench. It is an optional route; the normal global-fit continuation uses its saved outputs directly. The supplied fitted table is an example of this external layout.

The **three-component worked example** uses `examples/fitted_three_components.csv`. It contains 13 wavelengths, from 498 to 653 nm, pre-exponential factors, and the corresponding relative steady-state intensities. For this example only, the order is:

| Amplitude column | Lifetime (ns) |
|---|---:|
| A1 | 114.15 |
| A2 | 17.415 |
| A3 | 4.613 |

Choose **Example or external fitted table**, then choose the worked example within that route. Its embedded steady-state values allow the reconstruction to run at delays such as `0.025, 1, 10, 50, 100` ns. These are contributed fitted data for checking the calculation, not known ground-truth lifetimes. Do not apply those lifetimes to another file or reorder only the lifetimes while keeping amplitude columns fixed.

The [reconstruction interface example](docs/interface_spectra.png) shows the spectral workflow. To reproduce the supplied example without opening the app, run this command from the program folder in your activated environment:

```sh
python examples/run_spectra_example.py
```

The [worked example script](examples/run_spectra_example.py) saves both spectral matrices, metadata, and a TRES/TRANES plot under `example_results/spectral_reconstruction`. It uses the contributed fitted table directly; it does not refit the separate measured-decay example. The [notebook review](docs/notebook_review/NOTEBOOK_REVIEW.md) explains the supplied calculation, the input-format correction, and independent numerical checks.

### Time origin, normalization, and limits

Reconstructed time is the delay after excitation and the curves are the intrinsic model, not IRF-blurred count histograms. Reconvolution amplitudes use the model's excitation-origin convention. Tail-fit amplitudes refer to the first fitted tail channel; using them at excitation requires an explicitly known elapsed reference delay. The workbench needs that reference information before back-extrapolating a tail fit. Short components can become highly uncertain under such back-extrapolation.

Steady-state scaling provides relative spectral intensity unless the input spectrum and experimental factors are absolutely calibrated. With peak-normalized steady-state intensity, TRES values have relative-intensity-per-ns units. TRANES has inverse-nm units when area normalization is performed over wavelength in nm. It has unit area only over the supplied wavelength range.

Large reconstruction delays may lie beyond the measured kinetics and are model extrapolation. Five times the longest lifetime is a useful preview horizon, not a guarantee of negligible experimental signal. Spectral curve styles interpolate the calculated samples for display; they do not establish additional spectral resolution or a precise sub-grid peak position. Wavenumber conversion and its required density Jacobian are not implemented.

The spectral export uses wavelength rows and delay columns, matching the reconstructed-spectra layout. It includes both TRES and TRANES matrices, sampled steady-state values, fitted amplitudes used in reconstruction, and reconstruction/source metadata. This orientation differs from the raw fitting input, which has time rows and wavelength columns. See [DATA_FORMATS.md](DATA_FORMATS.md).

## Preparing your own files

The recommended TRES arrangement is one time bin per row, with time in column 1 and one decay per remaining column:

```csv
time_ns,500,520,540
0.000,2,0,1
0.025,1,3,2
0.050,0,1,1
```

These are illustrative counts. The first header cell is a label; the following header values are the actual wavelengths in nm. A separate shared IRF has two columns:

```csv
time_ns,counts
0.000,0
0.025,1
0.050,0
```

Use raw, nonnegative photon counts and retain zero-count channels. Select the calibrated time unit in the app; an ordinary CSV header does not silently override that selector. Shared IRF and decay bin widths must match. Preserve the measured IRF pulse and its relevant tail.

The app accepts tiny rounding differences in exported channel-time labels. It interprets the known fixed-width channel sequence within a tolerance of **0.0001 of one bin**, retains original labels in the import audit, and leaves every photon count in its original channel. This is not histogram rebinning. Missing channels or larger irregularities still require investigation.

Use UTF-8 text, a comma separator, a decimal point, and full numeric precision. For Origin worksheet export, put `time_ns` and the numeric wavelengths in the Long Name row; export that one header row, turn off extra short-name/unit/comment/index rows, and keep full precision enabled. Check the exported text in a plain-text editor. See [DATA_FORMATS.md](DATA_FORMATS.md) for complete examples, multi-curve assignment, missing endpoints, and file-specific mappings.

## What gets saved

The browser's download button produces a ZIP containing:

- `result.json`: parameters, settings, statistics, warnings, and local lifetime standard errors.
- `curves.csv`: observed counts, predicted counts, inclusion flags, and residuals.
- `amplitudes.csv`: component amplitudes and amplitude × lifetime.
- Import audit and source-curve mapping where applicable; wavelength mapping when supplied.

The reconstruction workspace has a separate spectral download containing `tres.csv`, `tranes.csv`, `sampled_steady_state.csv`, `normalized_steady_state.csv`, `fitted_amplitudes.csv`, `spectral_metadata.json`, and source metadata. Save both the fit result and the reconstruction result to preserve the complete analysis.

Keep the original measurement files alongside your result archive. A page reload or stopped Python process is not a substitute for downloading results. Loading files and fitting them does not overwrite the originals.

## Troubleshooting

| Symptom | What to check |
|---|---|
| `conda` is not recognized | Open Anaconda Prompt from Start, or the environment terminal from Navigator. |
| `python` is not recognized | Activate the installed environment in that terminal. |
| `No module named streamlit` | Check the active environment; use `python -m pip show streamlit` or `conda list streamlit` to inspect that same environment. |
| `app.py` cannot be found | Change to the inner extracted folder containing `app.py` and `requirements.txt`. |
| Browser shows an older page | Stop the old server with Ctrl+C and use the Local URL printed by the new process. |
| A required DLL is blocked by Application Control | Windows has blocked a compiled library. The Conda/MKL route uses an alternative official build but acceptance is not guaranteed. Record the blocked filename and obtain a trusted supported build; security-policy changes are not part of installation. |
| Wavelengths are requested | Supply a numeric wavelength header or one wavelength per decay column. Curve indices are not wavelengths. |
| Reconstruction asks for a steady-state spectrum | Supply a spectrum belonging to the fitted sample with wavelength coverage spanning every fitted wavelength. The measured example decays do not include that spectrum. |
| Reconstruction amplitude/lifetime mismatch | The global-fit route retains this pairing automatically. For optional external results, enter one lifetime per amplitude column in the same order. |
| Normalized wavelengths do not match | Supply one normalized row per global-fit wavelength, or choose the full-spectrum path to normalize and sample a denser spectrum. |
| Time-axis error | Check the selected unit, duplicate/missing channels, and actual channel calibration. Tiny export rounding is handled automatically. |
| Native binary file rejected | Export a text histogram containing the numeric time axis and raw counts. Renaming a binary file does not convert it. |

## Methods, limits, and advanced use

Read [SCIENTIFIC_NOTES.md](SCIENTIFIC_NOTES.md) before interpreting multicomponent fits. The current models are nonnegative parallel exponentials with tail fitting or IRF reconvolution. Polarized VV/VH import is supported, but rotational-anisotropy and general signed kinetic models are not implemented. Local standard errors are conditional approximations and exclude IRF and calibration uncertainty.

The command line provides `fit`, `global`, `inspect`, `curve-fit`, and `curve-global`:

```sh
python -m tcspc --help
python -m tcspc curve-global --help
```

The Python API provides `FitConfig`, `fit_decay`, `fit_global`, and `load_plotbox`. The spectral module supplies `load_fitted_spectra`, `load_steady_state`, and `reconstruct_spectra`. Array order for global counts is `[time, trace]`; fitted amplitudes are `[trace, component]`.

To run the included checks:

```sh
python -m unittest discover -s tests -v
```

The release's numerical and interface checks are recorded in [VALIDATION.md](VALIDATION.md).

See `LICENSE` for the code license and [CHANGELOG.md](CHANGELOG.md) for release changes. Example measurements were supplied for development and are distinct from the synthetic known-lifetime validation data.
