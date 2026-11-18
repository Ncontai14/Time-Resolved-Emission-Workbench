# Windows quick start — version 0.1.0

For the full experimental-user operating manual, open [USER_GUIDE.md](USER_GUIDE.md) or `USER_GUIDE.html`.

Use the entire extracted package. Find the folder containing `app.py`, `requirements.txt`, and `tcspc`.

## Already have a working environment?

Use it if the required dependencies are installed. Stop any running app with **Ctrl+C**, activate that environment, and launch `app.py` from this extracted folder. The installation steps below are for a first-time user.

## First installation with Anaconda

Open **Anaconda Prompt** from Start. Create the environment with this single-line command:

```bat
conda create -n tcspc_workbench --override-channels -c defaults python=3.12 numpy scipy matplotlib pandas threadpoolctl "streamlit>=1.51,<2" "blas=*=mkl"
```

Enter `y` when prompted. This MKL-specific command is for Windows on Intel/AMD 64-bit processors.

Activate and check it:

```bat
conda activate tcspc_workbench
python -c "import numpy, scipy.optimize, scipy.signal, matplotlib.pyplot, pandas, streamlit, threadpoolctl; print('Dependencies OK')"
```

Change to your actual program folder. Replace the example path:

```bat
cd /d "C:\Users\YourName\Downloads\TCSPC_Workbench"
python -m streamlit run app.py
```

Open the Local URL printed in the terminal, usually `http://localhost:8501`. Keep the terminal open. Press **Ctrl+C** to stop the app.

The Conda command installs all required dependencies. Use the separate pip/virtual-environment route in README only if you choose that installation method instead.

## Every later session

```bat
conda activate tcspc_workbench
cd /d "C:\Users\YourName\Downloads\TCSPC_Workbench"
python -m streamlit run app.py
```

Use your actual folder path. Once that environment is activated, `START_WORKBENCH_WINDOWS.bat` is an optional launch shortcut from the same terminal.

## Try the normal single-decay window first

Choose **Single decay → Bundled decay + IRF**. Use one component, IRF reconvolution and Poisson likelihood. In **Background, shift, and excitation**, fit the decay background, subtract an IRF background of `1.4` counts/bin, limit the IRF pulse to `3.8–8` ns, and fit its shared shift. These settings belong only to this supplied example.

In **Kinetics preview and fit interval**, move the two ends of **Fit interval handles (time bins)** or enter **First fitted time (ns)** and **Last fitted time (ns)**. The shaded region and guides show the channels selected for fitting. **Kinetics preview view → Selected fit window** enlarges that region; **Whole record** shows its position in the acquisition. These controls are available for both a normal single decay and a global dataset.

Press **Run fit** after changing the interval. The previous result remains until a new fit runs. Using the whole imported record and the example settings above gives approximately `2.549 ns`. Inspect residuals and download the fit-result ZIP before closing or refreshing.

## Try a measured TRES example

Choose **Global TRES → Upload ASCII data**. Upload `examples/measured_TRES.csv` and `examples/measured_IRF.csv`. Select **ns** and leave wavelengths blank because they are in the CSV header.

Use two components, initial lifetimes `2.5, 4.0`, Poisson likelihood, fitted backgrounds, IRF baseline `0.38`, IRF range `3.2–6.5` ns, and a fitted shared shift with bounds `−0.5, 0.5` ns. Leave periodic correction off. Inspect residuals and download the result ZIP.

Before running a global fit, use the linked fit-range slider and start/end boxes to select the channels entering the objective. Inspect the individual kinetics preview and shaded fitting region; a plot zoom alone is not a different fit range.

## Continue from the global fit to TRES/TRANES

The saved global fit supplies wavelengths, amplitudes, and matching lifetime components automatically. Click **Continue to TRES / TRANES** and keep **Fitted-decay source → Latest global fit**. Provide only the steady-state fluorescence information. Choose one of three paths:

- **Normalized values at fitted wavelengths**: one row for every fitted wavelength; values between 0 and 1 are matched and used directly.
- **Full steady-state spectrum**: upload wavelength versus fluorescence intensity, choose division by its maximum or 0–1 min–max normalization, then inspect the normalized spectrum sampled at the fitted wavelengths.
- **Enter normalized intensities manually**: fill the intensity column beside the fixed wavelengths and click Apply intensities; no file is needed.

### Enter normalized intensities directly

Choose **Steady-state input → Enter normalized intensities manually**. The table already contains the fitted wavelengths; that column is locked. Type or paste normalized fluorescence intensities into the intensity column and click **Apply intensities**. Fill every row with a finite value between 0 and 1; zeros are allowed, but the whole spectrum cannot be zero. Blank cells stay blank and must be completed.

Values are used directly, without interpolation or renormalization. Their maximum can be below 1 if the fitted wavelengths exclude the peak of a previously normalized spectrum. After editing a value, click **Apply intensities** again to update the spectra. The spectral ZIP records the entered values in `sampled_steady_state.csv` and identifies manual input in its metadata. No steady-state file upload is needed for this path.

The full spectrum must cover the fitted wavelength range. Inspect the matching table, enter reconstruction delays, and choose **Spectral curve style → Smooth curve (PCHIP)**, **Straight lines**, or **Cubic B-spline**. Plot style does not change the numerical reconstruction samples or exported matrices. A B-spline can overshoot between points.

The measured 51-wavelength example has no supplied matching steady-state spectrum. The optional three-component worked example has its own fitted amplitudes and embedded spectrum. Its A1, A2, A3 lifetimes are **114.15, 17.415, 4.613 ns**, belonging only to that separate example.

For your own files and troubleshooting, read [INSTALLATION.md](INSTALLATION.md) and [DATA_FORMATS.md](DATA_FORMATS.md). If a library import is blocked by Windows, record the exact error and blocked filename; the alternative Conda build is a test, not a guarantee of security-policy acceptance.
