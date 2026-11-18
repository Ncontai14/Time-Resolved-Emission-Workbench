# TCSPC Workbench

A local Python application for fluorescence lifetime fitting, global TCSPC analysis, and reconstruction of time-resolved emission spectra (TRES) and area-normalized spectra (TRANES).

**Application version: 0.1.0 — initial public release.** The repository includes the local application, numerical examples, operating guides, tests and contribution tools.

![Single-decay fitting interval with linked time fields and range handles](docs/interface_single_range.png)

## What it does

- Fits single decays by measured-IRF reconvolution or exponential tail fitting.
- Fits multiple wavelength decays with shared lifetimes and separate nonnegative amplitudes/backgrounds.
- Provides linked fit-window handles, numerical time limits and kinetics previews.
- Reconstructs TRES/TRANES using the global-fit result and steady-state fluorescence intensities.
- Accepts normalized intensity files, full spectra with normalization/sampling, or manually entered intensities at the fitted wavelengths.
- Offers PCHIP, straight-line and B-spline display connections while preserving calculated spectral samples.
- Exports numerical results, settings, diagnostics, wavelength mappings and reconstruction provenance.

Calculations run locally; Streamlit displays the interface in your browser. Ordinary local analysis needs no cloud account. The application includes a Python API and command-line tools.

## Start here

**Experimental users:** read the illustrated [User guide](USER_GUIDE.md). Its offline edition, `USER_GUIDE.html`, opens directly in a browser after extracting the whole package.

**Already have a compatible Python environment?** You can use it. Open a terminal with that environment activated, change to this repository's folder containing `app.py`, then run:

```bash
python -m streamlit run app.py
```

Use the Local URL printed in the terminal, usually `http://localhost:8501`. Keep that terminal open; Ctrl+C stops the app.

**First installation:** use the [Windows quick start](QUICK_START_WINDOWS.md) or the complete [installation instructions](INSTALLATION.md). The Windows guide includes an Anaconda/MKL option. A standard virtual-environment route is:

```bash
python -m venv .venv
```

Activate it with `source .venv/bin/activate` on macOS/Linux or `.venv\Scripts\activate.bat` in Windows Command Prompt, then:

```bash
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

Python 3.12 is the version used for local release validation; the declared minimum is 3.10. Dependency installation requires internet access. For Anaconda Navigator, PowerShell or an existing environment, follow [INSTALLATION.md](INSTALLATION.md) rather than mixing installation routes.

## Typical workflow

1. Load a decay and measured IRF, or a time-by-wavelength decay matrix.
2. Confirm units, component count, IRF preparation and the fitting interval.
3. Run the fit and inspect residuals, convergence and parameter stability.
4. For global spectral analysis, select **Continue to TRES / TRANES**. The wavelengths, amplitudes and lifetimes transfer automatically.
5. Supply steady-state intensities by file or manual table, choose reconstruction delays and download both the fit and reconstruction archives.

Use original photon counts for Poisson fitting. The [data-format guide](DATA_FORMATS.md) shows CSV layouts, units and source-curve selection. User measurements belong outside the tracked examples, for example in the ignored `data/` directory.

## Scientific scope and current limitations

- IRF and decay **bin widths must match** for reconvolution. Different acquisition windows can be supported; unequal bin widths are not automatically rebinned.
- One time-unit selector applies to both plain-text decay and IRF uploads. Convert them to the same units before import.
- Tail-fit reconstruction needs the known elapsed delay between excitation and the first fitted channel.
- The model uses nonnegative exponential components; signed rise components, dedicated anisotropy models and kinetic target analysis are not implemented.
- TRES reconstruction needs corresponding steady-state intensities. The bundled 51-wavelength measured example does not include its matching steady-state spectrum.
- Reconstruction uncertainties are not propagated. Spline connections are for display and do not establish finer measured spectral resolution.

See [Scientific notes](SCIENTIFIC_NOTES.md) for definitions and assumptions. A converged optimizer or a diagnostic near one does not establish a physical model.

## Validation and development

The v0.1.0 numerical/import/CLI suite contains **69 tests**. Additional release checks exercised manual table entry and the single/global interfaces in Streamlit and Chromium; details are in [VALIDATION.md](VALIDATION.md).

Run the numerical suite from the repository root:

```bash
python -m unittest discover -s tests -v
```

GitHub Actions is configured to run the suite on Linux and Windows with Python 3.10 and 3.12. Results of hosted runs become available after the repository is uploaded and Actions runs; no hosted status is claimed here.

[Contributing](CONTRIBUTING.md) explains development and numerical-change expectations. [Support](SUPPORT.md) describes useful bug reports. Keep private measurements, credentials and local environments out of commits; `.gitignore` covers designated data/output folders, not arbitrary files everywhere on the computer.

## Repository contents

| Path | Purpose |
|---|---|
| `app.py`, `tcspc/` | Application, numerical engine, importers and CLI |
| `tests/` | Numerical, import, export and CLI tests |
| `examples/` | Synthetic and contributed examples, with provenance |
| `example_results/` | Reference outputs reproducible with example scripts |
| `docs/` and root guides | Screenshots, scientific review and offline documentation |
| `.github/` | CI workflow and issue/pull-request templates |

## Licence, data and citation

Original project code/documentation use the [MIT licence](LICENSE). Contributed measured examples retain their original rights; the code licence does not grant rights to those measurements. Read [Example-data provenance](examples/README.md) before redistributing the full example collection publicly.

Use [CITATION.cff](CITATION.cff) for software citation. It currently credits the project contributors collectively; the maintainer should enter the intended author details before public release. No paper, repository URL or DOI is invented. A versioned repository can later be archived on Zenodo.

**Preparing your own GitHub upload?** Follow [GITHUB_UPLOAD.md](GITHUB_UPLOAD.md), including what to personalize and how to create the v0.1.0 release.
