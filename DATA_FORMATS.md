# Data formats and import

Time-Resolved Emission Workbench reads plain text tables and multi-curve ASCII exports. A file can contain a TCSPC histogram or a compatible emission intensity scan sampled at uniform time intervals. Recommended extensions are `.csv`, `.txt`, `.tsv`, or `.dat`.

## Decay and shared IRF tables

A plain decay or shared IRF has time in its first column and the measured signal in its second. A TCSPC example uses raw counts:

```csv
time_ns,counts
0.000,0
0.025,2
0.050,1
0.075,0
```

Use a decimal point. UTF-8 is recommended; plain ASCII characters are also valid UTF-8. CSV uses commas. The text reader can also handle tabs, semicolons, or whitespace when used consistently.

Choose the calibrated input unit explicitly in the app. Times and lifetimes are converted to ns internally. A plain header such as `time_ns` is a label rather than an automatic override of the selected unit.

Raw counts are nonnegative and include empty channels. Poisson analysis expects photon counts before smoothing, normalization, or background subtraction. Background is represented in the fitted model.

For a compatible non-counting measurement such as fluorescence upconversion, put the finite nonnegative intensity in the signal column. Choose unweighted least squares when approximately constant noise variance is a reasonable assumption. Keep the actual delay schedule when checking compatibility with the required uniform time grid. For reconvolution, place the measured temporal response in a separate file on a grid with the same spacing and time units.

## Raw time-by-wavelength fitting matrix

Rows are equally spaced time samples. Columns after time are wavelength-specific decays:

```csv
time_ns,500,520,540
0.000,1,2,0
0.025,0,1,3
0.050,2,0,1
0.075,1,1,0
```

Values shown are illustrative. The header values `500,520,540` are wavelengths in nm. Every count column needs a corresponding known wavelength. Alternatively, give the wavelength list in the app in exactly the count-column order.

A matrix with 51 wavelength decays has 52 CSV columns, including time. An IRF is a separate file, not a wavelength decay column. The Python global API uses counts shaped `[time, trace]`; amplitudes are returned as `[trace, component]`.

The packaged `measured_TRES.csv` has 1,999 time bins and 51 wavelengths: 490–590 nm in 2 nm steps. `measured_IRF.csv` is the matching shared IRF. The contributed example mapping is known; generic files must carry their own wavelength metadata.

## Origin export settings

Use a worksheet with time in column A and one decay in each subsequent column. Put `time_ns` in A's Long Name cell, then each actual numeric wavelength in the corresponding Long Name cell.

For Text/CSV export, use:

| Setting | Value |
|---|---|
| Encoding | UTF-8 |
| Separator | Comma |
| Decimal separator | Point |
| Include Long Name | On |
| Include Short Name, Units, Comments | Off for a single header row |
| Include row labels and indices | Off |
| Export with Full Precision | On |

Include the whole intended data range. Check selected-only export and hidden columns if a decay disappears. Inspect the exported file in a text editor: the wavelength labels must actually be present, not just visible as worksheet properties.

High precision prevents additional rounding; it cannot undo rounding already stored in imported coordinates. The app handles small timestamp rounding automatically. There is normally no need to regenerate a calibrated TCSPC time column before export.

## Channel-time rounding

An instrument may use equal-width channels while its printed time labels have small rounding differences. The importer checks whether the full axis agrees with one fixed-width channel sequence to within 0.0001 of a bin. When it does, the calculation uses that sequence and records the original labels and maximum adjustment.

Each observed value remains at its original channel. A missing internal channel or genuinely uneven acquisition times are rejected at import, including for tail fitting. Use the original acquisition settings to distinguish an export-rounding issue from a deliberate variable-step scan or a missing channel.

## Multi-curve ASCII bundles

Supported bundles begin with a marker and an explicit time-unit header, for example:

```text
PlotBox Export
Time[s]X	crv[0] []	crv[1] []
0	0	2
2.5e-11	1	3
5.0e-11	0	1
7.5e-11	0	2
```

Here `\t` represents a tab separator. The real file uses tab characters. `Time[s]X` explicitly declares seconds, so the bundle importer converts to ns automatically. Generic `crv[index]` names do not identify which curve is an IRF, wavelength, polarization channel, or laboratory time point.

In the interface, select the decay curve(s) and shared IRF explicitly. Indices are zero-based and exclude the time column. A range such as `1-51` includes both endpoints. A separate wavelength list is optional for an unlabeled series; its plots use curve indices when physical wavelengths are not supplied.

The command line supports inspection and selected-curve fits:

```sh
python -m tcspc inspect --input examples/curve_exports/TRES_IRF.dat
python -m tcspc curve-fit --input examples/curve_exports/decay_IRF.dat --decay-curve 1 --irf-curve 0 --irf-background 1.4 --irf-start 3.8 --irf-end 8 --fit-shift --out my_fit
python -m tcspc curve-global --input examples/curve_exports/TRES_IRF.dat --decay-curves 1-51 --irf-curve 0 --components 2 --tau 2.5 4.0 --irf-background 0.38 --irf-start 3.2 --irf-end 6.5 --fit-shift --out my_global
```

## Missing values and acquisition windows

In a bundle, `-` means absent data, not zero photon counts. It remains missing during import. Selecting curves can trim incomplete endpoints to retain their common complete interval. Missing values inside that interval are rejected. Endpoint rows with unknown timestamps are excluded and audited; times are not invented for them.

Plain tables require complete numeric observations. Resolve missing data explicitly before making a plain table; changing a missing value to zero changes the observation and its noise model.

The single-decay and global kinetics previews share linked channel handles, numerical first/last fitted times and display zoom controls. A selected display zoom does not crop the imported table; the chosen fitting interval is recorded in the fit settings.

Preserve the original full measurement and record any selected common interval. The fit objective uses only the chosen decay window, but reconvolution is evaluated using the available full imported excitation history before that window is applied.

The repetition-period model needs a complete, uniformly spaced, nonduplicated period: `number_of_bins × bin_width = period`. An exported array may contain a duplicate endpoint or lose final bins through missingness. Its displayed last timestamp is not enough to establish a complete period.

## Contributed raw example assignments

These assignments apply to the named files supplied with the package; the generic importer does not guess them from a filename.

| File in `examples/curve_exports/` | Zero-based curve roles |
|---|---|
| `decay_IRF.dat` | 0 IRF; 1 fluorescence decay |
| `decay.dat` | 0 long-lived RuBiPy decay; suitable for a tail fit |
| `TRES_IRF.dat` | 0 IRF; 1–51 decays at 490–590 nm in 2 nm steps |
| `decaytimecourse_IRF.dat` | 0 IRF; 1–10 repeated decays, nominal laboratory times 0, 15, …, 135 s |
| `anisodec.dat` | 0 VV; 1 VH |
| `anisodec_IRF.dat` | 0 IRF; 1 VV; 2 VH |
| `anisodec_IRF_magic.dat` | 0 IRF; 1 VV; 2 magic angle; 3 VH |

VV is parallel polarization and VH perpendicular polarization. The magic-angle intensity is curve 2 in the four-curve file. Coupled rotational-anisotropy fitting is not implemented, and a calibrated G-factor was not supplied. The polarized curves remain useful for inspecting the source data; a simple intensity fit of them is not a rotational-anisotropy measurement.

The repeated-decay laboratory times are nominal sequence times and exclude execution overhead. They are separate from the ns photon-delay axis and must not be labeled as wavelengths.

Native binary measurement/workspace formats are not imported. Use text exports exposing the numeric curves, units, and metadata. Keep the original acquisitions for provenance.

## Steady-state inputs for global-fit reconstruction

The saved global result supplies wavelength labels, pre-exponential amplitudes, and component lifetimes automatically. Supply the steady-state intensities in the manual table or a two-column numeric fluorescence file, for example:

```csv
wavelength_nm,intensity
490,0.12
500,0.35
510,0.80
520,1.00
```

These values are illustrative. Use wavelength in nm, UTF-8 text, decimal points, and a consistent comma, tab, or whitespace separator. Wavelengths must be positive, finite, and unique. Sorting wavelengths also sorts their paired intensities. The source measurement file remains unchanged.

### Already-normalized values at fitted wavelengths

Supply exactly one row for each global-fit wavelength, with finite normalized intensity from 0 to 1 and some positive signal. Rows can arrive in a different order. The workbench matches them to the global result and uses the intensity values directly. It does not renormalize, interpolate, or infer missing rows in this path. Missing wavelengths, duplicate rows, and unrelated extra wavelengths are rejected. Tiny coordinate-rounding differences are accepted only when matching is unambiguous.

A subset of a previously normalized emission spectrum may have a largest value below 1. That is allowed; the workbench does not force that subset's maximum to 1.

### Full steady-state spectrum

Supply the measured fluorescence spectrum with wavelength coverage spanning every global-fit wavelength. It may contain a denser or irregular wavelength grid. Choose one normalization:

| Method | Calculation on the full uploaded spectrum |
|---|---|
| Divide by maximum | `I / max(I)`; requires nonnegative intensities and a positive maximum. |
| 0–1 min–max | `(I - min(I)) / (max(I) - min(I))`; allows a finite negative raw offset but requires a nonconstant spectrum. |

Normalization precedes sampling at the fitted wavelengths. The workbench then uses linear interpolation inside the supplied coverage. It rejects extrapolation beyond the measured spectrum. A spectrum's maximum does not need to coincide with any global-fit wavelength.

Min–max normalization subtracts a constant before scaling. This changes the relative spectral shape and reconstructed TRANES when the original minimum is not zero. It is a selected normalization convention, not an automatic experimental baseline correction. Record any prior blank or response correction with the original measurement.

The preview and export record the normalized source spectrum, selected normalization, offset, divisor, wavelength-matching method, and intensities used at the fitted wavelengths. This fluorescence spectrum is separate from the raw photon histograms used for lifetime fitting.

## Compact fitted-amplitude table

This optional input supports externally fitted results and worked examples. A current global fit already provides these quantities and requires no fitted-table upload. The compact table is for reconstruction, not for refitting photon histograms. Each row is one wavelength:

```csv
wavelength_nm,A1,A2,A3,SteadyState
498,0,0,11239,0.0737495452699143
515,0,0,11608,0.715823095325772
532,10.17,378.9,12063,0.465667354828892
```

The first three rows above are from the contributed three-component example; they are not a complete spectrum. The complete file is `examples/fitted_three_components.csv`. It uses the original named final steady-state column, which the importer recognizes. Enter `114.15, 17.415, 4.613` ns for its A1, A2, A3 columns respectively. The file has 13 wavelengths from 498 through 653 nm; their spacing is not uniform, which is allowed for spectral integration.

A final steady-state column is optional. With a header, name it clearly, for example `SteadyState`, `SS`, or `Fluo_intensity_steadystate`; the remaining columns after wavelength are component amplitudes. Without embedded steady-state values, upload a separate two-column spectrum. In the browser, use a named header when the final column contains steady-state intensity. A headerless table treats all columns after wavelength as amplitudes; the Python API can instead explicitly set `has_steady_state=True` for that compact layout.

Use CSV, tab-separated text, or whitespace-delimited numeric text consistently, with decimal points and UTF-8 encoding. The uploaded reference notebook used a tab-only, no-header reader, while the supplied fitted table is comma-separated and has a header. The workbench's dedicated fitted-table importer handles that format difference. Wavelengths are finite, positive, unique values in nm; amplitudes are finite and nonnegative.

A normal fit export `amplitudes.csv` uses a long trace/component layout. It is not this compact wide table. Continue directly from the current global result in the interface, or explicitly reshape external results into the compact layout and retain the wavelength mapping.

Lifetime values must match the amplitude-column order. Sorting wavelengths reorders whole rows, including their amplitudes and intensities. It must not reorder component lifetimes independently. Current global results already retain the correct amplitude/lifetime pairing.

## Reconstruction output matrices

The spectral download contains `tres.csv` and `tranes.csv`. Both place **wavelengths in rows and reconstruction delays in columns**, for example:

```csv
wavelength_nm,0.025_ns,1_ns,10_ns
498,0.012,0.010,0.002
515,0.041,0.034,0.006
```

The numeric intensities here are schematic, not reference results. Delays are nonnegative ns after the declared excitation origin. Metadata records their exact values, lifetime order, normalization, amplitude reference delay, and any model-extrapolation qualifications. The files represent reconstructed intrinsic spectra, not newly measured count data.

`sampled_steady_state.csv` retains the steady-state values used at the reconstruction wavelengths. `normalized_steady_state.csv` retains the complete supplied spectrum after the selected normalization, before sampling. Both steady-state CSV files have two columns (wavelength and normalized intensity) so they can be read again as spectrum inputs. Original supplied intensities, normalization offset and divisor remain in `spectral_metadata.json`. `fitted_amplitudes.csv` retains the component amplitudes in the recorded reference-time convention; the reference delay is recorded in metadata. `spectral_metadata.json` records the mathematical conventions, steady-state normalization and matching, and scale. The browser ZIP also includes `reconstruction_source.json` with the selected input provenance and curve-style setting. Smooth plot connections are display interpolation; numerical matrices retain the original fitted wavelength samples and their trapezoidal TRANES area normalization. Keep these together with the lifetime-fit export and original spectrum. The output orientation is intentionally different from the raw fitting matrix; relabeling a pre-exponential-factor plot does not produce TRES.

For tail fits, supply the known delay of the amplitude reference time after excitation. The conversion to excitation-origin amplitudes is `A_at_excitation = A_at_reference × exp(reference_delay / tau)`. The original amplitude table remains unchanged. A positive steady-state intensity with zero fitted decay area is inconsistent and cannot be reconstructed; a spectrum with zero area cannot be normalized to TRANES.

## Manual steady-state table

In TRES/TRANES, **Enter normalized intensities manually** requires no file. Fit wavelengths are fixed; enter one normalized intensity (0–1) per row and select **Apply intensities**. Every row must be filled and at least one intensity must be positive. Values retain their entered scale. `sampled_steady_state.csv` and the reconstruction metadata preserve these inputs.
