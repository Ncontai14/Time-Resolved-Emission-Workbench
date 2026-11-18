"""Continue a global lifetime fit with steady-state-scaled TRES and TRANES."""
from __future__ import annotations

from io import BytesIO
import hashlib
import json
from pathlib import Path
import tempfile
import threading
import zipfile

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from .spectra import (SteadyStateSpectrum, load_fitted_spectra,
                      load_steady_state, reconstruct_spectra)
from .spectral_display import display_curve


@st.cache_resource(show_spinner=False)
def _local_plot_lock():
    return threading.RLock()


def _numbers(text, *, name, count=None):
    try:
        values = np.asarray([float(v) for v in text.replace(",", " ").split()])
    except ValueError as error:
        raise ValueError(f"{name} must be a comma- or space-separated list of numbers.") from error
    if not len(values) or not np.all(np.isfinite(values)):
        raise ValueError(f"Enter finite values for {name.lower()}.")
    if count is not None and len(values) != count:
        raise ValueError(f"Enter exactly {count} lifetimes, in the same order as the amplitude columns.")
    return values


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return str(value)


def _spectra_zip(result, source_metadata):
    with tempfile.TemporaryDirectory(prefix="tcspc_spectra_") as folder:
        result.export(folder)
        (Path(folder) / "reconstruction_source.json").write_text(
            json.dumps(source_metadata, indent=2, default=_json_default, allow_nan=False),
            encoding="utf-8")
        output = BytesIO()
        with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for file in sorted(Path(folder).iterdir()):
                if file.is_file():
                    archive.write(file, file.name)
        return output.getvalue()


def _plot_curve(ax, wavelengths, values, style, *, label=None, markersize=3):
    x, y = display_curve(wavelengths, values, style=style)
    line, = ax.plot(x, y, label=label)
    ax.plot(wavelengths, values, "o", ms=markersize, color=line.get_color(), label="_nolegend_")
    scale = max(float(np.max(np.abs(values))), np.finfo(float).tiny)
    return bool(np.min(y) < -1e-12*scale)


def render_spectra_workspace(fit_payload=None, plot_lock=None):
    """Use the latest global fit by default; external tables are an optional route."""
    lock = _local_plot_lock() if plot_lock is None else plot_lock() if callable(plot_lock) else plot_lock
    # A failed or incomplete new input must not leave an earlier result in use.
    st.session_state.pop("last_spectral_reconstruction", None)
    st.subheader("TRES and TRANES from your global fit")
    workflow = st.radio("Fitted-decay source", ["Latest global fit", "Example or external fitted table"],
                        key="spectra_workflow", horizontal=True)
    fit = None
    embedded_ss = None
    source_metadata = {"source_mode": workflow}
    try:
        if workflow == "Latest global fit":
            if fit_payload is None:
                st.info("Run a fit in Global TRES first, then select Continue to TRES / TRANES. Its wavelengths, component lifetimes and pre-exponential amplitudes are carried over automatically; only a steady-state fluorescence file is needed here.")
                return
            fit = fit_payload["result"]
            wavelengths = fit_payload.get("wavelengths")
            if wavelengths is None or len(wavelengths) < 2:
                st.info("The global result needs at least two confirmed wavelengths. Supply the wavelength mapping in Global TRES and run the fit again.")
                return
            amplitudes, lifetimes = fit.amplitudes, fit.lifetimes
            source_name = fit_payload.get("source", "Global fit")
            st.success(f"Connected to global fit: {source_name}")
            st.write(f"Using {len(wavelengths)} fitted wavelengths and {len(lifetimes)} lifetime component(s). No amplitude or lifetime upload is required.")
            st.write("Shared lifetimes (ns): " + ", ".join(f"{tau:.8g}" for tau in lifetimes))
            st.caption("This uses the latest completed global fit. Changing fit settings takes effect after Run fit is pressed again.")
            if not fit.diagnostics.get("converged", False):
                st.warning("The selected global fit did not converge. Review its fit diagnostics before interpreting the reconstructed spectra.")
            source_metadata.update({"source_name": source_name, "fit_parameters": fit.to_dict(),
                                    "fit_input_identity": fit_payload.get("data_id")})
            source_identity = str(fit_payload.get("data_id", source_name))
        else:
            source = st.selectbox("Example or fitted-table input", ["Three-component worked example", "Upload fitted amplitudes"], key="spectra_advanced_source")
            if source == "Three-component worked example":
                path = Path(__file__).resolve().parents[1] / "examples" / "fitted_three_components.csv"
                table = load_fitted_spectra(path)
                tau_default = "114.15, 17.415, 4.613"
                source_identity = path.name
                st.info("This separate worked example includes fitted amplitudes and its own steady-state intensities. For your measurements, use Latest global fit.")
            else:
                upload = st.file_uploader("Fitted amplitude table", type=["csv", "dat", "txt", "tsv"], key="spectra_fitted_upload")
                st.caption("Optional external route: wavelength_nm, A1, A2, … and an optional named final steady-state column.")
                if upload is None:
                    st.info("Upload a fitted table, or choose Latest global fit to use the workbench result directly.")
                    return
                table = load_fitted_spectra(BytesIO(upload.getvalue()))
                tau_default = ""
                source_identity = hashlib.sha256(upload.getvalue()).hexdigest()
            wavelengths, amplitudes, embedded_ss = table.wavelengths, table.amplitudes, table.steady_state
            source_metadata.update({"source_name": table.metadata.get("source_name"), "fitted_table_import": table.metadata})
            tau_text = st.text_input("Component lifetimes in amplitude-column order (ns)", value=tau_default,
                                     help="Keep the lifetime order paired with A1, A2, and subsequent columns.", key=f"spectra_taus_{source}")
            lifetimes = _numbers(tau_text, name="Lifetimes", count=amplitudes.shape[1])
            if np.any(lifetimes <= 0):
                raise ValueError("All component lifetimes must be positive.")

        with st.expander("Inspect inputs from the fit", expanded=False):
            frame = pd.DataFrame(amplitudes, columns=[f"A{k+1}" for k in range(amplitudes.shape[1])])
            frame.insert(0, "wavelength_nm", wavelengths)
            st.dataframe(frame, hide_index=True, width="stretch")
        input_key = hashlib.sha256((source_identity + repr(np.asarray(wavelengths).tolist())).encode()).hexdigest()[:16]
        steady_path = st.radio("Steady-state input", ["Normalized values at fitted wavelengths", "Full steady-state spectrum", "Enter normalized intensities manually"],
                               key=f"spectra_ss_path_{workflow}")
        manual = steady_path == "Enter normalized intensities manually"
        exact = manual or steady_path == "Normalized values at fitted wavelengths"
        use_embedded = not manual and embedded_ss is not None and st.checkbox("Use the steady-state column in this example/table", value=True, key=f"spectra_embedded_{input_key}")
        if manual:
            st.caption("Wavelengths come from the fit and are fixed. Type or paste one normalized intensity per row, then click Apply intensities. Values are used as entered; the largest value does not need to be 1.")
            saved_key = f"spectra_manual_values_{input_key}"
            initial = st.session_state.get(saved_key, [np.nan]*len(wavelengths))
            manual_frame = pd.DataFrame({"wavelength_nm": np.asarray(wavelengths),
                                         "normalized_intensity": initial})
            with st.form(f"spectra_manual_form_{input_key}"):
                entered = st.data_editor(manual_frame, hide_index=True, num_rows="fixed",
                    disabled=["wavelength_nm"], width="stretch", key=f"spectra_manual_editor_{input_key}",
                    column_config={"wavelength_nm": st.column_config.NumberColumn("Wavelength (nm)", format="%.6g"),
                                   "normalized_intensity": st.column_config.NumberColumn("Normalized fluorescence intensity", min_value=0.0, max_value=1.0, format="%.8g", required=True)})
                applied = st.form_submit_button("Apply intensities", type="primary")
            raw_values = entered["normalized_intensity"].tolist()
            if applied:
                st.session_state[saved_key] = raw_values
            missing = pd.isna(entered["normalized_intensity"])
            if missing.any():
                message = f"Enter a normalized intensity for every wavelength ({int(missing.sum())} blank row(s)), then click Apply intensities. Blank cells are not treated as zero."
                (st.error if applied else st.info)(message)
                return
            values = np.asarray(raw_values, dtype=float)
            if not np.all(np.isfinite(values)) or np.any(values < 0) or np.any(values > 1):
                raise ValueError("Manual intensities must be finite numbers between 0 and 1.")
            if not np.any(values > 0):
                raise ValueError("Enter at least one positive intensity; an all-zero spectrum cannot produce TRANES.")
            normalization, alignment = "none", "exact"
            steady = SteadyStateSpectrum(np.asarray(wavelengths), values,
                       {"source_name": "manually entered normalized intensities", "input_method": "manual table"})
            source_metadata["steady_state_source"] = "manual normalized intensity table"
            st.caption("Spectra use the last applied table. Click Apply intensities again after changing any entries. Entered values are included in the spectral download.")
        elif exact:
            st.caption("Upload two columns: wavelength_nm, normalized_intensity. Include exactly one row per fitted wavelength, in any order. Values are matched by wavelength and used directly; no rescaling or interpolation is applied.")
            normalization = "none"
            alignment = "exact"
            if use_embedded:
                steady = SteadyStateSpectrum(np.asarray(wavelengths), np.asarray(embedded_ss), {"source_name": "embedded steady-state column"})
            else:
                upload = st.file_uploader("Normalized steady-state intensities", type=["csv", "dat", "txt", "tsv"], key=f"spectra_exact_ss_{input_key}")
                if upload is None:
                    st.info("Upload the normalized steady-state values for these fitted wavelengths, or choose Full steady-state spectrum.")
                    return
                steady = load_steady_state(BytesIO(upload.getvalue()))
                source_metadata["steady_state_source"] = upload.name
            if np.any(steady.intensities < 0) or np.any(steady.intensities > 1+1e-9):
                raise ValueError("Normalized intensities must lie between 0 and 1. For unnormalized data, choose Full steady-state spectrum.")
        else:
            st.caption("Upload a two-column spectrum: wavelength_nm, intensity. Normalize the full spectrum first; values are then matched or linearly interpolated at the fitted wavelengths. The spectrum must cover all fitted wavelengths.")
            normalization = st.selectbox("Normalization method", ["peak", "minmax"],
                format_func={"peak": "Divide by maximum (I / max)", "minmax": "Scale minimum to 0 and maximum to 1"}.get,
                key="spectra_normalization_method")
            alignment = "interpolate"
            if normalization == "minmax":
                st.caption("Uses (I − minimum) / (maximum − minimum). This subtracts the minimum as a baseline before scaling.")
            if use_embedded:
                steady = SteadyStateSpectrum(np.asarray(wavelengths), np.asarray(embedded_ss), {"source_name": "embedded steady-state column"})
            else:
                upload = st.file_uploader("Steady-state fluorescence spectrum", type=["csv", "dat", "txt", "tsv"], key=f"spectra_full_ss_{input_key}")
                if upload is None:
                    st.info("Upload the steady-state fluorescence spectrum measured for this sample.")
                    return
                steady = load_steady_state(BytesIO(upload.getvalue()), allow_negative=normalization == "minmax")
                source_metadata["steady_state_source"] = upload.name
        if use_embedded:
            source_metadata["steady_state_source"] = "embedded fitted-table column"
        source_metadata.update({"steady_state_path": steady_path, "steady_state_import": steady.metadata})

        reference_delay = 0.0
        if fit is not None and fit.config.method == "tail":
            st.warning("Tail-fit amplitudes refer to the first fitted channel. Supply its known elapsed delay after excitation to reconstruct excitation-origin spectra.")
            if not st.checkbox("I know the tail reference delay after excitation", key="spectra_tail_confirm"):
                return
            reference_delay = st.number_input("Tail reference delay after excitation (ns)", min_value=0.0, value=0.0, format="%.8g", key="spectra_tail_delay")
        elif fit is None:
            origin = st.selectbox("Amplitude reference", ["At excitation", "At a later tail reference"], key="spectra_amplitude_origin")
            if origin == "At a later tail reference":
                reference_delay = st.number_input("Tail reference delay after excitation (ns)", min_value=0.0, value=0.0, format="%.8g", key="spectra_external_tail_delay")

        default_delays = (", ".join(f"{float(np.max(lifetimes))*fraction:.8g}" for fraction in [0, 0.25, 0.5, 1, 2])
                          if fit is not None else "0.025, 1, 10, 50, 100")
        delays = _numbers(st.text_input("Reconstruction delays (ns)", value=default_delays,
                                        key=f"spectra_times_{input_key}"), name="Reconstruction delays")
        if np.any(delays < 0) or len(np.unique(delays)) != len(delays) or len(delays) > 64:
            raise ValueError("Enter up to 64 unique, nonnegative reconstruction delays.")
        delays = np.sort(delays)
        styles = {"pchip": "Smooth curve (PCHIP)", "linear": "Straight lines", "bspline": "Cubic B-spline"}
        style = st.selectbox("Spectral curve style", list(styles), format_func=styles.get, key="spectra_curve_style")
        if style == "bspline":
            st.caption("B-splines can overshoot between measured wavelengths. Markers and exported matrices retain the calculated samples.")
            if len(wavelengths) < 4:
                st.caption(f"With {len(wavelengths)} fitted wavelengths, the B-spline uses degree {len(wavelengths)-1}.")
        preview_max = st.number_input("Maximum preview delay (ns)", min_value=0.001, value=max(float(delays[-1]), float(5*np.max(lifetimes)), 1.0), format="%.8g", key=f"spectra_preview_max_{workflow}")
        preview_delay = st.slider("Preview delay after excitation (ns)", min_value=0.0, max_value=float(preview_max), value=min(1.0, float(preview_max)), step=min(0.001, float(preview_max)/1000), format="%.5g", key=f"spectra_preview_delay_{workflow}")
        options = dict(steady_state=steady, steady_normalization=normalization,
                       steady_alignment=alignment, amplitude_reference_delay_ns=float(reference_delay))
        result = reconstruct_spectra(wavelengths, amplitudes, lifetimes, delays, **options)
        preview = reconstruct_spectra(wavelengths, amplitudes, lifetimes, [preview_delay], **options)
        result.metadata.update({"steady_state_input_path": steady_path, "display_curve_style": style,
                                "display_curve_points": 500, "sampled_spectra_modified_by_curve_style": False,
                                "display_curve_area_convention": "TRANES unit area applies to trapezoidal integration of fitted-wavelength samples; smooth display curves are not renormalized"})
        matched = pd.DataFrame({"fitted_wavelength_nm": result.wavelengths, "normalized_fluorescence_intensity": result.steady_state})
        st.subheader("Steady-state intensities matched to the global fit")
        st.caption("Matched directly by wavelength; supplied normalized values retained." if exact else "Full spectrum normalized before sampling at the fitted wavelengths.")
        st.dataframe(matched, hide_index=True, width="stretch")
        with st.expander("Inspect full steady-state spectrum and normalization", expanded=False):
            source_waves = np.asarray(result.metadata["steady_state_source_wavelengths_nm"])
            source_raw = np.asarray(result.metadata["steady_state_source_intensities"])
            source_normalized = np.asarray(result.metadata["steady_state_source_normalized_intensities"])
            st.dataframe(pd.DataFrame({"wavelength_nm": source_waves, "supplied_intensity": source_raw,
                                       "normalized_intensity": source_normalized}), hide_index=True, width="stretch")
            with lock:
                fig, ax = plt.subplots(figsize=(9, 3))
                ax.plot(source_waves, source_normalized, color="0.5", lw=1, label="Supplied normalized spectrum")
                ax.plot(result.wavelengths, result.steady_state, "o", label="Values matched to global fit")
                ax.set(xlabel="Wavelength (nm)", ylabel="Normalized fluorescence intensity")
                ax.legend(fontsize="small");fig.tight_layout();st.pyplot(fig, width="stretch");plt.close(fig)
        overshoot = False
        with lock:
            fig, axes = plt.subplots(1, 2, figsize=(12, 3.8))
            for ax, values, title, ylabel in [(axes[0], preview.tres[:, 0], "TRES", "Relative intensity / ns"),
                                              (axes[1], preview.tranes[:, 0], "TRANES", "Area-normalized intensity (1/nm)")]:
                overshoot |= _plot_curve(ax, preview.wavelengths, values, style, markersize=4)
                ax.set(xlabel="Wavelength (nm)", ylabel=ylabel, title=f"{title} at {preview_delay:.5g} ns")
            fig.tight_layout();st.pyplot(fig, width="stretch");plt.close(fig)
        with lock:
            fig, axes = plt.subplots(1, 2, figsize=(12, 4.2))
            for j, delay in enumerate(result.times):
                overshoot |= _plot_curve(axes[0], result.wavelengths, result.tres[:, j], style, label=f"{delay:.6g} ns")
                overshoot |= _plot_curve(axes[1], result.wavelengths, result.tranes[:, j], style, label=f"{delay:.6g} ns")
            axes[0].set(xlabel="Wavelength (nm)", ylabel="Relative intensity / ns", title="Reconstructed TRES")
            axes[1].set(xlabel="Wavelength (nm)", ylabel="Area-normalized intensity (1/nm)", title="TRANES")
            for ax in axes:ax.legend(fontsize="x-small", ncol=2)
            fig.tight_layout();st.pyplot(fig, width="stretch");plt.close(fig)
        if overshoot:
            st.warning("The selected B-spline dips below zero between wavelength samples. This is a display interpolation artifact; choose PCHIP to avoid it. Calculated samples are unchanged.")
        st.caption("Curve styles change the visual connections only. Markers and exported matrices use the fitted wavelengths. TRANES has unit trapezoidal area on that wavelength grid; a smooth display curve may have a different area.")
        area = np.trapezoid(result.tranes, result.wavelengths, axis=0) if hasattr(np, "trapezoid") else np.trapz(result.tranes, result.wavelengths, axis=0)
        st.dataframe(pd.DataFrame({"delay_ns": result.times, "TRANES_spectral_area": area}), hide_index=True, width="stretch")
        st.caption("These are relative, intrinsic model spectra. Reconstruction outside the measured post-excitation time window depends on model extrapolation; fit and steady-state uncertainties are not propagated.")
        with st.expander("Calculation and reconstruction metadata", expanded=False):
            st.latex(r"I(\lambda,t)=S(\lambda)\frac{\sum_k A_k(\lambda)e^{-t/\tau_k}}{\sum_k A_k(\lambda)\tau_k},\qquad \mathrm{TRANES}(\lambda,t)=\frac{I(\lambda,t)}{\int I(\lambda,t)\,d\lambda}")
            st.caption("The decay-area denominator ensures that the time integral at each wavelength equals its normalized steady-state intensity. Constant background and IRF convolution are excluded from reconstruction.")
            st.json(result.metadata)
        st.download_button("Download TRES + TRANES matrices", data=_spectra_zip(result, source_metadata), file_name="tres_tranes_reconstruction.zip", mime="application/zip", width="content", key="spectra_download", on_click="ignore")
        st.session_state["last_spectral_reconstruction"] = result
    except (ValueError, OSError, RuntimeError) as error:
        st.error(str(error))
