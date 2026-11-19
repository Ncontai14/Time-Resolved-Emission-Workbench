"""Local Streamlit interface. Run: streamlit run app.py"""
from __future__ import annotations

from io import BytesIO
from pathlib import Path
import hashlib
import json
import tempfile
import threading
from types import SimpleNamespace
import zipfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import pandas as pd
import streamlit as st

from tcspc import FitConfig, fit_decay, fit_global
from tcspc._release import WORKBENCH_NAME, __version__
from tcspc.io import load_decay, load_tres
from tcspc.plotbox import load_plotbox
from tcspc.cli import parse_curve_indices
from tcspc.spectra_ui import render_spectra_workspace

ROOT = Path(__file__).resolve().parent
st.set_page_config(page_title=WORKBENCH_NAME, page_icon="🔬", layout="wide")
st.title(WORKBENCH_NAME)
st.write("Fit measured decays, choose the fitting interval visually, and reconstruct TRES/TRANES from steady-state spectra.")
st.caption(f"Version {__version__} · local analysis · all fitted times and lifetimes are reported in ns")


@st.cache_resource(show_spinner=False)
def _plot_lock():
    # Matplotlib has shared parser/font state; Streamlit reruns can overlap.
    return threading.RLock()


def parse_lifetimes(text: str, components: int) -> list[float] | None:
    if not text.strip():
        return None
    values = [float(item) for item in text.replace(",", " ").split()]
    if len(values) != components or any(not np.isfinite(item) or item <= 0 for item in values):
        raise ValueError(f"Enter {components} positive lifetime value(s) in ns.")
    return values


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    return str(value)


def is_plotbox_content(content: bytes) -> bool:
    return content.lstrip(b"\xef\xbb\xbf\r\n \t").startswith(b"PlotBox Export")


def parse_wavelengths(text: str, expected: int):
    if not text.strip():
        return None
    values = np.asarray([float(item) for item in text.replace(",", " ").split()])
    if len(values) != expected or not np.all(np.isfinite(values)) or np.any(values <= 0) or len(np.unique(values)) != expected:
        raise ValueError(f"Enter {expected} finite, positive, unique wavelengths in nm, in the selected-curve order.")
    return values


def result_zip(result, wavelengths=None, import_metadata=None, curve_indices=None) -> bytes:
    with tempfile.TemporaryDirectory(prefix="tcspc_export_") as folder:
        result.export(Path(folder))
        if wavelengths is not None:
            pd.DataFrame({"trace": np.arange(1, len(wavelengths) + 1), "wavelength_nm": wavelengths}).to_csv(Path(folder) / "wavelengths.csv", index=False)
        if import_metadata:
            (Path(folder) / "import_metadata.json").write_text(json.dumps(import_metadata, indent=2, default=_json_default, allow_nan=False), encoding="utf-8")
        if curve_indices is not None:
            pd.DataFrame({"trace": np.arange(1, len(curve_indices) + 1), "source_curve_index": curve_indices}).to_csv(Path(folder) / "curve_mapping.csv", index=False)
        buffer = BytesIO()
        with zipfile.ZipFile(buffer, "w", compression=zipfile.ZIP_DEFLATED) as archive:
            for item in sorted(Path(folder).iterdir()):
                if item.is_file():
                    archive.write(item, item.name)
        return buffer.getvalue()



def _sync_range_slider(time, slider_key, start_key, end_key):
    lo, hi = st.session_state[slider_key]
    st.session_state[start_key] = float(time[int(lo)])
    st.session_state[end_key] = float(time[int(hi)])


def _sync_range_numbers(time, slider_key, start_key, end_key):
    lo = int(np.argmin(np.abs(time - st.session_state[start_key])))
    hi = int(np.argmin(np.abs(time - st.session_state[end_key])))
    # Handles cannot cross; a collapsed interval stays visible and is not fitted.
    hi = max(lo, hi)
    st.session_state[slider_key] = (lo, hi)
    st.session_state[start_key] = float(time[lo])
    st.session_state[end_key] = float(time[hi])


def _continue_to_spectra():
    st.session_state["analysis_mode"] = "TRES / TRANES"
    st.session_state["spectra_workflow"] = "Latest global fit"


with st.sidebar:
    st.header("Data")
    st.caption("Import → select a model → fit → inspect residuals → export")
    mode = st.radio("Analysis", ["Single decay", "Global TRES", "TRES / TRANES"], key="analysis_mode")

if mode == "TRES / TRANES":
    render_spectra_workspace(st.session_state.get("last_global_fit"), plot_lock=_plot_lock())
    st.stop()

with st.sidebar:
    source_options = ["Upload ASCII data"]
    if mode == "Single decay":
        source_options.insert(0, "Bundled decay + IRF")
    source = st.selectbox("Source", source_options)
    time_unit = st.selectbox("Time unit in uploaded files", ["ns", "ps", "us", "ms", "s"], disabled=source == "Bundled decay + IRF", help="Choose your calibrated input unit; time is converted to ns. The bundled example uses an explicit ns assumption.")
    st.caption("Upload nonnegative signals on a uniform time grid. Use raw photon counts for Poisson fits. PlotBox exports use the time unit in their Time[...] header.")
    decay_upload = None
    irf_upload = None
    wavelength_text = ""
    if source == "Upload ASCII data":
        decay_upload = st.file_uploader("Decay or TRES ASCII file", type=["dat", "txt", "csv", "tsv"], key=f"data_{mode}")
        irf_upload = st.file_uploader("IRF ASCII file (optional for tail fit)", type=["dat", "txt", "csv", "tsv"], key=f"irf_{mode}")
        st.caption("Decay or IRF columns: time, signal. TRES columns: time, decay at λ1, decay at λ2, ….")
        if mode == "Global TRES":
            wavelength_text = st.text_input("Wavelengths in nm (optional for PlotBox)", placeholder="500, 520, 540", help="One value per selected decay curve, in selection order. PlotBox crv indices are not wavelengths. Plain TRES files need a numeric header or this mapping.")
    with st.expander("First time using the workbench?"):
        st.write("1. Choose the bundled decay and IRF, or upload your own files with the correct time unit.")
        st.write("2. Select the shaded fitting interval, then press Run fit.")
        st.write("3. Inspect the fitted curve and residuals, then download the result ZIP.")
        st.caption("Open USER_GUIDE.html in the extracted package for a step-by-step experimental workflow. README.html explains installation and startup.")
    st.markdown("Import numeric text files (`.csv`, `.dat`, `.txt`, or `.tsv`). Convert proprietary binary files to text before loading them.")

input_data = None
input_irf = None
import_metadata = {}
source_curve_indices = None
input_wavelengths = None
try:
    if source == "Bundled decay + IRF":
        # Bundled text files contain no time-unit metadata; ns is an explicit assumption.
        # Confirm the acquisition calibration before interpreting these lifetimes physically.
        input_data = load_decay(ROOT / "examples" / "Decay.dat", time_unit="ns")
        input_irf = load_decay(ROOT / "examples" / "IRF.dat", time_unit="ns")
        st.info("Loaded the supplied Decay.dat and IRF.dat examples, assuming the time values are in ns. Confirm the instrument calibration; the files contain no unit metadata.")
        st.caption("For this example, inspect an IRF background near 1.4 counts/bin and a pulse range of 3.8–8 ns. Choose the range from your own measured IRF when using other data.")
    elif decay_upload is not None:
        content = decay_upload.getvalue()
        if is_plotbox_content(content):
            bundle = load_plotbox(BytesIO(content))
            st.info(f"Multi-curve PlotBox export detected: {len(bundle.labels)} curves. Header unit {bundle.source_time_unit!r} is converted to ns; the plain-file unit selector is ignored.")
            for message in bundle.warnings:
                st.warning(message)
            with st.sidebar:
                st.subheader("PlotBox curve assignment")
                st.caption("Curve indices are zero-based. Select the measured decay(s) and IRF explicitly; curve names alone do not identify their roles.")
                options = list(range(len(bundle.labels)))
                label = lambda index: f"{index}: {bundle.labels[index]}"
                irf_index = st.selectbox("IRF curve in this file", [None, *options], format_func=lambda index: "No same-file IRF" if index is None else label(index), key=f"same_irf_{mode}")
                if mode == "Single decay":
                    decay_indices = [st.selectbox("Decay curve", options, format_func=label, key=f"decay_curve_{mode}")]
                else:
                    selection_method = st.radio("Choose decay curves by", ["Individual curves", "Curve index ranges"], key="curve_selection_method")
                    if selection_method == "Individual curves":
                        decay_indices = st.multiselect("Decay curves", options, default=[], format_func=label, key=f"decay_curves_{mode}_manual")
                    else:
                        curve_ranges = st.text_input("Decay curve indices/ranges", placeholder="1-51 or 1,3,5-8", help="Zero-based indices; ranges include both endpoints. Exclude the IRF and every curve that is not part of this ordinary decay series.")
                        decay_indices = parse_curve_indices(curve_ranges, limit=len(options)) if curve_ranges.strip() else []
            import_metadata = {"source_file": decay_upload.name, "source_sha256": hashlib.sha256(content).hexdigest(), "bundle": bundle.metadata}
            if not decay_indices:
                st.info("Select at least one decay curve in the sidebar. Exclude the IRF from the decay selection.")
            else:
                if irf_index in decay_indices:
                    raise ValueError("The IRF curve must be excluded from the selected decay curves.")
                selected = list(decay_indices)
                tt, cc = bundle.select_curves(selected)
                input_data = SimpleNamespace(time_ns=tt, counts=cc[:, 0] if mode == "Single decay" else cc)
                source_curve_indices = selected
                input_wavelengths = parse_wavelengths(wavelength_text, len(selected)) if mode == "Global TRES" else None
                import_metadata["decay_selection"] = bundle.selection_summary(selected)
                st.caption(f"Selected source curves: {selected}. Kept {len(tt):,} calibrated time bins; missing endpoint values are excluded, never replaced by zero.")
                if irf_index is not None:
                    if irf_upload is not None:
                        raise ValueError("Choose either a same-file IRF or a separate IRF upload.")
                    it, ic = bundle.select_curves([irf_index])
                    input_irf = SimpleNamespace(time_ns=it, counts=ic[:, 0])
                    import_metadata["irf_selection"] = bundle.selection_summary([irf_index])
            st.caption("Polarized VV/VH traces require a dedicated anisotropy model and G-factor. Select a magic-angle decay for ordinary lifetime fitting when available; do not interpret a positive-mixture global fit of polarized channels as anisotropy analysis.")
            with st.expander("PlotBox import audit", expanded=False):
                st.json(json.loads(json.dumps(import_metadata, default=_json_default)))
        elif mode == "Single decay":
            input_data = load_decay(BytesIO(content), time_unit=time_unit)
            import_metadata["decay_time_axis"] = input_data.metadata
        else:
            wavelengths_input = None
            if wavelength_text.strip():
                wavelengths_input = [float(item) for item in wavelength_text.replace(",", " ").split()]
            input_data = load_tres(BytesIO(content), time_unit=time_unit, wavelengths=wavelengths_input)
            import_metadata["decay_time_axis"] = input_data.metadata
            input_wavelengths = input_data.wavelengths
        if irf_upload is not None:
            irf_content = irf_upload.getvalue()
            if is_plotbox_content(irf_content):
                irf_bundle = load_plotbox(BytesIO(irf_content))
                with st.sidebar:
                    index = st.selectbox("IRF curve in separate upload", list(range(len(irf_bundle.labels))), format_func=lambda i: f"{i}: {irf_bundle.labels[i]}")
                it, ic = irf_bundle.select_curves([index])
                input_irf = SimpleNamespace(time_ns=it, counts=ic[:, 0])
                import_metadata["separate_irf"] = {"source_file": irf_upload.name, "source_sha256": hashlib.sha256(irf_content).hexdigest(), "bundle": irf_bundle.metadata, "selection": irf_bundle.selection_summary([index])}
                for message in irf_bundle.warnings:
                    st.warning(message)
            else:
                input_irf = load_decay(BytesIO(irf_content), time_unit=time_unit)
                import_metadata["irf_time_axis"] = input_irf.metadata
except (ValueError, OSError) as error:
    st.error(str(error))
    input_data = None

if input_data is None:
    st.info("Upload a decay or TRES matrix to begin. The bundled single-decay example is available in the sidebar.")
    st.stop()

t = input_data.time_ns
y = input_data.counts
data_identity = hashlib.sha256(decay_upload.getvalue()).hexdigest()[:12] if decay_upload is not None else "bundled"
range_identity = f"{data_identity}_{len(t)}_{t[0]:.14g}_{t[-1]:.14g}"
st.caption("Small rounding differences in exported time labels are handled automatically. Signal values stay in their original channels.")
st.caption(f"Loaded {len(t):,} time bins, {1 if y.ndim == 1 else y.shape[1]} decay(s), {t[0]:g}–{t[-1]:g} ns. Bin width: {np.median(np.diff(t)):g} ns.")

with st.expander("Inspect input", expanded=False):
    if y.ndim == 1:
        frame = pd.DataFrame({"time_ns": t, "counts": y})
    else:
        frame = pd.DataFrame(y, columns=[f"trace_{i}" for i in range(y.shape[1])])
        frame.insert(0, "time_ns", t)
    st.dataframe(frame.head(100), width="stretch", hide_index=True)

left, right = st.columns([1, 2])
with left:
    st.subheader("Fit settings")
    method = st.selectbox("Method", ["reconvolution", "tail"], format_func=lambda v: "IRF reconvolution" if v == "reconvolution" else "Exponential tail fit")
    if method == "tail":
        st.caption("Set the first fitted time after the IRF-affected rise; tail fits do not model that rise.")
    components = st.number_input("Exponential components", min_value=1, max_value=6, value=1, step=1)
    statistic = st.selectbox("Fit statistic", ["poisson", "neyman", "uniform"], format_func=lambda v: {"poisson": "Poisson likelihood (raw counts)", "neyman": "Observed-count weighted least squares", "uniform": "Unweighted least squares"}[v], help="Use Poisson likelihood for raw photon counts. Unweighted least squares is suitable for intensities with roughly constant noise variance.")
    default_start = float(t[0])
    if method == "tail":
        aggregate = y if y.ndim == 1 else y.sum(axis=1)
        peak_index = int(np.argmax(aggregate))
        default_start = float(t[min(len(t) - 2, peak_index + 3)])
    range_base = f"{mode}_{method}_{range_identity}"
    slider_key = f"fit_bins_{range_base}"
    start_key = f"fit_start_{range_base}"
    end_key = f"fit_end_{range_base}"
    if slider_key not in st.session_state:
        default_index = int(np.argmin(np.abs(t-default_start)))
        st.session_state[slider_key] = (default_index, len(t)-1)
        st.session_state[start_key] = float(t[default_index])
        st.session_state[end_key] = float(t[-1])
    range_args = (t, slider_key, start_key, end_key)
    fit_start = st.number_input("First fitted time (ns)", min_value=float(t[0]), max_value=float(t[-1]), format="%.6g", key=start_key, on_change=_sync_range_numbers, args=range_args)
    fit_end = st.number_input("Last fitted time (ns)", min_value=float(t[0]), max_value=float(t[-1]), format="%.6g", key=end_key, on_change=_sync_range_numbers, args=range_args)
    # The handles are displayed beside the kinetics preview. Their callback
    # updates these time fields before the next run, so button gating can use
    # the same selected channels without cropping the imported data.
    fit_bin_range = st.session_state[slider_key]
    selected_bins = int(fit_bin_range[1]-fit_bin_range[0]+1)
    st.caption("Choose the fitting interval with these time fields or the two handles above the kinetics preview.")
    if selected_bins < 4:
        st.warning("Select at least four time bins before fitting.")
    tau_text = st.text_input("Initial lifetimes (ns; optional)", placeholder="Example: 0.5, 3.0")
    fixed_tau_text = st.text_input("Fixed lifetimes (ns; optional)", placeholder="Leave blank to optimize lifetimes", help="If supplied, fixes every component lifetime and overrides initial lifetime guesses.")
    with st.expander("Background, shift, and excitation"):
        fit_background = st.checkbox("Fit a background for each decay", value=True)
        background = st.number_input("Decay background (input units)", min_value=0.0, value=0.0, format="%.6g", help="Starting or fixed background, depending on the checkbox.")
        irf_background = st.number_input("Subtract IRF background (IRF units)", min_value=0.0, value=0.0, format="%.6g", help="Subtract this constant, restrict to the selected pulse range if enabled, then clip negative IRF values to zero and normalize.")
        limit_irf = st.checkbox("Limit IRF pulse range", value=False, disabled=method == "tail", help="Exclude baseline outside the measured pulse. Clipping a noisy baseline over a long window can otherwise create a false broad IRF pedestal.")
        default_irf_start = 3.8 if source == "Bundled decay + IRF" else float(input_irf.time_ns[0]) if input_irf is not None else float(t[0])
        default_irf_end = 8.0 if source == "Bundled decay + IRF" else float(input_irf.time_ns[-1]) if input_irf is not None else float(t[-1])
        irf_start = st.number_input("First included IRF pulse time (ns)", value=default_irf_start, format="%.6g", disabled=not limit_irf or method == "tail", key=f"irf_start_{mode}_{data_identity}")
        irf_end = st.number_input("Last included IRF pulse time (ns)", value=default_irf_end, format="%.6g", disabled=not limit_irf or method == "tail", key=f"irf_end_{mode}_{data_identity}")
        fit_shift = st.checkbox("Fit one shared IRF shift", value=False, disabled=method == "tail")
        shift = st.number_input("IRF shift (ns)", value=0.0, format="%.6g", disabled=method == "tail", help="Positive values delay the IRF.")
        shift_low = st.number_input("Shift lower bound (ns)", value=-1.0, format="%.6g", disabled=method == "tail")
        shift_high = st.number_input("Shift upper bound (ns)", value=1.0, format="%.6g", disabled=method == "tail")
        periodic = st.checkbox("Use repetition-period correction", value=False, disabled=method == "tail")
        period = st.number_input("Laser repetition period (ns)", min_value=0.000001, value=max(float(len(t) * np.median(np.diff(t))), 1.0), format="%.6g", disabled=not periodic)
    with st.expander("Optimizer"):
        n_starts = st.number_input("Starting-point attempts", min_value=1, max_value=30, value=3, step=1)
        max_nfev = st.number_input("Evaluation limit per attempt", min_value=100, max_value=100000, value=1000, step=100)
    run = st.button("Run fit", type="primary", width="stretch", disabled=selected_bins < 4)

with right:
    if y.ndim > 1:
        with _plot_lock():
            fig, ax = plt.subplots(figsize=(9, 3.6))
            display_values = np.maximum(y.T, 0.5)
            im = ax.imshow(display_values, aspect="auto", origin="lower", extent=[t[0], t[-1], 0, y.shape[1]], norm=LogNorm(vmin=0.5, vmax=max(1.0, float(display_values.max()))), interpolation="nearest")
            ax.set(xlabel="Time (ns)", ylabel="Trace index", title="Measured time-by-wavelength signals")
            fig.colorbar(im, ax=ax, label="Signal (input units)")
            fig.tight_layout()
            st.pyplot(fig, width="stretch")
            plt.close(fig)
        st.caption("Wavelength labels are preserved when supplied. Otherwise amplitudes are shown against source curve indices; no wavelength mapping is assumed.")

    st.subheader("Kinetics preview and fit interval")
    fit_bin_range = st.slider("Fit interval handles (time bins)", min_value=0, max_value=len(t)-1, step=1, key=slider_key, on_change=_sync_range_slider, args=range_args,
                              help="Drag the two handles to select the first and last fitted channels. The time fields and plot update together. Signal values stay in their original channels.")
    st.caption(f"Selected X-axis interval: {fit_start:.6g}–{fit_end:.6g} ns · {selected_bins:,} channels. The green region marks the channels that will be fitted.")
    st.caption("Drag the left and right handles to choose the start and end, or type times in the boxes on the left. Switch to Selected fit window to inspect that part of the curve more closely.")

    if y.ndim > 1:
        def preview_label(index):
            if index == "sum":
                return "Sum of selected traces"
            if input_wavelengths is not None:
                return f"{input_wavelengths[index]:g} nm"
            if source_curve_indices is not None:
                return f"Source curve {source_curve_indices[index]}"
            return f"Trace {index+1}"
        preview_trace = st.selectbox("Kinetics preview trace", ["sum", *range(y.shape[1])], format_func=preview_label, key=f"preview_trace_{range_base}")
        preview_counts = y.sum(axis=1) if preview_trace == "sum" else y[:, preview_trace]
        preview_name = preview_label(preview_trace)
    else:
        preview_counts = y
        preview_name = "Measured decay"
    preview_view = st.radio("Kinetics preview view", ["Whole record", "Selected fit window"], horizontal=True, key=f"preview_view_{range_base}")
    with _plot_lock():
        fig, ax = plt.subplots(figsize=(9, 3.9))
        ax.plot(t, np.maximum(preview_counts, 0.5), ".", ms=2, color="#365d80", label=preview_name)
        if input_irf is not None:
            scale = float(np.max(preview_counts))/max(float(np.max(input_irf.counts)), 1.0)
            ax.plot(input_irf.time_ns, np.maximum(input_irf.counts*scale, 0.5), lw=0.8, alpha=0.7, label="IRF (scaled for display)")
            if limit_irf and method == "reconvolution":
                ax.axvspan(irf_start, irf_end, color="#efba8b", alpha=0.25, label="IRF pulse range")
        ax.axvspan(fit_start, fit_end, color="#89bea4", alpha=0.28, label="Included in fit")
        ax.axvline(fit_start, color="#1b7751", lw=1.1)
        ax.axvline(fit_end, color="#1b7751", lw=1.1)
        ax.set(yscale="log", xlabel="Time (ns)", ylabel="Signal (input units)",
               title="Fit interval on measured decay" if y.ndim == 1 else "Global-fit interval on measured kinetics")
        if preview_view == "Selected fit window" and fit_end > fit_start:
            ax.set_xlim(fit_start, fit_end)
        ax.legend(fontsize="small")
        fig.tight_layout()
        st.pyplot(fig, width="stretch")
        plt.close(fig)
    if y.ndim > 1:
        st.caption("The green region is fitted for every selected decay. Changing the preview trace or its zoom does not change which decays are fitted.")
    else:
        st.caption("The green region is fitted when you press Run fit. Changing the zoom only changes the view; the handles and time fields select the fitted channels.")
    st.caption("Reconvolution uses the full imported time grid before applying the fit interval. Values below 0.5 appear at 0.5 on logarithmic plots. The fit uses the full precision of the signal.")
    st.write("Global fits share component lifetimes and fit separate nonnegative amplitudes and backgrounds for each selected trace." if mode == "Global TRES" else "Inspect the residuals as well as the fitted lifetime. A smooth-looking decay can still need a different model or fit range.")

if run:
    try:
        if method == "reconvolution" and input_irf is None:
            raise ValueError("Upload an IRF for reconvolution, or choose a tail fit.")
        if fit_start >= fit_end:
            raise ValueError("The first fitted time must be before the last fitted time.")
        fixed_lifetimes = parse_lifetimes(fixed_tau_text, int(components))
        initial_lifetimes = None if fixed_lifetimes else parse_lifetimes(tau_text, int(components))
        if shift_low >= shift_high:
            raise ValueError("The shift lower bound must be less than the upper bound.")
        if limit_irf and method == "reconvolution" and irf_start >= irf_end:
            raise ValueError("The IRF pulse range must begin before it ends.")
        config = FitConfig(
            n_components=int(components), method=method, statistic=statistic,
            fit_start=float(fit_start), fit_end=float(fit_end),
            initial_lifetimes=None if fixed_lifetimes else initial_lifetimes,
            fixed_lifetimes=fixed_lifetimes, fit_background=fit_background,
            background=float(background), irf_background=float(irf_background),
            irf_start=float(irf_start) if limit_irf and method == "reconvolution" else None,
            irf_end=float(irf_end) if limit_irf and method == "reconvolution" else None,
            fit_shift=bool(fit_shift and method == "reconvolution"),
            shift=float(shift) if method == "reconvolution" else 0.0,
            shift_bounds=(float(shift_low), float(shift_high)),
            repetition_period=float(period) if periodic and method == "reconvolution" else None,
            n_starts=int(n_starts), max_nfev=int(max_nfev),
        )
        kwargs = {"config": config}
        if input_irf is not None and method == "reconvolution":
            kwargs.update(irf=input_irf.counts, irf_time=input_irf.time_ns)
        with st.spinner("Fitting decays…"):
            fitter = fit_decay if y.ndim == 1 else fit_global
            result = fitter(t, y, **kwargs)
        wavelengths = None if y.ndim == 1 else input_wavelengths
        st.session_state["last_fit"] = {
            "result": result, "wavelengths": wavelengths,
            "source": "Supplied example" if source == "Bundled decay + IRF" else decay_upload.name,
            "method": method, "statistic": statistic,
            "import_metadata": import_metadata, "curve_indices": source_curve_indices,
            "data_id": data_identity,
            "analysis_mode": mode,
        }
        if mode == "Global TRES":
            st.session_state["last_global_fit"] = st.session_state["last_fit"]
    except (ValueError, RuntimeError, OSError) as error:
        st.error(str(error))

payload = st.session_state.get("last_fit")
if payload is not None:
    result = payload["result"]
    wavelengths = payload["wavelengths"]
    curve_indices = payload.get("curve_indices")
    st.divider()
    st.subheader("Most recent fit")
    st.caption(f"Source: {payload['source']} · {payload['method']} · {payload['statistic']}. Changed settings take effect when you run another fit.")
    used_bins = np.asarray(result.fit_mask, dtype=bool)
    st.caption(f"Last fit used {result.time[used_bins][0]:.6g}–{result.time[used_bins][-1]:.6g} ns ({int(used_bins.sum()):,} channels per trace). Current preview changes are applied only when Run fit is pressed.")
    diagnostics = result.to_dict().get("diagnostics", {})
    metric_columns = st.columns(3)
    metric_columns[0].metric("Convergence", "Converged" if diagnostics.get("converged") else "Review fit")
    pearson_value = diagnostics.get("reduced_pearson_chi_square")
    deviance_value = diagnostics.get("reduced_poisson_deviance")
    metric_columns[1].metric("Reduced Pearson χ²", f"{pearson_value:.4g}" if pearson_value is not None else "Unavailable")
    metric_columns[2].metric("Reduced Poisson deviance", f"{deviance_value:.4g}" if deviance_value is not None else "Unavailable")
    st.caption("Pearson χ² and Poisson deviance are count-based diagnostics. For analogue intensities, assess the residuals in signal units using the unweighted least-squares objective when the noise variance is roughly constant.")
    for warning in result.warnings:
        st.warning(warning)
    stderr = result.lifetime_stderr
    table = {"Component": np.arange(1, len(result.lifetimes) + 1), "Lifetime (ns)": result.lifetimes}
    if stderr is not None:
        table["Approximate standard error (ns)"] = stderr
    st.dataframe(pd.DataFrame(table), width="stretch", hide_index=True)
    st.caption("Standard errors are local approximations conditional on the measured IRF and selected model. They do not include instrument calibration uncertainty.")
    trace_count = result.observed.shape[1]
    trace_index = 0
    if trace_count > 1:
        trace_index = st.slider("Inspect fitted trace", min_value=0, max_value=trace_count - 1, value=0)
        if wavelengths is not None:
            st.caption(f"Selected wavelength (nm): {wavelengths[trace_index]}")
        elif curve_indices is not None:
            st.caption(f"Selected source curve index: {curve_indices[trace_index]}")
    fit_mask = np.asarray(result.fit_mask, dtype=bool)
    with _plot_lock():
        fig, (ax, rx) = plt.subplots(2, 1, figsize=(11, 5.4), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
        ax.plot(result.time, np.maximum(result.observed[:, trace_index], 0.5), ".", ms=2, label="Measured")
        ax.plot(result.time[fit_mask], np.maximum(result.fitted[fit_mask, trace_index], 0.5), lw=1.4, label="Fit")
        ax.set(yscale="log", ylabel="Signal (input units)")
        ax.legend()
        rx.axhline(0, color="0.5", lw=0.8)
        rx.plot(result.time, result.residuals[:, trace_index], lw=0.6)
        rx.set(xlabel="Time (ns)", ylabel="Fit residual")
        fig.tight_layout()
        st.pyplot(fig, width="stretch")
        plt.close(fig)
    st.caption("Residuals use the selected fit statistic and cover the fitted interval. The logarithmic display places values below 0.5 at 0.5.")
    autocorrelations = diagnostics.get("residual_autocorrelation", [])
    if trace_index < len(autocorrelations):
        acf = np.asarray(autocorrelations[trace_index], dtype=float)
        with _plot_lock():
            fig, ax = plt.subplots(figsize=(11, 2.3))
            ax.axhline(0, color="0.5", lw=0.8)
            ax.plot(np.arange(len(acf)), acf, "o-", ms=3, lw=0.8)
            ax.set(xlabel="Lag (time bins)", ylabel="Autocorrelation", title="Residual autocorrelation")
            fig.tight_layout()
            st.pyplot(fig, width="stretch")
            plt.close(fig)
    if trace_count > 1:
        amplitude_frame = pd.DataFrame(result.amplitudes, columns=[f"A{i+1}" for i in range(result.amplitudes.shape[1])])
        amplitude_coordinates = wavelengths if wavelengths is not None else curve_indices if curve_indices is not None else np.arange(trace_count)
        amplitude_frame.insert(0, "wavelength_nm" if wavelengths is not None else "source_curve_index", amplitude_coordinates)
        st.write("Pre-exponential amplitude spectra (not reconstructed TRES)")
        st.dataframe(amplitude_frame, width="stretch", hide_index=True)
        with _plot_lock():
            fig, ax = plt.subplots(figsize=(11, 3))
            coordinates = amplitude_coordinates
            for component, lifetime in enumerate(result.lifetimes):
                ax.plot(coordinates, result.amplitudes[:, component], "o-", ms=3, label=f"τ = {lifetime:.5g} ns")
            ax.set(xlabel="Wavelength (nm)" if wavelengths is not None else "Source curve index", ylabel="Pre-exponential amplitude")
            ax.legend()
            fig.tight_layout()
            st.pyplot(fig, width="stretch")
            plt.close(fig)
        st.caption("These are fitted pre-exponential amplitudes, not steady-state-scaled TRES. Reconstruction below requires a measured steady-state spectrum. Amplitudes are constrained nonnegative; signed rise components and kinetic target analysis are not implemented.")
    with st.expander("Fit diagnostics"):
        st.json(result.to_dict().get("diagnostics", {}))
    st.download_button("Download fit results (JSON + CSV ZIP)", data=result_zip(result, wavelengths, payload.get("import_metadata"), curve_indices), file_name="tcspc_fit_results.zip", mime="application/zip", width="content")

    if payload.get("analysis_mode") == "Global TRES":
        st.subheader("Next: TRES and TRANES")
        st.write("Continue with these wavelengths, lifetimes and pre-exponential amplitudes. Only the matching steady-state fluorescence intensities are needed next.")
        st.button("Continue to TRES / TRANES", on_click=_continue_to_spectra,
                  width="content", key="continue_to_spectra")
