"""Steady-state-constrained emission reconstruction from fitted exponentials.

This module evaluates intrinsic (IRF-free) model decays. Decay-associated
amplitudes alone are not time-resolved emission spectra. All times use ns and
all wavelength coordinates use nm; these units must be supplied correctly.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import re

import numpy as np
from scipy.integrate import trapezoid

from .io import DataFormatError, _read_source, _tokens, _numeric_row


@dataclass(frozen=True)
class FittedSpectra:
    wavelengths: np.ndarray
    amplitudes: np.ndarray
    steady_state: np.ndarray | None
    component_names: tuple[str, ...]
    metadata: dict


@dataclass(frozen=True)
class SteadyStateSpectrum:
    wavelengths: np.ndarray
    intensities: np.ndarray
    metadata: dict


@dataclass(frozen=True)
class SpectralResult:
    wavelengths: np.ndarray
    times: np.ndarray
    lifetimes: np.ndarray
    amplitudes: np.ndarray
    steady_state: np.ndarray
    tres: np.ndarray
    tranes: np.ndarray
    metadata: dict
    warnings: list[str]

    def export(self, directory):
        """Save matrices with wavelength rows and one column per requested delay."""
        path = Path(directory)
        path.mkdir(parents=True, exist_ok=True)
        header = "wavelength_nm," + ",".join(f"{t:.17g}_ns" for t in self.times)
        for name, values in [("tres.csv", self.tres), ("tranes.csv", self.tranes)]:
            np.savetxt(path / name, np.column_stack([self.wavelengths, values]),
                       delimiter=",", header=header, comments="", fmt="%.15g")
        np.savetxt(path / "sampled_steady_state.csv",
                   np.column_stack([self.wavelengths, self.steady_state]),
                   delimiter=",", header="wavelength_nm,steady_state_intensity",
                   comments="", fmt="%.15g")
        if "steady_state_source_normalized_intensities" in self.metadata:
            np.savetxt(path / "normalized_steady_state.csv",
                       np.column_stack([
                           self.metadata["steady_state_source_wavelengths_nm"],
                           self.metadata["steady_state_source_normalized_intensities"],
                       ]), delimiter=",",
                       header="wavelength_nm,normalized_intensity",
                       comments="", fmt="%.15g")
        np.savetxt(path / "fitted_amplitudes.csv",
                   np.column_stack([self.wavelengths, self.amplitudes]),
                   delimiter=",", header="wavelength_nm," + ",".join(
                       f"A{i+1}" for i in range(len(self.lifetimes))),
                   comments="", fmt="%.15g")
        (path / "spectral_metadata.json").write_text(
            json.dumps({**self.metadata, "warnings": self.warnings},
                       indent=2, allow_nan=False) + "\n", encoding="utf-8")
        return path


def _table(source, delimiter=None):
    """Read small spectral tables; two rows suffice (unlike a TCSPC histogram)."""
    text, name = _read_source(source)
    headers, rows = [], []
    width = None
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(("#", "%", "//", ";")):
            continue
        tokens = _tokens(line, delimiter)
        if any(token.lower() in {"nan", "+nan", "-nan", "inf", "+inf", "-inf", "infinity", "+infinity", "-infinity"} for token in tokens):
            raise DataFormatError(f"Line {number}: spectral values must be finite.")
        values = _numeric_row(tokens)
        if values is None:
            if tokens and _numeric_row(tokens[:1]) is not None:
                raise DataFormatError(f"Line {number}: a wavelength row contains missing or nonnumeric values.")
            if rows:
                raise DataFormatError(f"Line {number}: expected numeric spectral data; missing cells are unsupported.")
            headers.append(tokens)
            continue
        if width is None:
            width = len(values)
        if len(values) != width:
            raise DataFormatError(f"Line {number}: spectral column count changed.")
        rows.append(values)
    if len(rows) < 2 or width is None:
        raise DataFormatError("A spectrum needs at least two numeric wavelength rows.")
    return np.asarray(rows, dtype=float), headers, name


def _wavelength_order(wavelengths):
    waves = np.asarray(wavelengths, dtype=float)
    if waves.ndim != 1 or len(waves) < 2 or not np.all(np.isfinite(waves)) or np.any(waves <= 0):
        raise ValueError("Supply at least two finite positive wavelengths in nm.")
    if len(np.unique(waves)) != len(waves):
        raise ValueError("Wavelengths must be unique; duplicate spectral rows need resolving before reconstruction.")
    order = np.argsort(waves)
    return waves[order].copy(), order


def load_fitted_spectra(source, *, has_steady_state=None, delimiter=None):
    """Read ``wavelength, A1, ..., AK[, steady_state]`` from text.

    A named last column containing ``steadystate`` or ``ss`` is recognized.
    Without a header, all columns after wavelength are amplitudes unless
    ``has_steady_state=True`` is explicit. Lifetime order is supplied separately
    and must match the amplitude-column order.
    """
    table, headers, name = _table(source, delimiter)
    if table.shape[1] < 2:
        raise DataFormatError("Fitted input needs wavelength and amplitude columns.")
    header = next((h for h in reversed(headers) if len(h) == table.shape[1]), None)
    if has_steady_state is None:
        last = re.sub(r"[^a-z0-9]", "", header[-1].lower()) if header else ""
        has_steady_state = "steadystate" in last or last in {"ss", "iss", "ssintensity"}
    if has_steady_state and table.shape[1] < 3:
        raise DataFormatError("Combined fitted input needs wavelength, at least one amplitude, and steady-state intensity.")
    stop = -1 if has_steady_state else None
    amplitudes = table[:, 1:stop]
    if not np.all(np.isfinite(amplitudes)) or np.any(amplitudes < 0):
        raise DataFormatError("This reconstruction requires finite nonnegative parallel-model amplitudes.")
    waves, order = _wavelength_order(table[:, 0])
    steady = table[:, -1].copy() if has_steady_state else None
    if steady is not None and (np.any(steady < 0) or not np.all(np.isfinite(steady))):
        raise DataFormatError("Steady-state intensities must be finite and nonnegative.")
    component_names = tuple(header[1:stop]) if header else tuple(f"A{i+1}" for i in range(amplitudes.shape[1]))
    return FittedSpectra(waves, amplitudes[order].copy(),
                         None if steady is None else steady[order], component_names,
                         {"source_name": name, "n_components": amplitudes.shape[1],
                          "has_steady_state": bool(has_steady_state),
                          "component_order": list(component_names), "wavelength_unit": "nm"})


load_fitted_table = load_fitted_spectra


def load_steady_state(source, *, delimiter=None, allow_negative=False):
    """Read exactly two columns: wavelength in nm and steady-state intensity.

    ``allow_negative=True`` permits raw signed intensities for an explicitly
    selected min-max normalization step. Reconstruction still rejects negative
    values when another normalization is selected.
    """
    table, _, name = _table(source, delimiter)
    if table.shape[1] != 2:
        raise DataFormatError("Steady-state input needs exactly two columns: wavelength_nm, intensity.")
    waves, order = _wavelength_order(table[:, 0])
    intensity = table[:, 1][order].copy()
    if not np.all(np.isfinite(intensity)):
        raise DataFormatError("Steady-state intensities must be finite.")
    if not allow_negative and (np.any(intensity < 0) or not np.any(intensity > 0)):
        raise DataFormatError("Steady-state intensity must be finite, nonnegative, and positive somewhere.")
    return SteadyStateSpectrum(waves, intensity,
                              {"source_name": name, "wavelength_unit": "nm",
                               "intensity_unit": "supplied relative units",
                               "negative_raw_values_allowed": bool(allow_negative)})


def reconstruct_spectra(wavelengths, amplitudes, lifetimes, times, steady_state=None,
                        *, steady_wavelengths=None, steady_intensities=None,
                        steady_normalization="none", amplitude_reference_delay_ns=0,
                        steady_alignment="interpolate"):
    """Reconstruct TRES and wavelength-area-normalized TRANES.

    ``I(lam,t) = SS(lam) * sum(A_k*exp(-t/tau_k))/sum(A_k*tau_k)``.
    Amplitudes are columns, wavelengths rows, and lifetimes retain column order.
    With ``steady_alignment='interpolate'``, separate/dense steady-state input is
    linearly interpolated within coverage; extrapolation is rejected. ``'exact'``
    requires one unique matching source row per fitted wavelength, within an
    absolute 1e-6 nm tolerance, and does not interpolate. Input normalization
    precedes sampling. Min-max normalization subtracts the source minimum and
    divides by its range; signed raw intensities are accepted only in that mode.
    The reference delay is zero for excitation-origin amplitudes. Tail amplitudes
    need the elapsed time from excitation to their reference (in ns).
    """
    waves, order = _wavelength_order(wavelengths)
    aa = np.asarray(amplitudes, dtype=float)
    tau = np.asarray(lifetimes, dtype=float)
    tt = np.asarray(times, dtype=float)
    if aa.ndim == 1:
        aa = aa[:, None]
    if aa.ndim != 2 or aa.shape[0] != len(waves) or aa.shape[1] < 1:
        raise ValueError("Amplitudes must have shape [wavelength, component].")
    if not np.all(np.isfinite(aa)) or np.any(aa < 0):
        raise ValueError("Reconstruction requires finite nonnegative amplitudes.")
    if tau.ndim != 1 or len(tau) != aa.shape[1] or not np.all(np.isfinite(tau)) or np.any(tau <= 0):
        raise ValueError("Supply one finite positive lifetime in ns for each amplitude column, in the same order.")
    if tt.ndim != 1 or len(tt) < 1 or not np.all(np.isfinite(tt)) or np.any(tt < 0):
        raise ValueError("Reconstruction delays must be finite, nonnegative times in ns.")
    if len(np.unique(tt)) != len(tt):
        raise ValueError("Reconstruction delays must be unique.")
    if steady_normalization not in {"none", "peak", "area", "minmax"}:
        raise ValueError("Steady-state normalization must be 'none', 'peak', 'area', or 'minmax'.")
    if steady_alignment not in {"interpolate", "exact"}:
        raise ValueError("Steady-state alignment must be 'interpolate' or 'exact'.")
    delay = float(amplitude_reference_delay_ns)
    if not np.isfinite(delay) or delay < 0:
        raise ValueError("Amplitude-reference delay after excitation must be finite and nonnegative.")
    aa = aa[order].copy()
    exact_indices = None
    matching_tolerance = 1e-6
    if isinstance(steady_state, SteadyStateSpectrum):
        if steady_wavelengths is not None or steady_intensities is not None:
            raise ValueError("Supply only one steady-state spectrum source.")
        steady_wavelengths = steady_state.wavelengths
        steady_intensities = steady_state.intensities
        steady_state = None
    if steady_state is not None:
        if steady_wavelengths is not None or steady_intensities is not None:
            raise ValueError("Supply aligned steady-state values or a separate spectrum, not both.")
        source_waves = waves.copy()
        source_ss = np.asarray(steady_state, dtype=float)
        if source_ss.shape != (len(waves),):
            raise ValueError("Supply one steady-state value per fitted wavelength.")
        source_ss = source_ss[order].copy()
        interpolated = False
    else:
        if steady_wavelengths is None or steady_intensities is None:
            raise ValueError("Steady-state intensity is required. Amplitude spectra alone do not provide TRES.")
        source_waves, ss_order = _wavelength_order(steady_wavelengths)
        source_ss = np.asarray(steady_intensities, dtype=float)
        if source_ss.shape != (len(source_waves),):
            raise ValueError("Steady-state wavelengths and intensities have different lengths.")
        source_ss = source_ss[ss_order].copy()
        if steady_alignment == "exact":
            if len(source_waves) != len(waves):
                raise ValueError("Exact steady-state matching requires exactly one row per fitted wavelength.")
            left = np.searchsorted(source_waves, waves-matching_tolerance, side="left")
            right = np.searchsorted(source_waves, waves+matching_tolerance, side="right")
            if np.any(right-left != 1) or len(np.unique(left)) != len(waves):
                raise ValueError("Exact steady-state matching requires one unambiguous row at every fitted wavelength (within 1e-6 nm).")
            exact_indices = left
            interpolated = False
        else:
            if waves[0] < source_waves[0] or waves[-1] > source_waves[-1]:
                raise ValueError("Steady-state spectrum must cover every fitted wavelength; extrapolation is not performed.")
            interpolated = not np.array_equal(waves, source_waves)
    if not np.all(np.isfinite(source_ss)):
        raise ValueError("Steady-state intensities must be finite.")
    if steady_normalization != "minmax" and (np.any(source_ss < 0) or not np.any(source_ss > 0)):
        raise ValueError("Steady-state intensities must be finite, nonnegative, and positive somewhere; signed raw values require min-max normalization.")
    offset = float(source_ss.min()) if steady_normalization == "minmax" else 0.0
    normalizer = (float(source_ss.max()) if steady_normalization == "peak" else
                  float(trapezoid(source_ss, source_waves)) if steady_normalization == "area" else
                  float(source_ss.max()-offset) if steady_normalization == "minmax" else 1.0)
    if not np.isfinite(normalizer) or normalizer <= 0:
        raise ValueError("The steady-state normalization denominator must be positive and finite.")
    normalized_source = (source_ss-offset) / normalizer
    if exact_indices is not None:
        sampled = normalized_source[exact_indices].copy()
    elif steady_state is not None or not interpolated:
        sampled = normalized_source.copy()
    else:
        sampled = np.interp(waves, source_waves, normalized_source)
    if not np.any(sampled > 0):
        raise ValueError("The supplied steady-state spectrum is zero at every fitted wavelength.")

    # Per-wavelength common amplitude scale cancels from the formula. Log-space
    # weights also permit a tail-reference conversion without overflowing A0.
    log_weights = np.full_like(aa, -np.inf)
    positive = aa > 0
    with np.errstate(over="ignore", divide="ignore"):
        reference_correction = delay / tau
    if not np.all(np.isfinite(reference_correction)):
        raise ValueError("Amplitude-reference correction exceeds numerical range; check lifetimes and delay.")
    log_weights[positive] = np.log(aa[positive])
    log_weights += reference_correction
    row_max = np.max(log_weights, axis=1)
    supported = np.isfinite(row_max)
    if np.any((~supported) & (sampled > 0)):
        raise ValueError("A positive steady-state intensity has no supported fitted decay amplitude at that wavelength.")
    weights = np.zeros_like(aa)
    weights[supported] = np.exp(log_weights[supported] - row_max[supported, None])
    area_fit = weights @ tau
    model = np.exp(-tt[None, :] / tau[:, None])
    tres = np.zeros((len(waves), len(tt)))
    tres[supported] = sampled[supported, None] * ((weights[supported] @ model) / area_fit[supported, None])
    spectral_areas = trapezoid(tres, waves, axis=0)
    if np.any(~np.isfinite(spectral_areas)) or np.any(spectral_areas <= 0):
        raise ValueError("Requested delays produce zero or non-finite spectral area. Use earlier times or check inputs.")
    tranes = tres / spectral_areas[None, :]
    warnings = ["Reconstructed intensities are model-based relative spectra, not an absolute radiometric calibration.",
                "Reconstruction outside the measured post-excitation delay range is model extrapolation."]
    if delay:
        warnings.append("Tail amplitudes were extrapolated back to excitation using the supplied reference delay; early-time reconstruction depends on that model extrapolation.")
    metadata = {
        "software": "TCSPC Workbench 0.1.0 spectral reconstruction",
        "formula": "I(lambda,t)=SS(lambda)*sum(Ak*exp(-t/tauk))/sum(Ak*tauk)",
        "tranes_definition": "TRES divided by its trapezoidal area over the fitted wavelength grid, separately at each delay",
        "wavelength_unit": "nm", "time_unit": "ns",
        "tres_unit": {"none": "supplied steady-state intensity units / ns",
                      "peak": "peak-normalized steady-state intensity / ns",
                      "minmax": "min-max-normalized steady-state intensity / ns",
                      "area": "1/(nm ns) over the supplied steady-state wavelength interval"}[steady_normalization],
        "tranes_unit": "1/nm over the selected wavelength interval",
        "lifetimes_ns_in_amplitude_column_order": tau.tolist(),
        "reconstruction_delays_ns": tt.tolist(),
        "wavelengths_nm": waves.tolist(),
        "steady_state_normalization": steady_normalization,
        "steady_state_normalization_offset": offset,
        "steady_state_normalization_divisor": normalizer,
        "steady_state_alignment": steady_alignment,
        "steady_state_matching_tolerance_nm": matching_tolerance if steady_alignment == "exact" else None,
        "steady_state_source_wavelengths_nm": source_waves.tolist(),
        "steady_state_source_intensities": source_ss.tolist(),
        "steady_state_source_normalized_intensities": normalized_source.tolist(),
        "sampled_steady_state_intensities": sampled.tolist(),
        "steady_state_interpolation": "none; one-to-one wavelength match" if steady_alignment == "exact" else "linear within measured coverage" if interpolated else "provided at fitted wavelengths",
        "spectral_extrapolation": False,
        "spectral_areas_before_tranes_normalization": spectral_areas.tolist(),
        "amplitude_reference_delay_ns": delay,
        "amplitude_columns": [f"A{i+1}" for i in range(len(tau))],
        "spectral_smoothing_applied": False,
        "steady_state_is_absolute_calibration": False,
        "spectral_interval_nm": [float(waves[0]), float(waves[-1])],
    }
    return SpectralResult(waves, tt.copy(), tau.copy(), aa, sampled, tres, tranes,
                          metadata, warnings)
