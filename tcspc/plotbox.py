"""Read PlotBox ASCII exports without guessing curve identities.

PlotBox uses ``-`` for absent cells. Missing counts remain NaN until the caller
selects curves; selection can trim incomplete endpoints but never fills counts
or removes internal bins. Header-declared time units are converted to ns.
"""
from __future__ import annotations

from dataclasses import dataclass
import csv
import re

import numpy as np

from .io import DataFormatError, _read_source


_UNIT_SCALE = {"ps": 1e-3, "ns": 1.0, "us": 1e3, "ms": 1e6, "s": 1e9}
_TIME_HEADER = re.compile(r"^Time\s*\[\s*(ps|ns|us|µs|μs|ms|s)\s*\]\s*X$", re.I)


def _unit(value):
    normalized = str(value).strip().lower().replace("µ", "u").replace("μ", "u")
    if normalized not in _UNIT_SCALE:
        raise DataFormatError("Choose ps, ns, us, ms, or s as a time unit.")
    return normalized


@dataclass(frozen=True)
class CurveBundle:
    """Unassigned raw curves sharing a calibrated, canonical time grid.

    ``counts`` may contain NaN for absent cells; only ``select_curves`` produces
    engine-ready complete data. Raw arrays include rows without time stamps.
    """

    time_ns: np.ndarray
    counts: np.ndarray
    labels: tuple[str, ...]
    source_time_unit: str
    name: str
    warnings: list[str]
    metadata: dict
    original_time_ns: np.ndarray
    raw_time: np.ndarray
    raw_counts: np.ndarray

    def _selection(self, indices, missing_policy):
        if missing_policy not in ("intersection", "raise"):
            raise DataFormatError("missing_policy must be 'intersection' or 'raise'.")
        try:
            selected = list(indices)
        except TypeError as exc:
            raise DataFormatError("Select a sequence of zero-based curve indices.") from exc
        if not selected or any(isinstance(i, (bool, np.bool_)) or not isinstance(i, (int, np.integer)) for i in selected):
            raise DataFormatError("Select at least one integer curve index.")
        if len(set(selected)) != len(selected) or min(selected) < 0 or max(selected) >= self.counts.shape[1]:
            raise DataFormatError("Curve indices must be unique and within the bundle's curve range.")
        counts = self.counts[:, selected]
        complete = np.all(np.isfinite(counts), axis=1)
        if missing_policy == "raise" and not np.all(complete):
            raise DataFormatError("Selected curves contain missing counts ('-'); choose endpoint intersection or another curve.")
        valid = np.flatnonzero(complete)
        if len(valid) < 4:
            raise DataFormatError("Selected curves need at least four complete time bins.")
        first, last = int(valid[0]), int(valid[-1])
        if not np.all(complete[first:last+1]):
            raise DataFormatError("Missing counts occur inside the selected curves' time range. Internal bins cannot be dropped or filled for reconvolution.")
        return selected, counts, first, last

    def select_curves(self, indices, *, missing_policy="intersection"):
        """Return ``(time_ns, counts[time, curve])`` with missing endpoints trimmed.

        Curve indices are zero-based positions, preserving their requested order.
        Internal missing bins are rejected. No missing value becomes a zero.
        See ``selection_summary`` for an audit of omitted endpoint observations.
        """
        _, counts, first, last = self._selection(indices, missing_policy)
        return self.time_ns[first:last+1].copy(), counts[first:last+1].copy()

    def selection_summary(self, indices, *, missing_policy="intersection"):
        selected, counts, first, last = self._selection(indices, missing_policy)
        omit = np.ones(len(self.time_ns), dtype=bool)
        omit[first:last+1] = False
        omitted_counts = counts[omit]
        return {
            "curve_indices": [int(i) for i in selected],
            "curve_labels": [self.labels[i] for i in selected],
            "missing_policy": missing_policy,
            "retained_time_bins": last-first+1,
            "leading_endpoint_rows_trimmed": first,
            "trailing_endpoint_rows_trimmed": len(self.time_ns)-last-1,
            "finite_count_cells_omitted_at_endpoints": int(np.isfinite(omitted_counts).sum()),
            "counts_omitted_at_endpoints": float(np.nansum(omitted_counts)),
            "fit_start_available_ns": float(self.time_ns[first]),
            "fit_end_available_ns": float(self.time_ns[last]),
            "global_import_dropped_timestamp_rows": self.metadata["dropped_timestamp_data_rows"],
        }


def is_plotbox_export(source):
    """Return whether a text input begins with the PlotBox ASCII marker."""
    text, _ = _read_source(source)
    return text.lstrip().startswith("PlotBox Export")


def load_plotbox(source, *, time_unit=None, grid_tolerance_bins=1e-4):
    """Read a PlotBox ASCII file into a :class:`CurveBundle`.

    The explicit ``Time[s]X`` header establishes seconds for the supplied files.
    A conflicting optional ``time_unit`` is rejected rather than silently
    reinterpreting the data. Small coordinate quantization is canonicalized to
    a uniform grid anchored at the original first timestamp, with its maximum
    displacement recorded. The default tolerance is 0.0001 of a bin.
    """
    text, name = _read_source(source)
    lines = [(i+1, line.strip()) for i, line in enumerate(text.splitlines()) if line.strip()]
    if len(lines) < 3 or lines[0][1] != "PlotBox Export":
        raise DataFormatError("Expected a 'PlotBox Export' ASCII file.")
    header = [cell.strip() for cell in next(csv.reader([lines[1][1]], delimiter="\t"))]
    if len(header) < 2:
        raise DataFormatError("PlotBox export needs one time column and at least one curve column.")
    match = _TIME_HEADER.fullmatch(header[0])
    if match is None:
        raise DataFormatError("Expected an explicit time header such as 'Time[s]X'.")
    source_unit = _unit(match.group(1))
    if time_unit is not None and _unit(time_unit) != source_unit:
        raise DataFormatError(f"Selected time unit {time_unit!r} conflicts with header unit {source_unit!r}.")
    labels = tuple(header[1:])
    if any(not label for label in labels) or len(set(labels)) != len(labels):
        raise DataFormatError("Curve column labels must be nonempty and unique.")
    rows = []
    for number, line in lines[2:]:
        fields = next(csv.reader([line], delimiter="\t"))
        if len(fields) != len(header):
            raise DataFormatError(f"Line {number}: expected {len(header)} columns, found {len(fields)}.")
        parsed = []
        for field in fields:
            if field.strip() == "-":
                parsed.append(np.nan)
                continue
            try:
                value = float(field)
            except ValueError as exc:
                raise DataFormatError(f"Line {number}: expected a numeric cell or '-' for missing data.") from exc
            if not np.isfinite(value):
                raise DataFormatError(f"Line {number}: non-finite numeric values are unsupported; only '-' denotes missing cells.")
            parsed.append(value)
        rows.append(parsed)
    raw = np.asarray(rows, dtype=float)
    if len(raw) < 4:
        raise DataFormatError("At least four numeric time rows are required.")
    if np.any(raw[:, 1:][np.isfinite(raw[:, 1:])] < 0):
        raise DataFormatError("Raw photon counts cannot be negative.")
    known_time = np.isfinite(raw[:, 0])
    positions = np.flatnonzero(known_time)
    if len(positions) < 4:
        raise DataFormatError("At least four calibrated timestamps are required.")
    if not np.all(known_time[positions[0]:positions[-1]+1]):
        raise DataFormatError("Missing timestamps occur inside the recorded interval; a complete uniform time axis is required.")
    original = raw[known_time, 0] * _UNIT_SCALE[source_unit]
    if np.any(np.diff(original) <= 0):
        raise DataFormatError("Time must be strictly ascending, without duplicate bins.")
    if not np.isfinite(grid_tolerance_bins) or not 0 < grid_tolerance_bins <= 0.01:
        raise DataFormatError("grid_tolerance_bins must be positive and at most 0.01 of one bin.")
    index = np.arange(len(original), dtype=float)
    step = float(np.dot(index, original-original[0]) / np.dot(index, index))
    canonical = original[0] + index * step
    maximum_displacement = float(np.max(np.abs(original-canonical)))
    displacement_bins = maximum_displacement / step
    if displacement_bins > grid_tolerance_bins:
        raise DataFormatError(f"Time grid is nonuniform: canonicalization would move a timestamp by {displacement_bins:.6g} bins, above {grid_tolerance_bins:g}.")
    counts = raw[known_time, 1:].copy()
    missing_rows = np.flatnonzero(~known_time)
    unknown_time_counts = raw[~known_time, 1:]
    warnings = []
    if source_unit != "ns":
        warnings.append(f"Time header explicitly declares {source_unit}; imported times were converted to ns.")
    if maximum_displacement > step*1e-10:
        warnings.append(f"Rounded export timestamps were canonicalized to a uniform grid; maximum correction {maximum_displacement:.6g} ns ({displacement_bins:.6g} bins). Original timestamps are retained.")
    if len(missing_rows):
        warnings.append(f"Excluded {len(missing_rows)} endpoint row(s) without timestamps, containing {int(np.isfinite(unknown_time_counts).sum())} finite count cells. Raw rows are retained for audit; no time was invented.")
    if np.any(~np.isfinite(counts)):
        warnings.append("Some curve endpoint counts are missing ('-'). Select curves to trim their incomplete endpoints; missing counts are never replaced by zero.")
    warnings.append("Curve labels do not identify IRF, polarization, or wavelength. Assign roles and wavelength metadata explicitly.")
    metadata = {
        "format": "PlotBox ASCII",
        "source_name": name,
        "header": header,
        "source_time_unit": source_unit,
        "output_time_unit": "ns",
        "raw_data_rows": len(raw),
        "retained_timestamp_rows": int(known_time.sum()),
        "curve_count": len(labels),
        "dropped_timestamp_data_rows": missing_rows.tolist(),
        "dropped_timestamp_finite_count_cells": int(np.isfinite(unknown_time_counts).sum()),
        "dropped_timestamp_counts_sum": float(np.nansum(unknown_time_counts)),
        "missing_count_cells_per_curve": np.sum(~np.isfinite(counts), axis=0).astype(int).tolist(),
        "bin_width_ns": step,
        "canonical_grid_anchor_ns": float(original[0]),
        "max_timestamp_adjustment_ns": maximum_displacement,
        "max_timestamp_adjustment_bins": displacement_bins,
        "grid_tolerance_bins": grid_tolerance_bins,
        "curve_roles": "unassigned",
        "wavelengths_nm": None,
    }
    return CurveBundle(canonical, counts, labels, source_unit, name, warnings, metadata,
                       original.copy(), raw[:, 0].copy(), raw[:, 1:].copy())

