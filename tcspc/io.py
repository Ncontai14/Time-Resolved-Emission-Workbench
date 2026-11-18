"""Validated text import for raw TCSPC histograms and TRES matrices.

Time units are explicit, never inferred from labels. Native binary files
are deliberately unsupported: use an ASCII export with time in the first column.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import csv
import io
import re
from typing import BinaryIO, TextIO

import numpy as np
from .time_axis import channel_grid


class DataFormatError(ValueError):
    """An input does not describe a supported raw photon-counting dataset."""


@dataclass(frozen=True)
class DecayData:
    time_ns: np.ndarray
    counts: np.ndarray
    name: str = "Decay"
    metadata: dict = field(default_factory=dict)


@dataclass(frozen=True)
class TRESData:
    time_ns: np.ndarray
    counts: np.ndarray
    wavelengths: np.ndarray
    metadata: dict = field(default_factory=dict)


_UNITS_TO_NS = {"ps": 1e-3, "ns": 1.0, "us": 1e3, "ms": 1e6, "s": 1e9}
_NATIVE_EXTENSIONS = {".etc", ".ete", ".etr", ".etw", ".exe"}
_NUMERIC = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][+-]?\d+)?$")


def _time_scale(time_unit: str) -> float:
    unit = time_unit.strip().lower().replace("µ", "u").replace("μ", "u")
    try:
        return _UNITS_TO_NS[unit]
    except KeyError as exc:
        raise DataFormatError("Choose an explicit time unit: ps, ns, us, ms, or s.") from exc


def _read_source(source: str | Path | bytes | BinaryIO | TextIO) -> tuple[str, str]:
    if isinstance(source, (str, Path)):
        path = Path(source)
        name = path.name
        content = path.read_bytes()
    elif isinstance(source, bytes):
        name, content = "Uploaded data", source
    elif hasattr(source, "read"):
        name = Path(str(getattr(source, "name", "Uploaded data"))).name
        # Do not change an upload's current position for subsequent plots/imports.
        if hasattr(source, "getvalue"):
            content = source.getvalue()
        else:
            position = source.tell() if hasattr(source, "tell") else None
            content = source.read()
            if position is not None and hasattr(source, "seek"):
                source.seek(position)
    else:
        raise TypeError("source must be a path, bytes, or a readable file-like object")
    if Path(name).suffix.lower() in _NATIVE_EXTENSIONS:
        raise DataFormatError(
            "Native binary files and executables cannot be imported. Export the "
            "curve(s) as ASCII/text with time in the first column, then load that export."
        )
    if isinstance(content, bytes):
        if b"\x00" in content[:256]:
            raise DataFormatError(
                "This is binary data, not a text histogram. Export the measurement "
                "curve(s) as ASCII/text with time in the first column."
            )
        try:
            content = content.decode("utf-8-sig")
        except UnicodeDecodeError as exc:
            raise DataFormatError("Use a UTF-8/ASCII text export with decimal points.") from exc
    return content.lstrip("\ufeff"), name


def _tokens(line: str, delimiter: str | None) -> list[str]:
    if delimiter == "whitespace":
        return line.split()
    selected = delimiter
    if selected is None:
        # Tabs take priority, since commas may be decimal separators in a bad export.
        selected = "\t" if "\t" in line else ";" if ";" in line else "," if "," in line else None
    if selected is None:
        return line.split()
    if selected not in (",", ";", "\t", " "):
        raise DataFormatError("delimiter must be comma, semicolon, tab, space, or 'whitespace'.")
    return [token.strip() for token in next(csv.reader([line], delimiter=selected))]


def _numeric_row(tokens: list[str]) -> list[float] | None:
    if not tokens or any(not _NUMERIC.fullmatch(token) for token in tokens):
        return None
    return [float(token.replace("D", "E").replace("d", "e")) for token in tokens]


def _table(source, delimiter: str | None) -> tuple[np.ndarray, list[list[str]], str]:
    text, name = _read_source(source)
    if text.lstrip().startswith("PlotBox Export"):
        raise DataFormatError(
            "This is a PlotBox export. Use tcspc.load_plotbox() to read its "
            "header time units, missing cells, and multiple curves, then select curve roles explicitly."
        )
    rows: list[list[float]] = []
    headers: list[list[str]] = []
    width: int | None = None
    for number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(("#", "//", "%", ";")):
            continue
        # Trailing comments are accepted only after a whitespace-separated # marker.
        line = re.split(r"\s+#", line, maxsplit=1)[0].strip()
        tokens = _tokens(line, delimiter)
        if any(
            token.strip().lower() in {"nan", "+nan", "-nan", "inf", "+inf", "-inf", "infinity"}
            for token in tokens
        ):
            raise DataFormatError(f"Line {number}: non-finite values are not allowed.")
        numeric = _numeric_row(tokens)
        if numeric is None:
            # A comma within a tab/semicolon-separated token is an ambiguous decimal comma.
            if any(re.fullmatch(r"[+-]?\d+,\d+(?:[eE][+-]?\d+)?", token) for token in tokens):
                raise DataFormatError(f"Line {number}: decimal commas are unsupported; export decimal points.")
            if rows:
                raise DataFormatError(f"Line {number}: expected numeric data; check missing values or delimiters.")
            headers.append(tokens)
            continue
        if width is None:
            width = len(numeric)
        if len(numeric) != width:
            raise DataFormatError(
                f"Line {number}: found {len(numeric)} columns; expected {width}. "
                "Check missing values and the selected delimiter."
            )
        rows.append(numeric)
    if len(rows) < 3:
        raise DataFormatError("At least three numeric rows are required.")
    return np.asarray(rows, dtype=float), headers, name


def _validate(time_ns: np.ndarray, counts: np.ndarray) -> tuple[np.ndarray, dict]:
    if not np.all(np.isfinite(time_ns)) or not np.all(np.isfinite(counts)):
        raise DataFormatError("Time and raw counts must contain only finite numbers.")
    if np.any(counts < 0):
        raise DataFormatError("Raw photon counts cannot be negative; import unsubtracted count histograms.")
    dt = np.diff(time_ns)
    if np.any(dt <= 0):
        raise DataFormatError("Time must be strictly ascending, with no duplicate bins.")
    try:
        grid, audit = channel_grid(time_ns)
    except ValueError as error:
        raise DataFormatError(str(error)) from error
    if not np.any(counts > 0):
        raise DataFormatError("The dataset contains no positive photon counts.")
    return grid, audit


def load_decay(
    source,
    *,
    time_unit: str = "ns",
    time_column: int | None = 0,
    count_column: int = 1,
    delimiter: str | None = None,
    bin_width: float | None = None,
    name: str | None = None,
) -> DecayData:
    """Read a text histogram; ``bin_width`` is in ``time_unit`` for count-only input.

    ``time_column=None`` selects count-only input and requires a positive bin width.
    For a one-column file use ``count_column=0``. Column indices are zero-based.
    Extra numeric columns are retained only when explicitly selected.
    """
    scale = _time_scale(time_unit)
    table, _, source_name = _table(source, delimiter)
    if count_column < 0 or count_column >= table.shape[1]:
        raise DataFormatError(f"Count column {count_column} is outside the {table.shape[1]} input columns.")
    if time_column is None:
        if bin_width is None or not np.isfinite(bin_width) or bin_width <= 0:
            raise DataFormatError("Count-only input requires a finite, positive bin_width in the selected time unit.")
        time = np.arange(table.shape[0], dtype=float) * bin_width * scale
    else:
        if time_column < 0 or time_column >= table.shape[1] or time_column == count_column:
            raise DataFormatError("Select different valid time and count columns (zero-based indices).")
        time = table[:, time_column] * scale
    counts = table[:, count_column].copy()
    time, audit = _validate(time, counts)
    return DecayData(time, counts, name or source_name, audit)


def load_tres(
    source,
    *,
    time_unit: str = "ns",
    wavelengths=None,
    delimiter: str | None = None,
) -> TRESData:
    """Read rows ``time, counts_at_lambda1, counts_at_lambda2, ...``.

    Supply wavelengths in nm explicitly, or use a header such as
    ``time_ns,500,520,540``. Wavelengths are never invented from channel indices.
    """
    scale = _time_scale(time_unit)
    table, headers, _ = _table(source, delimiter)
    if table.shape[1] < 2:
        raise DataFormatError("TRES input needs a time column and at least one wavelength column.")
    n_trace = table.shape[1] - 1
    if wavelengths is None:
        for tokens in reversed(headers):
            if len(tokens) != table.shape[1]:
                continue
            candidate = _numeric_row(tokens[1:])
            if candidate is not None:
                wavelengths = candidate
                break
    if wavelengths is None:
        raise DataFormatError(
            "Provide wavelengths in nm or a header such as 'time_ns,500,520,540'. "
            "The first column must contain time."
        )
    try:
        labels = np.asarray(wavelengths, dtype=float)
    except (TypeError, ValueError) as exc:
        raise DataFormatError("Wavelengths must be a numeric sequence in nm.") from exc
    if labels.ndim != 1 or len(labels) != n_trace:
        raise DataFormatError(f"Expected {n_trace} wavelengths, one per count column.")
    if not np.all(np.isfinite(labels)) or np.any(labels <= 0) or len(np.unique(labels)) != len(labels):
        raise DataFormatError("Wavelengths must be finite, positive, and unique.")
    time, counts = table[:, 0] * scale, table[:, 1:].copy()
    time, audit = _validate(time, counts)
    if np.any(np.sum(counts, axis=0) <= 0):
        raise DataFormatError("Every TRES channel must contain positive photon counts; remove empty channels.")
    return TRESData(time, counts, labels.copy(), audit)
