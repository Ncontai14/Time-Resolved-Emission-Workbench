"""Connections between calculated spectral samples for display only.

The interpolation here does not alter reconstruction matrices, wavelength-area
normalization, or exported sampled values. Cubic B-splines may overshoot their
nodes, including below zero; those values are preserved for an honest display.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import PchipInterpolator, make_interp_spline


def display_curve(wavelengths, values, style="pchip", points=500):
    """Return a plotted spectral connection on a dense wavelength grid.

    ``style`` is ``'linear'``, ``'pchip'``, or ``'bspline'``. Original wavelengths
    are included in the returned grid. Interpolation is confined to the measured
    wavelength interval. Inputs must be one-dimensional, strictly increasing
    positive wavelengths and finite nonnegative intensity samples.
    """
    waves = np.asarray(wavelengths, dtype=float)
    yy = np.asarray(values, dtype=float)
    if waves.ndim != 1 or len(waves) < 2 or not np.all(np.isfinite(waves)) or np.any(waves <= 0) or np.any(np.diff(waves) <= 0):
        raise ValueError("Display wavelengths must contain at least two finite, positive, strictly increasing samples.")
    if yy.shape != waves.shape or not np.all(np.isfinite(yy)) or np.any(yy < 0):
        raise ValueError("Display intensities must be finite, nonnegative, and match the wavelength samples.")
    if style not in {"linear", "pchip", "bspline"}:
        raise ValueError("Spectral curve style must be 'linear', 'pchip', or 'bspline'.")
    if isinstance(points, (bool, np.bool_)) or not isinstance(points, (int, np.integer)) or points < 2:
        raise ValueError("Display sampling needs an integer number of points of at least two.")
    grid = np.unique(np.concatenate([np.linspace(waves[0], waves[-1], points), waves]))
    if style == "linear":
        connected = np.interp(grid, waves, yy)
    elif style == "pchip":
        connected = PchipInterpolator(waves, yy, extrapolate=False)(grid)
    else:
        connected = make_interp_spline(waves, yy, k=min(3, len(waves)-1))(grid)
    # Preserve the calculated nodes exactly; only the connections are inferred.
    connected[np.searchsorted(grid, waves)] = yy
    return grid, connected
