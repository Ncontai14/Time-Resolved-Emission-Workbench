"""Interpret rounded channel labels without resampling photon histograms."""
import numpy as np


def channel_grid(time, *, min_bins=3):
    original = np.asarray(time, dtype=float)
    if original.ndim != 1 or len(original) < min_bins or not np.all(np.isfinite(original)):
        raise ValueError(f"Time must contain at least {min_bins} finite channel labels.")
    if np.any(np.diff(original) <= 0):
        raise ValueError("Time must be strictly ascending, without duplicate channel labels.")
    index = np.arange(len(original), dtype=float)
    width = float(np.dot(index, original-original[0]) / np.dot(index, index))
    grid = original[0] + index*width
    displacement = float(np.max(np.abs(original-grid)))
    fraction = displacement/width
    if fraction > 1e-4:
        raise ValueError(
            "The time labels differ from one fixed-width channel sequence by "
            f"{fraction:.6g} bins, beyond the export-rounding tolerance (0.0001 bins). "
            "Check for missing channels, an incorrect time column, or nonuniform acquisition. "
            "Photon counts have not been altered."
        )
    adjusted = fraction > 1e-10
    return (grid if adjusted else original.copy()), {
        "bin_width_ns": width, "rounded_channel_labels": adjusted,
        "max_timestamp_adjustment_ns": displacement if adjusted else 0.0,
        "max_timestamp_adjustment_bins": fraction if adjusted else 0.0,
        "original_time_ns": original.tolist(),
        "photon_counts_resampled": False,
    }
