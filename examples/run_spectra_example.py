"""Reconstruct the contributed three-component TRES/TRANES example.

Run from the project folder: python examples/run_spectra_example.py
The steady-state values belong to this fitted table only, not the separate
51-wavelength measured-decay example. No fitting is performed by this script.
"""
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid

from tcspc.spectra import load_fitted_spectra, reconstruct_spectra


def main():
    fitted = load_fitted_spectra(ROOT / "examples" / "fitted_three_components.csv")
    # Preserve the notebook's amplitude-column order: A1, A2, A3.
    lifetimes_ns = [114.15, 17.415, 4.613]
    delays_ns = [0, 0.025, 1, 10, 50, 100]
    result = reconstruct_spectra(
        fitted.wavelengths, fitted.amplitudes, lifetimes_ns, delays_ns,
        steady_state=fitted.steady_state, steady_normalization="none",
    )
    output = result.export(ROOT / "example_results" / "spectral_reconstruction")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))
    for j, delay in enumerate(result.times):
        axes[0].plot(result.wavelengths, result.tres[:, j], "o-", ms=3,
                     label=f"{delay:g} ns")
        axes[1].plot(result.wavelengths, result.tranes[:, j], "o-", ms=3,
                     label=f"{delay:g} ns")
    axes[0].set(xlabel="Wavelength (nm)", ylabel="Relative intensity / ns",
                title="Reconstructed TRES")
    axes[1].set(xlabel="Wavelength (nm)", ylabel="Area-normalized intensity (1/nm)",
                title="TRANES")
    for ax in axes:
        ax.legend(fontsize="small", ncol=2)
    fig.tight_layout()
    fig.savefig(output / "tres_tranes.png", dpi=180)
    plt.close(fig)
    print(f"Exported matrices, metadata and plot to {output}")
    print("TRANES spectral areas:", trapezoid(result.tranes, result.wavelengths, axis=0))


if __name__ == "__main__":
    main()
