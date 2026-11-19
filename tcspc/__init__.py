"""Time-Resolved Emission Workbench for exponential decay fitting and TRES."""
from ._release import __version__
from .core import FitConfig, FitResult, fit_decay, fit_global, predict
from .plotbox import CurveBundle, is_plotbox_export, load_plotbox

__all__ = ["FitConfig", "FitResult", "fit_decay", "fit_global", "predict",
           "CurveBundle", "is_plotbox_export", "load_plotbox"]

from .spectra import (FittedSpectra, SteadyStateSpectrum, SpectralResult,
                      load_fitted_spectra, load_steady_state, reconstruct_spectra)
__all__ += ["FittedSpectra", "SteadyStateSpectrum", "SpectralResult",
            "load_fitted_spectra", "load_steady_state", "reconstruct_spectra"]
