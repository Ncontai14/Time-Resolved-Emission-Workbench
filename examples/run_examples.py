"""Reproduce the bundled single-decay, synthetic and measured global fits.

Run from the project folder: python examples/run_examples.py
"""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np
from tcspc import FitConfig, fit_decay, fit_global
from tcspc.io import load_decay, load_tres


def main():
    examples = ROOT / "examples"
    out = ROOT / "example_results"
    data = load_decay(examples / "Decay.dat")
    irf = load_decay(examples / "IRF.dat")
    single = fit_decay(data.time_ns, data.counts, irf=irf.counts,
                       irf_time=irf.time_ns,
                       config=FitConfig(irf_background=1.4, irf_start=3.8,
                                        irf_end=8, fit_shift=True))
    single.export(out / "single_decay")
    summaries = {"single_decay": single.lifetimes.tolist()}
    for label, data_name, irf_name, config in [
        ("synthetic_global", "synthetic_TRES.csv", "synthetic_IRF.csv",
         FitConfig(n_components=2, initial_lifetimes=(0.7, 3.3))),
        ("measured_global", "measured_TRES.csv", "measured_IRF.csv",
         FitConfig(n_components=2, initial_lifetimes=(2.5, 4.0),
                   irf_background=0.38, irf_start=3.2, irf_end=6.5,
                   fit_shift=True, shift_bounds=(-0.5, 0.5))),
    ]:
        data = load_tres(examples / data_name)
        irf = load_decay(examples / irf_name)
        result = fit_global(data.time_ns, data.counts, irf=irf.counts,
                            irf_time=irf.time_ns, config=config)
        result.export(out / label)
        np.savetxt(out / label / "wavelengths.csv",
                   np.column_stack([np.arange(1, len(data.wavelengths)+1), data.wavelengths]),
                   delimiter=",", header="trace,wavelength_nm", comments="", fmt="%g")
        (out / label / "import_metadata.json").write_text(
            json.dumps({"decay_time_axis": data.metadata,
                        "irf_time_axis": irf.metadata}, indent=2)+"\n")
        summaries[label] = result.lifetimes.tolist()
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
