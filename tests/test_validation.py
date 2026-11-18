import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from tcspc import FitConfig, fit_decay, fit_global, predict
from tcspc.io import load_decay, load_tres


class RegressionChecks(unittest.TestCase):
    def setUp(self):
        self.time = np.arange(100) * 0.1
        self.irf = np.zeros(100)
        self.irf[5] = 1
        self.decay = np.convolve(self.irf, 800*np.exp(-self.time/1.8))[:100] + 2

    def test_irf_support_removes_out_of_pulse_noise_without_cropping_decay(self):
        noisy = self.irf.copy()
        noisy[70] = 0.03
        r = fit_decay(self.time, self.decay, irf=noisy,
                      config=FitConfig(statistic="uniform", irf_start=0.4, irf_end=0.6, n_starts=1))
        self.assertAlmostEqual(r.lifetimes[0], 1.8, places=5)
        self.assertTrue(np.all(r.fit_mask))

    def test_invalid_shift_bounds_raise_value_error(self):
        with self.assertRaises(ValueError):
            fit_decay(self.time, self.decay, irf=self.irf,
                      config=FitConfig(fit_shift=True, shift_bounds=(0,)))

    def test_plain_loaders_reject_plotbox_instead_of_silently_scaling_seconds_as_ns(self):
        source = b"PlotBox Export\nTime[s]X\tcrv[0] []\n0\t1\n1e-9\t2\n2e-9\t3\n3e-9\t1\n"
        for loader in [load_decay, load_tres]:
            with self.assertRaisesRegex(ValueError, "load_plotbox"):
                loader(source)

    def test_components_can_vanish_at_some_wavelengths_while_remaining_globally_supported(self):
        time = np.arange(200) * 0.08
        irf = np.zeros(len(time))
        irf[8] = 1
        counts = np.column_stack([
            np.convolve(irf, amplitude*np.exp(-time/tau))[:len(time)] + 2
            for tau, amplitude in [(0.8, 2000)]*4 + [(3.2, 4000)]*4
        ])
        result = fit_global(time, counts, irf=irf,
                            config=FitConfig(n_components=2, statistic="uniform",
                                             initial_lifetimes=(0.8, 3.2), n_starts=1))
        np.testing.assert_allclose(result.lifetimes, [0.8, 3.2], rtol=1e-4)
        self.assertFalse(any("unsupported" in note for note in result.warnings))

    def test_fixed_parameters_do_not_mutate_initial_guess_and_export_is_json(self):
        initial = np.array([2.0])
        r = fit_decay(self.time, self.decay, irf=self.irf,
                      config=FitConfig(initial_lifetimes=initial, fixed_lifetimes=(1.8,), n_starts=1))
        self.assertEqual(initial[0], 2.0)
        with tempfile.TemporaryDirectory() as folder:
            r.export(folder)
            content = (Path(folder)/"result.json").read_text()
            decoded = json.loads(content, parse_constant=lambda x: self.fail(x))
            self.assertEqual(decoded["time_unit"], "ns")
            self.assertTrue((Path(folder)/"curves.csv").exists())


if __name__ == "__main__":
    unittest.main()
