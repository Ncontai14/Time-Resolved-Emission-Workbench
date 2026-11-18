"""Independent checks of steady-state-constrained TRES/TRANES reconstruction."""
from pathlib import Path
import io
import json
import tempfile
import unittest

import numpy as np
from scipy.integrate import simpson, trapezoid

from tcspc.spectra import load_fitted_spectra, load_steady_state, reconstruct_spectra


class SpectralReconstructionChecks(unittest.TestCase):
    def setUp(self):
        source = Path(__file__).resolve().parents[1]/"examples"/"fitted_three_components.csv"
        self.fitted = load_fitted_spectra(source)
        self.tau = np.array([114.15, 17.415, 4.613])
        self.times = np.array([0, 0.025, 1, 10, 50, 100])

    def reconstruct(self, **kwargs):
        f = self.fitted
        return reconstruct_spectra(f.wavelengths, f.amplitudes, self.tau,
                                   self.times, f.steady_state, **kwargs)

    def test_uploaded_three_component_data_matches_independent_notebook_formula(self):
        f = self.fitted
        self.assertEqual(f.amplitudes.shape, (13, 3))
        self.assertEqual(f.component_names, ("A1", "A2", "A3"))
        expected = np.array([
            [ss*sum(a*np.exp(-time/self.tau))/sum(a*self.tau) for time in self.times]
            for a, ss in zip(f.amplitudes, f.steady_state)
        ])
        result = self.reconstruct()
        np.testing.assert_allclose(result.tres, expected, rtol=1e-12, atol=1e-18)
        np.testing.assert_allclose(result.tranes, expected/trapezoid(expected, f.wavelengths, axis=0), rtol=1e-12)
        np.testing.assert_allclose(trapezoid(result.tranes, result.wavelengths, axis=0), 1, atol=1e-14)

    def test_time_integral_recovers_steady_state(self):
        f = self.fitted
        time = np.linspace(0, 20*self.tau.max(), 20001)
        result = reconstruct_spectra(f.wavelengths, f.amplitudes, self.tau, time, f.steady_state)
        np.testing.assert_allclose(simpson(result.tres, x=time, axis=1), f.steady_state, rtol=2e-7)

    def test_per_wavelength_amplitude_scale_cancels(self):
        f = self.fitted
        scales = np.geomspace(1e-40, 1e40, len(f.wavelengths))[:, None]
        result = reconstruct_spectra(f.wavelengths, f.amplitudes*scales, self.tau, self.times, f.steady_state)
        np.testing.assert_allclose(result.tres, self.reconstruct().tres, rtol=1e-12)

    def test_peak_normalization_and_linear_sampling_of_full_spectrum(self):
        source = load_steady_state(b"wavelength_nm,intensity\n490,2\n500,4\n510,8\n520,6\n530,0\n")
        result = reconstruct_spectra([500, 510, 525], np.ones((3,1)), [2], [0, 1],
                                     steady_state=source, steady_normalization="peak")
        np.testing.assert_allclose(result.steady_state, [0.5, 1, 0.375])
        np.testing.assert_allclose(result.tres[:,0], np.array([0.5, 1, 0.375])/2)
        self.assertEqual(result.metadata["steady_state_normalization_divisor"], 8)
        self.assertFalse(result.metadata["spectral_extrapolation"])
        with self.assertRaisesRegex(ValueError, "cover"):
            reconstruct_spectra([480, 510], np.ones((2,1)), [2], [0], steady_state=source)

    def test_global_steady_state_scale_changes_tres_but_not_tranes(self):
        f = self.fitted
        original = self.reconstruct()
        scaled = reconstruct_spectra(f.wavelengths, f.amplitudes, self.tau, self.times, f.steady_state*5)
        peak = self.reconstruct(steady_normalization="peak")
        area = self.reconstruct(steady_normalization="area")
        np.testing.assert_allclose(scaled.tres, original.tres*5, rtol=1e-12)
        for result in [scaled, peak, area]:
            np.testing.assert_allclose(result.tranes, original.tranes, rtol=1e-12)

    def test_sorting_preserves_alignment_and_does_not_mutate_inputs(self):
        f = self.fitted
        waves = f.wavelengths[::-1].copy()
        aa = f.amplitudes[::-1].copy()
        before = aa.copy()
        result = reconstruct_spectra(waves, aa, self.tau, self.times, f.steady_state[::-1])
        np.testing.assert_allclose(result.tres, self.reconstruct().tres, rtol=1e-12)
        np.testing.assert_array_equal(aa, before)

    def test_tail_reference_conversion_recovers_excitation_model(self):
        f = self.fitted
        delay = 8.0
        aa_at_reference = f.amplitudes*np.exp(-delay/self.tau)
        result = reconstruct_spectra(f.wavelengths, aa_at_reference, self.tau,
                                     self.times, f.steady_state,
                                     amplitude_reference_delay_ns=delay)
        np.testing.assert_allclose(result.tres, self.reconstruct().tres, rtol=1e-12)
        self.assertEqual(result.metadata["amplitude_reference_delay_ns"], delay)

    def test_delimiter_headers_and_explicit_headerless_combined_table(self):
        for sep in [",", "\t", " "]:
            text = sep.join(["wavelength", "A1", "SS"])+"\n"+"\n".join(sep.join(map(str,row)) for row in [(510,2,0.2),(500,1,0.1)])
            result = load_fitted_spectra(text.encode())
            np.testing.assert_array_equal(result.wavelengths, [500,510])
            np.testing.assert_allclose(result.steady_state, [.1,.2])
        bare = b"500 1 .1\n510 2 .2\n"
        self.assertEqual(load_fitted_spectra(bare).amplitudes.shape[1], 2)
        self.assertEqual(load_fitted_spectra(bare, has_steady_state=True).amplitudes.shape[1], 1)
        bad = b"wave,A1,SS\n500,1,-\n510,2,.2\n520,2,.3\n"
        with self.assertRaises(ValueError):
            load_fitted_spectra(bad)

    def test_invalid_and_unsupported_inputs_raise_informative_errors(self):
        cases = [
            ([500,500], [[1],[1]], [2], [0], [1,1]),
            ([500,510], [[-1],[1]], [2], [0], [1,1]),
            ([500,510], [[1],[1]], [0], [0], [1,1]),
            ([500,510], [[1],[1]], [2], [-1], [1,1]),
            ([500,510], [[0],[1]], [2], [0], [1,1]),
            ([500,510], [[1],[1]], [2], [0], [0,0]),
            ([500,510], [[1],[1]], [2], [0], None),
        ]
        for args in cases:
            with self.subTest(args=args):
                with self.assertRaises(ValueError):
                    reconstruct_spectra(*args)
        with self.assertRaisesRegex(ValueError,"zero or non-finite"):
            reconstruct_spectra([500,510], [[1],[1]], [1], [1e6], [1,1])

    def test_zero_steady_state_and_zero_amplitudes_are_explicitly_valid(self):
        result = reconstruct_spectra([500,510,520], [[0],[1],[1]], [2], [0,1], [0,1,1])
        np.testing.assert_array_equal(result.tres[0], 0)
        np.testing.assert_allclose(trapezoid(result.tranes,result.wavelengths,axis=0),1)

    def test_export_contains_tres_tranes_and_normalization_audit(self):
        result = self.reconstruct()
        with tempfile.TemporaryDirectory() as folder:
            result.export(folder)
            path = Path(folder)
            table = np.loadtxt(path/"tres.csv",delimiter=",",skiprows=1)
            np.testing.assert_allclose(table[:,1:], result.tres, rtol=1e-12)
            metadata = json.loads((path/"spectral_metadata.json").read_text())
            self.assertEqual(metadata["lifetimes_ns_in_amplitude_column_order"], self.tau.tolist())
            self.assertFalse(metadata["spectral_smoothing_applied"])
            self.assertTrue((path/"tranes.csv").exists())


if __name__ == "__main__":
    unittest.main()
