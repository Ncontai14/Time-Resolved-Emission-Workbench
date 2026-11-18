"""Checks for the global-fit -> steady-state -> spectral display workflow.

These use independently generated decays and invariants of intensity matching;
display interpolation is checked separately from the reconstructed matrices.
"""
from __future__ import annotations

import unittest
import tempfile
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

from tcspc import FitConfig, fit_global
from tcspc.spectra import load_steady_state, reconstruct_spectra
from tcspc.spectral_display import display_curve


class GlobalSpectralWorkflowChecks(unittest.TestCase):
    def test_exported_normalized_spectrum_can_be_loaded_again(self):
        result = reconstruct_spectra(
            [500, 510, 520], [[1], [2], [1]], [2], [0, 1],
            steady_wavelengths=[490, 500, 510, 520, 530],
            steady_intensities=[-2, 4, 10, 7, -2],
            steady_normalization="minmax")
        with tempfile.TemporaryDirectory() as directory:
            result.export(directory)
            normalized = load_steady_state(Path(directory)/"normalized_steady_state.csv")
            matched = load_steady_state(Path(directory)/"sampled_steady_state.csv")
        np.testing.assert_array_equal(normalized.wavelengths, [490, 500, 510, 520, 530])
        np.testing.assert_allclose(normalized.intensities, [0, .5, 1, .75, 0])
        reused = reconstruct_spectra(
            result.wavelengths, result.amplitudes, result.lifetimes, result.times,
            steady_state=matched, steady_alignment="exact")
        np.testing.assert_allclose(reused.tres, result.tres, rtol=1e-13)
        np.testing.assert_allclose(reused.tranes, result.tranes, rtol=1e-13)

    def test_global_fit_outputs_connect_to_shuffled_normalized_intensities(self):
        # Build independent bin-averaged decays for three wavelength channels.
        # The wavelength and uploaded steady-state orders deliberately differ.
        time = np.arange(240) * 0.05
        tau = np.array([0.8, 3.2])
        waves = np.array([540.0, 500.0, 520.0])
        amplitudes = np.array([[1200, 4200], [4500, 700], [2400, 2000]])
        bin_average = tau / 0.05 * (1 - np.exp(-0.05 / tau))
        kernel = np.exp(-time[:, None] / tau)
        counts = kernel @ (amplitudes * bin_average).T + 2.0
        irf = np.zeros_like(time)
        irf[0] = 1.0
        fitted = fit_global(time, counts, irf=irf, config=FitConfig(
            n_components=2, fixed_lifetimes=tuple(tau), statistic="uniform",
            n_starts=1, fit_background=False, background=2.0,
        ))
        np.testing.assert_allclose(fitted.amplitudes, amplitudes, rtol=1e-8)
        steady = load_steady_state(
            b"wavelength_nm,intensity\n520,0.8\n540,0.6\n500,0.25\n")
        delays = np.array([0.0, 1.0, 8.0])
        result = reconstruct_spectra(
            waves, fitted.amplitudes, fitted.lifetimes, delays,
            steady_state=steady, steady_alignment="exact",
            steady_normalization="none",
        )
        np.testing.assert_array_equal(result.wavelengths, [500, 520, 540])
        # No peak renormalization: a sampled subset is allowed to peak below1.
        np.testing.assert_array_equal(result.steady_state, [.25, .8, .6])
        expected = np.array([
            [ss * sum(a * np.exp(-t / tau)) / sum(a * tau) for t in delays]
            for a, ss in zip(amplitudes[[1, 2, 0]], [.25, .8, .6])
        ])
        np.testing.assert_allclose(result.tres, expected, rtol=1e-9)
        np.testing.assert_allclose(
            trapezoid(result.tranes, result.wavelengths, axis=0), 1, atol=1e-14)
        self.assertEqual(result.metadata["steady_state_alignment"], "exact")

    def test_full_raw_spectrum_is_peak_normalized_before_sampling(self):
        source = load_steady_state(
            b"wavelength_nm,intensity\n490,10\n500,30\n510,50\n520,100\n530,10\n")
        result = reconstruct_spectra(
            [505, 515], [[1], [1]], [2], [0, 1], steady_state=source,
            steady_normalization="peak", steady_alignment="interpolate")
        # The source peak is outside the requested wavelength subset.
        np.testing.assert_allclose(result.steady_state, [.4, .75])
        np.testing.assert_allclose(result.tres[:, 0], [.2, .375])
        self.assertEqual(result.metadata["steady_state_normalization_divisor"], 100)
        np.testing.assert_allclose(
            result.metadata["steady_state_source_normalized_intensities"],
            [.1, .3, .5, 1, .1])

    def test_full_raw_spectrum_minmax_uses_source_extrema_before_sampling(self):
        result = reconstruct_spectra(
            [505, 515], [[1], [1]], [2], [0, 1],
            steady_wavelengths=[490, 500, 510, 520, 530],
            steady_intensities=[10, 30, 50, 100, 10],
            steady_normalization="minmax")
        np.testing.assert_allclose(result.steady_state, [1 / 3, 13 / 18])
        self.assertEqual(result.metadata["steady_state_normalization_offset"], 10)
        self.assertEqual(result.metadata["steady_state_normalization_divisor"], 90)
        np.testing.assert_allclose(
            result.metadata["steady_state_source_normalized_intensities"],
            [0, 2 / 9, 4 / 9, 1, 0])

    def test_minmax_is_invariant_to_positive_affine_raw_intensity_scale(self):
        waves = np.array([500, 510, 530])
        raw = np.array([-5.0, 7.0, 4.0])
        kwargs = dict(wavelengths=waves, amplitudes=[[1], [2], [3]],
                      lifetimes=[2.0], times=[0.0, 3.0],
                      steady_wavelengths=waves, steady_normalization="minmax")
        first = reconstruct_spectra(**kwargs, steady_intensities=raw)
        changed_scale = reconstruct_spectra(**kwargs, steady_intensities=3 * raw + 11)
        np.testing.assert_allclose(first.tres, changed_scale.tres, rtol=1e-13)
        np.testing.assert_allclose(first.tranes, changed_scale.tranes, rtol=1e-13)

    def test_signed_raw_intensities_require_explicit_minmax(self):
        content = b"wavelength_nm,intensity\n500,-7\n510,-3\n520,-5\n"
        with self.assertRaises(ValueError):
            load_steady_state(content)
        raw = load_steady_state(content, allow_negative=True)
        result = reconstruct_spectra(
            [500, 510, 520], [[1], [1], [1]], [2], [0, 1],
            steady_state=raw, steady_normalization="minmax", steady_alignment="exact")
        np.testing.assert_allclose(result.steady_state, [0, 1, .5])
        self.assertEqual(result.metadata["steady_state_normalization_offset"], -7)
        for normalization in ["none", "peak", "area"]:
            with self.subTest(normalization=normalization), self.assertRaises(ValueError):
                reconstruct_spectra(
                    [500, 510, 520], [[1], [1], [1]], [2], [0],
                    steady_state=raw, steady_normalization=normalization)

    def test_constant_intensity_cannot_be_minmax_normalized(self):
        for intensity in [0.0, 2.0, -1.0]:
            with self.subTest(intensity=intensity), self.assertRaises(ValueError):
                reconstruct_spectra(
                    [500, 520], [[1], [1]], [2], [0],
                    steady_state=[intensity, intensity], steady_normalization="minmax")

    def test_exact_matching_rejects_missing_extra_duplicate_and_shifted_rows(self):
        cases = [
            ([500, 540], [.2, .4]),
            ([490, 500, 520, 540], [.1, .2, .3, .4]),
            ([500, 520, 520], [.2, .3, .4]),
            ([500, 520.01, 540], [.2, .3, .4]),
        ]
        for waves, intensity in cases:
            with self.subTest(waves=waves), self.assertRaises(ValueError):
                reconstruct_spectra(
                    [500, 520, 540], [[1], [1], [1]], [2], [0],
                    steady_wavelengths=waves, steady_intensities=intensity,
                    steady_alignment="exact")

    def test_exact_matching_tolerates_rounding_without_interpolation(self):
        result = reconstruct_spectra(
            [500, 520], [[1], [1]], [2], [0],
            steady_wavelengths=[520 + 5e-7, 500 - 5e-7],
            steady_intensities=[.7, .2], steady_alignment="exact")
        np.testing.assert_array_equal(result.wavelengths, [500, 520])
        np.testing.assert_array_equal(result.steady_state, [.2, .7])
        self.assertFalse(result.metadata["spectral_extrapolation"])

    def test_exact_matching_rejects_ambiguous_nearby_wavelengths(self):
        with self.assertRaises(ValueError):
            reconstruct_spectra(
                [500, 500 + 5e-7], [[1], [1]], [2], [0],
                steady_wavelengths=[500, 500 + 5e-7], steady_intensities=[.2, .7],
                steady_alignment="exact")


class SpectralDisplayChecks(unittest.TestCase):
    def test_all_curve_styles_retain_nodes_and_stay_within_wavelength_coverage(self):
        waves = np.array([500.0, 510.0, 520.0, 530.0])
        values = np.array([0.0, 1.0, 0.0, 0.0])
        for style in ["linear", "pchip", "bspline"]:
            with self.subTest(style=style):
                x, y = display_curve(waves, values, style=style, points=127)
                self.assertEqual(x[0], waves[0])
                self.assertEqual(x[-1], waves[-1])
                self.assertTrue(np.all(np.diff(x) > 0))
                np.testing.assert_allclose(y[np.searchsorted(x, waves)], values, atol=1e-13)

    def test_shape_preserving_style_avoids_negative_artifacts_while_bspline_is_unclipped(self):
        waves = np.array([500.0, 510.0, 520.0, 530.0])
        values = np.array([0.0, 1.0, 0.0, 0.0])
        _, pchip = display_curve(waves, values, style="pchip", points=500)
        _, spline = display_curve(waves, values, style="bspline", points=500)
        self.assertGreaterEqual(pchip.min(), -1e-14)
        self.assertLessEqual(pchip.max(), 1 + 1e-14)
        # The cubic interpolant through these nodes is negative near525nm.
        self.assertLess(spline.min(), -.25)

    def test_display_selection_does_not_mutate_spectral_matrices_or_area_normalization(self):
        result = reconstruct_spectra(
            [500, 510, 520, 530], [[1, 3], [3, 1], [1, 2], [2, 1]],
            [1.0, 4.0], [0, 5], steady_state=[0, 1, .1, .2])
        original_tres, original_tranes = result.tres.copy(), result.tranes.copy()
        for style in ["linear", "pchip", "bspline"]:
            for j in range(len(result.times)):
                display_curve(result.wavelengths, result.tres[:, j], style=style)
                display_curve(result.wavelengths, result.tranes[:, j], style=style)
        np.testing.assert_array_equal(result.tres, original_tres)
        np.testing.assert_array_equal(result.tranes, original_tranes)
        np.testing.assert_allclose(
            trapezoid(result.tranes, result.wavelengths, axis=0), 1, atol=1e-14)

    def test_two_wavelengths_have_a_linear_connection_for_every_style(self):
        for style in ["linear", "pchip", "bspline"]:
            with self.subTest(style=style):
                x, y = display_curve([500, 530], [.2, .8], style=style, points=17)
                np.testing.assert_allclose(y, .2 + .6 * (x - 500) / 30, atol=1e-13)

    def test_invalid_display_coordinates_values_and_style_are_rejected(self):
        cases = [([510, 500], [1, 2]), ([500, 500], [1, 2]),
                 ([500, 510], [-1, 2]), ([500, 510], [np.nan, 2])]
        for waves, values in cases:
            with self.subTest(waves=waves, values=values), self.assertRaises(ValueError):
                display_curve(waves, values)
        with self.assertRaises(ValueError):
            display_curve([500, 510], [1, 2], style="unknown")


if __name__ == "__main__":
    unittest.main()
