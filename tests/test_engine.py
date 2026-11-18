"""Independent numerical checks for the fitting engine.

Synthetic truth is built here with numpy convolution, never the engine predictor.
Count-noise tests use fixed random seeds and deliberately resolvable components.
"""
from __future__ import annotations

import unittest
import numpy as np

from tcspc import FitConfig, fit_decay, fit_global, predict


def sampled_irf(time, center=1.4, width=0.12):
    h = np.exp(-0.5 * ((time - center) / width) ** 2)
    return h / h.sum()


def independent_decay(time, irf, lifetimes, amplitudes, background=0.0):
    """Unaveraged discrete decay, whose amplitudes absorb bin averaging."""
    elapsed = np.arange(len(time)) * (time[1] - time[0])
    signal = sum(a * np.exp(-elapsed / tau)
                 for a, tau in zip(amplitudes, lifetimes))
    return np.convolve(irf / np.sum(irf), signal)[:len(time)] + background


class EngineNumerics(unittest.TestCase):
    def setUp(self):
        self.time = np.arange(400) * 0.05
        self.irf = sampled_irf(self.time)

    def test_predict_is_linear_causal_convolution(self):
        tau, amplitude, background = 1.8, 700.0, 3.0
        dt = self.time[1] - self.time[0]
        bin_average = tau / dt * (-np.expm1(-dt / tau))
        expected = independent_decay(self.time, self.irf, [tau],
                                     [amplitude * bin_average], background)
        got = predict(self.time, [tau], [amplitude], backgrounds=background,
                      irf=self.irf)
        self.assertEqual(got.shape, (len(self.time), 1))
        np.testing.assert_allclose(got[:, 0], expected, rtol=1e-10, atol=1e-9)

    def test_noiseless_single_decay_recovers_lifetime_and_background(self):
        counts = independent_decay(self.time, self.irf, [2.35], [5000], 7.0)
        result = fit_decay(self.time, counts, irf=self.irf, config=FitConfig(
            statistic='uniform', initial_lifetimes=(1.6,), n_starts=1))
        self.assertAlmostEqual(result.lifetimes[0], 2.35, delta=2e-4)
        self.assertAlmostEqual(result.backgrounds[0], 7.0, delta=2e-3)
        np.testing.assert_allclose(result.fitted[:, 0], counts, rtol=1e-5, atol=1e-3)

    def test_sparse_poisson_counts_with_zeros_remain_finite(self):
        rng = np.random.default_rng(20261116)
        mean = independent_decay(self.time, self.irf, [1.7], [75], 0.02)
        counts = rng.poisson(mean)
        self.assertTrue(np.any(counts == 0))
        result = fit_decay(self.time, counts, irf=self.irf, config=FitConfig(
            statistic='poisson', initial_lifetimes=(2.0,), n_starts=1))
        self.assertTrue(np.isfinite(result.fitted).all())
        self.assertTrue(np.isfinite(result.residuals[result.fit_mask]).all())
        self.assertAlmostEqual(result.lifetimes[0], 1.7, delta=0.18)
        self.assertGreaterEqual(result.backgrounds[0], 0.0)

    def test_global_two_component_lifetimes_are_shared(self):
        lifetimes = np.array([0.65, 3.4])
        amplitudes = np.array([[8000, 900], [2500, 5500], [800, 9000]])
        backgrounds = np.array([2.0, 5.0, 8.0])
        counts = np.column_stack([
            independent_decay(self.time, self.irf, lifetimes, a, b)
            for a, b in zip(amplitudes, backgrounds)])
        result = fit_global(self.time, counts, irf=self.irf, config=FitConfig(
            n_components=2, statistic='uniform', initial_lifetimes=(0.8, 3.0),
            n_starts=1))
        order = np.argsort(result.lifetimes)
        np.testing.assert_allclose(result.lifetimes[order], lifetimes, rtol=5e-4)
        np.testing.assert_allclose(result.backgrounds, backgrounds, atol=0.02)
        self.assertEqual(result.amplitudes.shape, (3, 2))
        self.assertEqual(result.fitted.shape, counts.shape)
        np.testing.assert_allclose(result.fitted, counts, rtol=1e-4, atol=0.02)

    def test_noisy_global_fit_with_distinct_measured_irfs(self):
        rng = np.random.default_rng(42)
        irfs = np.column_stack([
            sampled_irf(self.time, 1.2, 0.1),
            sampled_irf(self.time, 1.4, 0.15),
            sampled_irf(self.time, 1.0, 0.12)])
        mixtures = [[10000, 1000], [3000, 6000], [800, 9000]]
        counts = np.column_stack([
            rng.poisson(independent_decay(self.time, irfs[:, j],
                                          [0.7, 3.3], a, 3 + j))
            for j, a in enumerate(mixtures)])
        result = fit_global(self.time, counts, irf=irfs, config=FitConfig(
            n_components=2, statistic='poisson',
            initial_lifetimes=(0.5, 4.0), n_starts=1))
        np.testing.assert_allclose(np.sort(result.lifetimes), [0.7, 3.3],
                                   rtol=0.025)
        self.assertTrue(np.isfinite(result.lifetime_stderr).all())
        self.assertTrue(result.diagnostics['converged'])

    def test_fit_window_masks_after_full_reconvolution(self):
        # Start on the broad IRF's falling edge, where cropping before
        # convolution would discard a substantial excitation contribution.
        irf = sampled_irf(self.time, center=2.0, width=0.3)
        counts = independent_decay(self.time, irf, [2.2], [4000], 2.5)
        result = fit_decay(self.time, counts, irf=irf, config=FitConfig(
            statistic='uniform', fit_start=2.2, fit_end=12.0,
            initial_lifetimes=(1.8,), n_starts=1,
            fit_background=False, background=2.5))
        self.assertAlmostEqual(result.lifetimes[0], 2.2, delta=2e-4)
        selected = (self.time >= 2.2) & (self.time <= 12.0)
        np.testing.assert_array_equal(result.fit_mask, selected)
        np.testing.assert_allclose(result.fitted[selected, 0], counts[selected],
                                   rtol=1e-5, atol=1e-3)
        self.assertTrue(np.isnan(result.residuals[~selected]).all())

    def test_partial_fixed_lifetimes_stay_fixed(self):
        counts = independent_decay(self.time, self.irf, [0.7, 3.1],
                                   [5000, 5000], 3.0)
        result = fit_decay(self.time, counts, irf=self.irf, config=FitConfig(
            n_components=2, statistic='uniform',
            initial_lifetimes=(0.7, 2.7), fixed_lifetimes=(0.7, None),
            n_starts=1))
        np.testing.assert_allclose(np.sort(result.lifetimes), [0.7, 3.1], rtol=5e-4)
        self.assertTrue(np.any(result.lifetimes == 0.7))

    def test_positive_irf_shift_delays_response(self):
        shift = 0.175
        shifted_irf = np.interp(self.time - shift, self.time, self.irf,
                                left=0.0, right=0.0)
        counts = independent_decay(self.time, shifted_irf, [2.0], [6000], 4.0)
        result = fit_decay(self.time, counts, irf=self.irf, config=FitConfig(
            statistic='uniform', initial_lifetimes=(1.7,),
            fit_shift=True, shift=0.0, shift_bounds=(-0.5, 0.5), n_starts=1))
        self.assertAlmostEqual(result.shift, shift, delta=3e-4)
        self.assertAlmostEqual(result.lifetimes[0], 2.0, delta=3e-4)

    def test_tail_amplitude_refers_to_first_fitted_channel(self):
        start, tau, amplitude = 2.0, 2.6, 3000.0
        counts = amplitude * np.exp(-(self.time - start) / tau) + 6.0
        result = fit_decay(self.time, counts, config=FitConfig(
            method='tail', statistic='uniform', fit_start=start,
            initial_lifetimes=(2.0,), n_starts=1))
        self.assertAlmostEqual(result.lifetimes[0], tau, delta=3e-4)
        self.assertAlmostEqual(result.amplitudes[0, 0], amplitude, delta=0.03)
        self.assertAlmostEqual(result.backgrounds[0], 6.0, delta=0.01)

    def test_irf_background_subtraction_is_explicit(self):
        counts = independent_decay(self.time, self.irf, [1.9], [4500], 5.0)
        observed_irf = self.irf * 10000 + 2.0
        result = fit_decay(self.time, counts, irf=observed_irf, config=FitConfig(
            statistic='uniform', initial_lifetimes=(1.5,),
            irf_background=2.0, n_starts=1))
        self.assertAlmostEqual(result.lifetimes[0], 1.9, delta=3e-4)

    def test_periodic_prediction_matches_independent_circular_sum(self):
        time = np.arange(80) * 0.05
        dt, tau, amplitude = 0.05, 2.3, 600.0
        period = len(time) * dt
        irf = sampled_irf(time, center=0.7, width=0.1)
        decay = (amplitude * np.exp(-time / tau)
                 * tau / dt * (-np.expm1(-dt / tau))
                 / (-np.expm1(-period / tau)))
        expected = sum(irf[j] * np.roll(decay, j) for j in range(len(time)))
        got = predict(time, [tau], [amplitude], irf=irf,
                      repetition_period=period)
        np.testing.assert_allclose(got[:, 0], expected, rtol=1e-9, atol=1e-8)

    def test_invalid_inputs_are_rejected(self):
        counts = independent_decay(self.time, self.irf, [2.0], [2000], 2)
        invalid_counts = counts.copy()
        invalid_counts[5] = -1.0
        with self.assertRaises(ValueError):
            fit_decay(self.time, invalid_counts, irf=self.irf)
        irregular = self.time.copy()
        irregular[20] += 0.01
        with self.assertRaises(ValueError):
            fit_decay(irregular, counts, irf=self.irf)
        with self.assertRaises(ValueError):
            fit_decay(self.time, counts, irf=self.irf, config=FitConfig(
                repetition_period=15.0))
        with self.assertRaises(ValueError):
            fit_decay(self.time, counts, config=FitConfig(method='reconvolution'))


if __name__ == '__main__':
    unittest.main()
