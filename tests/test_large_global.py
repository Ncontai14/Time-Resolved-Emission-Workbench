"""Sparse global optimization equivalence and tail-domain regressions."""
from __future__ import annotations

import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

from tcspc import FitConfig, fit_decay, fit_global, predict
from tcspc import core


class SparseGlobalChecks(unittest.TestCase):
    @staticmethod
    def sample(time, irf, tau, amplitude, background):
        # Independent unaveraged discrete model: fitted amplitudes absorb the
        # known bin-average factor, while lifetimes/likelihood stay identical.
        lag = np.arange(len(time))*(time[1]-time[0])
        decay = sum(a*np.exp(-lag/t) for a, t in zip(amplitude, tau))
        return np.convolve(irf/irf.sum(), decay)[:len(time)] + background

    def test_sparse_and_dense_likelihood_and_uncertainty_agree(self):
        rng = np.random.default_rng(360)
        time = np.arange(240)*0.05
        irf = np.exp(-0.5*((time-1.1)/0.1)**2)
        amplitudes = [[1000+300*j, 4000-350*j] for j in range(8)]
        counts = np.column_stack([
            rng.poisson(self.sample(time, irf, [0.65, 2.8], a, 2+j/3))
            for j, a in enumerate(amplitudes)])
        config = FitConfig(n_components=2, initial_lifetimes=(0.8, 3.2),
                           n_starts=1)
        sparse = fit_global(time, counts, irf=irf, config=config)
        with patch.object(core, '_SPARSE_TRACE_THRESHOLD', 10**9), \
             patch.object(core, '_DENSE_JACOBIAN_ENTRY_LIMIT', 10**15):
            dense = fit_global(time, counts, irf=irf, config=config)
        self.assertEqual(sparse.diagnostics['jacobian_storage'], 'sparse CSR')
        self.assertEqual(dense.diagnostics['jacobian_storage'], 'dense')
        np.testing.assert_allclose(sparse.lifetimes, dense.lifetimes, rtol=2e-5)
        np.testing.assert_allclose(sparse.lifetime_stderr, dense.lifetime_stderr,
                                   rtol=2e-4)
        np.testing.assert_allclose(sparse.amplitudes, dense.amplitudes, rtol=2e-4)
        self.assertAlmostEqual(sparse.diagnostics['objective'],
                               dense.diagnostics['objective'], delta=1e-5)

    def test_sparse_poisson_zero_bins_recover_shared_lifetime(self):
        rng = np.random.default_rng(2026)
        time = np.arange(300)*0.06
        irf = np.exp(-0.5*((time-1.0)/0.1)**2)
        counts = np.column_stack([
            rng.poisson(self.sample(time, irf, [1.7], [80+10*j], 0.01))
            for j in range(8)])
        self.assertTrue(np.any(counts == 0))
        result = fit_global(time, counts, irf=irf, config=FitConfig(
            initial_lifetimes=(2.0,), n_starts=1))
        self.assertAlmostEqual(result.lifetimes[0], 1.7, delta=0.08)
        self.assertTrue(np.isfinite(result.fitted).all())
        self.assertTrue(np.isfinite(result.residuals).all())
        self.assertTrue(result.diagnostics['converged'])

    def test_tail_undefined_rise_is_nan_and_exports_as_blank(self):
        time = np.arange(500)*0.1
        start, tau = 35.0, 0.05
        counts = np.full(len(time), 2.0)
        included = time >= start
        counts[included] += 3000*np.exp(-(time[included]-start)/tau)
        with np.errstate(over='raise', invalid='raise'):
            prediction = predict(time, [tau], [3000], backgrounds=2,
                                 method='tail', fit_start=start)
            result = fit_decay(time, counts, config=FitConfig(
                method='tail', fit_start=start, statistic='uniform',
                initial_lifetimes=(tau,), fixed_lifetimes=(tau,),
                n_starts=1))
        self.assertTrue(np.isnan(prediction[~included]).all())
        self.assertTrue(np.isnan(result.fitted[~included]).all())
        self.assertTrue(np.isfinite(result.fitted[included]).all())
        with tempfile.TemporaryDirectory() as directory:
            result.export(directory)
            with (Path(directory)/'curves.csv').open(newline='') as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows[0]['fit_1'], '')
            self.assertNotEqual(rows[np.flatnonzero(included)[0]]['fit_1'], '')


if __name__ == '__main__':
    unittest.main()
