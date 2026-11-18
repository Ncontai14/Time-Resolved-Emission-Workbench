import unittest
import numpy as np
from tcspc import FitConfig, fit_decay
from tcspc.io import load_decay, load_tres
from tcspc.time_axis import channel_grid


class TimeRoundingChecks(unittest.TestCase):
    def test_rounded_export_loads_without_changing_counts(self):
        t = (np.arange(2001)*2.5e-11).astype(np.float32).astype(float)*1e9
        c = np.arange(len(t)) % 17
        content = '\n'.join(f'{x:.16g}\t{y}' for x,y in zip(t,c)).encode()
        loaded = load_decay(content)
        np.testing.assert_array_equal(loaded.counts, c)
        np.testing.assert_allclose(loaded.metadata['original_time_ns'], t, atol=1e-13)
        self.assertTrue(loaded.metadata['rounded_channel_labels'])
        self.assertFalse(loaded.metadata['photon_counts_resampled'])
        tres = load_tres(content, wavelengths=[500])
        np.testing.assert_array_equal(tres.counts[:,0], c)

    def test_direct_fit_accepts_rounded_channel_labels(self):
        ideal = np.arange(400)*0.025
        rounded = (ideal*1e-9).astype(np.float32).astype(float)*1e9
        before = rounded.copy()
        counts = 1000*np.exp(-ideal/1.7)+2
        result = fit_decay(rounded, counts, config=FitConfig(method='tail', statistic='uniform', n_starts=1))
        self.assertAlmostEqual(result.lifetimes[0], 1.7, places=5)
        np.testing.assert_array_equal(rounded, before)
        np.testing.assert_array_equal(result.observed[:,0], counts)
        self.assertTrue(result.diagnostics['time_axis']['rounded_channel_labels'])

    def test_missing_or_genuinely_irregular_channels_are_identified(self):
        ideal = np.arange(200)*0.025
        missing = np.delete(ideal,100)
        irregular = ideal.copy(); irregular[100] += .001
        for bad in [missing, irregular]:
            with self.assertRaisesRegex(ValueError, 'missing channels'):
                channel_grid(bad)

    def test_calibrated_grid_retains_labels_and_offset(self):
        ideal = -2+np.arange(200)*.025
        got, audit = channel_grid(ideal)
        np.testing.assert_array_equal(got, ideal)
        self.assertFalse(audit['rounded_channel_labels'])


if __name__ == '__main__':
    unittest.main()
