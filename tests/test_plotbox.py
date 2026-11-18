from io import BytesIO
from pathlib import Path
import unittest

import numpy as np

from tcspc.plotbox import is_plotbox_export, load_plotbox
from tcspc.io import DataFormatError


def export(rows, unit="s"):
    return (f"PlotBox Export\nTime[{unit}]X\tcrv[0] []\tcrv[1] []\n" +
            "\n".join("\t".join(str(v) for v in row) for row in rows)).encode()


class PlotBoxImportTests(unittest.TestCase):
    def test_seconds_are_explicit_and_roles_not_invented(self):
        data = export([(i*0.025e-9, i+2, 10-i) for i in range(8)])
        bundle = load_plotbox(data)
        np.testing.assert_allclose(bundle.time_ns, np.arange(8)*0.025)
        self.assertEqual(bundle.source_time_unit, "s")
        self.assertEqual(bundle.metadata["curve_roles"], "unassigned")
        self.assertIsNone(bundle.metadata["wavelengths_nm"])
        with self.assertRaises(DataFormatError):
            load_plotbox(data, time_unit="ns")

    def test_missing_endpoints_remain_missing_and_are_reported(self):
        rows = [(i*0.025e-9, i+2, 10-i) for i in range(7)]
        rows[-1] = (rows[-1][0], "-", 4)
        bundle = load_plotbox(export(rows))
        self.assertTrue(np.isnan(bundle.counts[-1, 0]))
        time, counts = bundle.select_curves([1, 0])
        self.assertEqual(counts.shape, (6, 2))
        self.assertEqual(counts[0, 0], 10)
        audit = bundle.selection_summary([0, 1])
        self.assertEqual(audit["trailing_endpoint_rows_trimmed"], 1)
        self.assertEqual(audit["counts_omitted_at_endpoints"], 4)
        with self.assertRaises(DataFormatError):
            bundle.select_curves([0], missing_policy="raise")

    def test_unknown_endpoint_timestamps_preserve_raw_counts(self):
        rows = [(i*0.025e-9, 1, 2) for i in range(6)] + [("-", "-", 13)]
        bundle = load_plotbox(export(rows))
        self.assertEqual(len(bundle.time_ns), 6)
        self.assertEqual(bundle.raw_counts[-1, 1], 13)
        self.assertEqual(bundle.metadata["dropped_timestamp_counts_sum"], 13)
        self.assertEqual(bundle.metadata["dropped_timestamp_data_rows"], [6])

    def test_internal_missing_counts_or_timestamps_are_rejected(self):
        rows = [(i*0.025e-9, 1, 2) for i in range(8)]
        rows[3] = (rows[3][0], "-", 2)
        bundle = load_plotbox(export(rows))
        with self.assertRaises(DataFormatError):
            bundle.select_curves([0])
        self.assertEqual(len(bundle.select_curves([1])[0]), 8)
        rows[3] = ("-", 1, 2)
        with self.assertRaises(DataFormatError):
            load_plotbox(export(rows))

    def test_float32_coordinates_canonicalize_but_real_gaps_do_not(self):
        original = (np.arange(2001)*0.025e-9).astype(np.float32).astype(float)
        rows = [(t, 1, 2) for t in original]
        bundle = load_plotbox(export(rows))
        np.testing.assert_allclose(np.diff(bundle.time_ns), bundle.metadata["bin_width_ns"], atol=1e-13)
        np.testing.assert_array_equal(bundle.original_time_ns, original*1e9)
        self.assertGreater(bundle.metadata["max_timestamp_adjustment_ns"], 0)
        self.assertLess(bundle.metadata["max_timestamp_adjustment_bins"], 1e-4)
        del rows[100]
        with self.assertRaises(DataFormatError):
            load_plotbox(export(rows))

    def test_malformed_negative_and_nonfinite_cells_rejected(self):
        valid = [(i*0.025e-9, 1, 2) for i in range(6)]
        for bad in [-1, "nan", "Inf", "nonsense"]:
            rows = valid.copy()
            rows[2] = (rows[2][0], bad, 2)
            with self.assertRaises(DataFormatError):
                load_plotbox(export(rows))
        with self.assertRaises(DataFormatError):
            load_plotbox(b"PlotBox Export\nTime[s]X\ta\tb\n0\t1\n")

    def test_upload_position_preserved_and_detection(self):
        stream = BytesIO(export([(i, 1, 2) for i in range(6)], unit="ns"))
        stream.seek(9)
        self.assertTrue(is_plotbox_export(stream))
        self.assertEqual(stream.tell(), 9)
        bundle = load_plotbox(stream)
        self.assertEqual(stream.tell(), 9)
        with self.assertRaises(DataFormatError):
            bundle.select_curves([0, 0])
        with self.assertRaises(DataFormatError):
            bundle.select_curves([2])

    def test_supplied_exports_import_and_select_all_curves(self):
        root = Path(__file__).resolve().parents[1]/"examples"/"curve_exports"
        if not root.exists():
            self.fail("Bundled instrument exports are missing from the package")
        expected = {
            "anisodec.dat": (1501, 2, 1501),
            "anisodec_IRF.dat": (1601, 3, 1601),
            "anisodec_IRF_magic.dat": (1601, 4, 1600),
            "decay.dat": (940, 1, 940),
            "decay_IRF.dat": (1601, 2, 1601),
            "decaytimecourse_IRF.dat": (1599, 11, 1599),
            "TRES_IRF.dat": (2001, 52, 1999),
        }
        for name, (n, m, common) in expected.items():
            with self.subTest(name=name):
                b = load_plotbox(root/name)
                self.assertEqual(b.counts.shape, (n, m))
                time, counts = b.select_curves(range(m))
                self.assertEqual(len(time), common)
                self.assertTrue(np.isfinite(counts).all())


if __name__ == "__main__":
    unittest.main()
