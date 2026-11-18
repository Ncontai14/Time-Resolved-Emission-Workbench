import io
import unittest

import numpy as np

from tcspc.io import DataFormatError, load_decay, load_tres


class TextImportTests(unittest.TestCase):
    def test_header_units_are_explicit_and_upload_is_not_consumed(self):
        upload = io.BytesIO(b"Instrument export\ntime_ns\tcounts\n0\t3\n25\t7\n50\t4\n")
        upload.seek(2)
        result = load_decay(upload, time_unit="ps")
        np.testing.assert_allclose(result.time_ns, [0, 0.025, 0.05])
        np.testing.assert_array_equal(result.counts, [3, 7, 4])
        self.assertEqual(upload.tell(), 2)

    def test_csv_and_semicolon_and_whitespace(self):
        for separator in [",", ";", "\t", " "]:
            with self.subTest(separator=separator):
                data = "\n".join(separator.join(map(str, row)) for row in [(0, 2), (1, 3), (2, 1)])
                self.assertEqual(load_decay(data.encode()).counts.sum(), 6)

    def test_count_only_requires_calibration(self):
        data = b"2\n3\n1\n"
        with self.assertRaises(DataFormatError):
            load_decay(data, time_column=None, count_column=0)
        result = load_decay(data, time_column=None, count_column=0, bin_width=0.1)
        np.testing.assert_allclose(result.time_ns, [0, 0.1, 0.2])

    def test_tres_header_and_explicit_wavelengths(self):
        data = b"time_ns,500,520\n0,2,1\n0.1,3,2\n0.2,1,1\n"
        tres = load_tres(data)
        self.assertEqual(tres.counts.shape, (3, 2))
        np.testing.assert_array_equal(tres.wavelengths, [500, 520])
        result = load_tres(b"0 2 1\n1 3 2\n2 1 1", wavelengths=[510, 530])
        np.testing.assert_array_equal(result.wavelengths, [510, 530])

    def test_invalid_rows_and_counts_are_rejected(self):
        cases = [
            b"0 1\n1 2\n2 -1",  # subtracted counts
            b"0 1\n1 2\n1 3",  # duplicate bin
            b"0 1\n1 2\n3 3",  # nonuniform bins
            b"0 1\n1 nan\n2 3",  # missing data must not be skipped
            b"0 1\n1 \n2 3",  # malformed row
            b"0 1\n1 bad\n2 3",
            b"0 0\n1 0\n2 0",
            b"0,0\t1\n1,0\t2\n2,0\t3",  # decimal comma
        ]
        for data in cases:
            with self.subTest(data=data):
                with self.assertRaises(DataFormatError):
                    load_decay(data)

    def test_native_binary_has_actionable_error(self):
        with self.assertRaisesRegex(DataFormatError, "ASCII"):
            load_decay(b"Native Container\x00\x00")

    def test_missing_wavelengths_are_not_fabricated(self):
        with self.assertRaisesRegex(DataFormatError, "wavelength"):
            load_tres(b"0 2 1\n1 3 2\n2 1 1")

    def test_nonfinite_and_duplicate_wavelengths(self):
        for wavelengths in ([500, 500], [500, np.nan], [500], [-1, 500]):
            with self.subTest(wavelengths=wavelengths):
                with self.assertRaises(DataFormatError):
                    load_tres(b"0 2 1\n1 3 2\n2 1 1", wavelengths=wavelengths)


if __name__ == "__main__":
    unittest.main()
