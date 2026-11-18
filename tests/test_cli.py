"""Public command names and metadata remain generic and portable."""
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
import json
from pathlib import Path
import tempfile
import unittest

from tcspc.cli import build_parser, main


ROOT = Path(__file__).resolve().parents[1]


class CommandLineChecks(unittest.TestCase):
    def test_public_command_names_are_generic(self):
        parser = build_parser()
        commands = next(action for action in parser._actions if action.dest == 'command')
        self.assertEqual(set(commands.choices),
                         {'fit', 'global', 'inspect', 'curve-fit', 'curve-global'})
        for name in commands.choices:
            self.assertIn(name, parser.format_help())

    def test_inspect_preserves_generic_format_and_source_curve_labels(self):
        with tempfile.TemporaryDirectory() as folder:
            target = Path(folder)/'inventory.json'
            with redirect_stdout(StringIO()):
                code = main(['inspect', '--input', str(ROOT/'examples/curve_exports/decay_IRF.dat'),
                             '--out', str(target)])
            self.assertEqual(code, 0)
            result = json.loads(target.read_text())
            self.assertEqual(result['import_metadata']['format'], 'PlotBox ASCII')
            self.assertEqual(len(result['curves']), 2)
            self.assertIn('crv[0]', result['curves'][0]['label'])

    def test_selected_curve_fit_exports_mapping_and_import_audit(self):
        with tempfile.TemporaryDirectory() as folder:
            with redirect_stdout(StringIO()), redirect_stderr(StringIO()):
                code = main(['curve-fit', '--input', str(ROOT/'examples/curve_exports/decay_IRF.dat'),
                             '--decay-curve', '1', '--irf-curve', '0', '--irf-background', '1.4',
                             '--irf-start', '3.8', '--irf-end', '8', '--fit-shift',
                             '--starts', '1', '--out', folder])
            self.assertEqual(code, 0)
            result = json.loads((Path(folder)/'result.json').read_text())
            self.assertAlmostEqual(result['lifetimes_ns'][0], 2.5491557, places=5)
            self.assertTrue((Path(folder)/'curve_mapping.csv').exists())
            metadata = json.loads((Path(folder)/'import_metadata.json').read_text())
            self.assertEqual(metadata['bundle']['format'], 'PlotBox ASCII')
            self.assertEqual(metadata['decay_selection']['curve_indices'], [1])


if __name__ == '__main__':
    unittest.main()
