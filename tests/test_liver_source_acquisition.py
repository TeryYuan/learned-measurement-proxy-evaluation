"""Bounded tests for new acquisition boundary; no science or network fixtures."""
import hashlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
import zipfile
from contextlib import redirect_stderr

spec = importlib.util.spec_from_file_location('fetch_liver_sources', Path(__file__).resolve().parents[1]/'tools/fetch_liver_sources.py')
fetch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fetch)


def workbook(names=('version 1', 'version 2')):
    b = io.BytesIO()
    with zipfile.ZipFile(b, 'w') as z:
        z.writestr('[Content_Types].xml', '<Types/>')
        z.writestr('xl/workbook.xml', '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"><sheets>'+''.join(f'<sheet name="{n}"/>' for n in names)+'</sheets></workbook>')
        z.writestr('xl/worksheets/sheet1.xml', '<worksheet/>')
    data = b.getvalue()
    return data, {'key': 'FDA_DILIrank', 'filename': 'source.xlsx', 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest()}


class TestSourceBoundary(unittest.TestCase):
    def test_html_rejected(self):
        with self.assertRaisesRegex(fetch.SourceError, 'html_instead'):
            fetch.validate_workbook(b'<html>error</html>', {})

    def test_latest_cannot_replace_historical_version(self):
        data, record = workbook()
        with self.assertRaisesRegex(fetch.SourceError, 'historical_version_mismatch'):
            fetch.validate_workbook(data+b'changed', record)

    def test_fda_requires_independent_v1_sheet(self):
        data, record = workbook(('version 2',))
        with self.assertRaisesRegex(fetch.SourceError, 'required_sheets_missing'):
            fetch.validate_workbook(data, record)

    def test_valid_exact_workbook(self):
        data, record = workbook()
        fetch.validate_workbook(data, record)

    def test_local_missing_is_specific(self):
        _, record = workbook()
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(fetch.SourceError, 'required_source_missing'):
                fetch.acquire(record, Path(tmp))

    def test_existing_output_never_overwritten(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stderr(io.StringIO()) as err:
            sentinel = Path(tmp)/'keep';sentinel.write_text('original')
            self.assertEqual(fetch.main(['--output-dir', tmp]), 2)
            self.assertEqual(sentinel.read_text(), 'original')
            self.assertIn('output_directory_exists', err.getvalue())
            self.assertNotIn(tmp, err.getvalue())

    def test_missing_local_input_reports_no_private_path_or_traceback(self):
        with tempfile.TemporaryDirectory() as tmp, redirect_stderr(io.StringIO()) as err:
            output = Path(tmp)/'output'
            self.assertEqual(fetch.main(['--input-dir', tmp, '--output-dir', str(output), '--sources', 'SI004']), 2)
            self.assertIn('SI004: required_source_missing', err.getvalue())
            self.assertNotIn(tmp, err.getvalue());self.assertNotIn('Traceback', err.getvalue())
            self.assertFalse(output.exists())


if __name__ == '__main__':
    unittest.main()
