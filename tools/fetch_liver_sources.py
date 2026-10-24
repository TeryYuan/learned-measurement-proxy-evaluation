#!/usr/bin/env python3
"""Acquire exact historical source workbooks; never regenerate saved scores.

Only the three explicitly listed official files are fetched. With --input-dir,
validate user-supplied copies instead of accessing the network. Source terms
remain applicable; see DATA_ACCESS_AND_REPLAY.md and THIRD_PARTY_NOTICES.md.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
from pathlib import Path
import sys
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'sources/liver_source_downloads.json'


class SourceError(Exception):
    """Public errors contain fixed codes and source keys only."""


def validate_workbook(data: bytes, record: dict) -> None:
    if data.lstrip().lower().startswith((b'<!doctype html', b'<html')):
        raise SourceError('html_instead_of_workbook')
    if not data.startswith(b'PK\x03\x04'):
        raise SourceError('not_xlsx')
    if len(data) != record['bytes'] or hashlib.sha256(data).hexdigest() != record['sha256']:
        raise SourceError('historical_version_mismatch')
    try:
        with zipfile.ZipFile(io.BytesIO(data)) as z:
            if z.testzip() is not None:
                raise SourceError('xlsx_crc_failure')
            for name in ('[Content_Types].xml', 'xl/workbook.xml'):
                ET.fromstring(z.read(name))
            wb = ET.fromstring(z.read('xl/workbook.xml'))
            ns = {'s': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
            sheets = [s.attrib['name'] for s in wb.findall('s:sheets/s:sheet', ns)]
            if not sheets or (record['key'] == 'FDA_DILIrank' and not {'version 1', 'version 2'} <= set(sheets)):
                raise SourceError('required_sheets_missing')
            if not any(n.startswith('xl/worksheets/sheet') and n.endswith('.xml') for n in z.namelist()):
                raise SourceError('worksheet_content_missing')
    except (zipfile.BadZipFile, KeyError, ET.ParseError) as error:
        raise SourceError('invalid_xlsx_structure') from error


def acquire(record: dict, input_dir: Path | None) -> bytes:
    if input_dir is not None:
        source = input_dir / record['filename']
        if not source.is_file():
            raise SourceError('required_source_missing')
        if source.is_symlink():
            raise SourceError('source_link_rejected')
        with source.open('rb') as f:
            data = f.read(record['bytes'] + 1)
    else:
        try:
            request = urllib.request.Request(record['url'], headers={'User-Agent': 'Measurement-proxy-source-audit/1.0'})
            with urllib.request.urlopen(request, timeout=60) as response:
                if response.status != 200 or not response.geturl().startswith('https://'):
                    raise SourceError('source_http_failure')
                data = response.read(record['bytes'] + 1)
        except (urllib.error.URLError, TimeoutError, OSError) as error:
            raise SourceError('source_unavailable_use_exact_local_copy') from error
    validate_workbook(data, record)
    return data


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, help='Optional directory containing legally obtained exact historical workbooks; disables network.')
    parser.add_argument('--output-dir', type=Path, required=True, help='New directory; existing directories are rejected.')
    parser.add_argument('--sources', nargs='+', choices=['SI002', 'SI004', 'FDA_DILIrank'], default=['SI002', 'SI004', 'FDA_DILIrank'])
    args = parser.parse_args(argv)
    if args.output_dir.exists() or args.output_dir.is_symlink():
        print('ERROR: output_directory_exists', file=sys.stderr)
        return 2
    records = json.loads(MANIFEST.read_text(encoding='utf-8'))['sources']
    data = {}
    for record in records:
        if record['key'] not in args.sources:
            continue
        try:
            data[record['filename']] = acquire(record, args.input_dir)
        except (SourceError, OSError) as error:
            code = str(error) if isinstance(error, SourceError) else 'source_read_failure'
            print(f"ERROR: {record['key']}: {code}; see sources/liver_source_downloads.json for the exact filename and version.", file=sys.stderr)
            return 2
    try:
        args.output_dir.mkdir(parents=True, exist_ok=False)
        for filename, content in data.items():
            with (args.output_dir / filename).open('xb') as f:
                f.write(content)
        report = {'status': 'PASS', 'source_keys': args.sources, 'historical_versions_exact': True,
                  'origin': 'user_supplied_exact_copies' if args.input_dir else 'official_downloads',
                  'saved_predictions_created': False, 'sampling_arrays_created': False}
        (args.output_dir / 'acquisition_report.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    except OSError:
        print('ERROR: output_write_failure; inspect the selected output directory before retrying.', file=sys.stderr)
        return 2
    print(f'PASS: {len(data)} exact historical workbooks validated. No predictions, labels or sampling arrays were inferred.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
