"""Fast AC0 tests: standard library only; do not require Hail or Docker.

Run from project root: python -m unittest discover -s tests -p 'test_ac0_header.py'
"""

import ast
import csv
import gzip
import json
import logging
import os
import re
import tempfile
import types
import unittest
from pathlib import Path

# Load the real AC0 definitions embedded in M4 without importing Hail/Spark.
# The full M4 module intentionally needs Hail, but header/report checks do not.
M4_SOURCE = Path(__file__).resolve().parents[1] / 'AF_hail_pipeline' / 'M4_af_annotation.py'
parsed = ast.parse(M4_SOURCE.read_text())
ac0_names = {
    'list_input_vcfs', 'inspect_vcf_header', 'inspect_input_formats',
    'new_report', '_format_schema_problems', '_count_incomplete_genotypes',
    '_flag_when_ac_zero', 'maybe_flag_ac0', 'write_ac0_reports',
    'append_to_workflow_csv',
}
constant_names = {
    'REQUIRED_FORMAT_FIELDS', 'GQ_MIN', 'DIPLOID_DP_MIN',
    'HAPLOID_DP_MIN', 'HET_AB_MIN', 'FORMAT_DECLARATION',
}
body = []
for node in parsed.body:
    if isinstance(node, ast.FunctionDef) and node.name in ac0_names:
        body.append(node)
    elif isinstance(node, ast.Assign) and any(
        isinstance(target, ast.Name) and target.id in constant_names
        for target in node.targets
    ):
        body.append(node)
namespace = {
    'csv': csv, 'gzip': gzip, 'json': json, 'logging': logging,
    'os': os, 're': re, 'Path': Path,
}
exec(compile(ast.Module(body=body, type_ignores=[]), str(M4_SOURCE), 'exec'), namespace)
ac0 = types.SimpleNamespace(**namespace)


class AC0HeaderTests(unittest.TestCase):
    def write_vcf(self, folder, name, fmt):
        header = ['##fileformat=VCFv4.2']
        for field in fmt:
            header.append(f'##FORMAT=<ID={field},Number=1,Type=String,Description="{field}">')
        header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1')
        data = '\n'.join(header) + '\n'
        path = folder / name
        if name.endswith(('.gz', '.bgz')):
            with gzip.open(path, 'wt') as out:
                out.write(data)
        else:
            path.write_text(data)
        return path

    def test_full_and_missing_headers_per_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.write_vcf(root, 'complete.vcf', ['GT', 'AD', 'DP', 'GQ'])
            self.write_vcf(root, 'incomplete.vcf.bgz', ['GT', 'GQ'])
            (root / 'complete.vcf.tbi').write_text('index')
            audited = ac0.inspect_input_formats(root)
            self.assertEqual(len(audited), 2)
            by_name = {Path(row['file']).name: row for row in audited}
            self.assertEqual(by_name['complete.vcf']['missing_fields'], [])
            self.assertEqual(by_name['incomplete.vcf.bgz']['missing_fields'], ['AD', 'DP'])

    def test_no_sample_columns_counts_as_unavailable(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            path = self.write_vcf(root, 'sites.vcf', ['GT', 'AD', 'DP', 'GQ'])
            content = path.read_text().replace('\tFORMAT\tS1', '')
            path.write_text(content)
            report = ac0.inspect_input_formats(root)
            self.assertIn('SAMPLE_COLUMNS', report[0]['missing_fields'])

    def test_enabled_missing_fields_is_nonfatal_even_without_hail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.write_vcf(root, 'incomplete.vcf.gz', ['GT'])
            original_matrix = object()
            returned_matrix, report = ac0.maybe_flag_ac0(original_matrix, True, root)
            self.assertIs(returned_matrix, original_matrix)
            self.assertEqual(report['status'], 'SKIPPED_MISSING_FORMAT_FIELDS')
            self.assertEqual(report['missing_fields'], ['AD', 'DP', 'GQ'])
            ac0.write_ac0_reports(report, str(root / 'ac0_status'))
            parsed = json.loads((root / 'ac0_status.json').read_text())
            self.assertEqual(parsed['status'], 'SKIPPED_MISSING_FORMAT_FIELDS')
            self.assertIn('GQ', (root / 'ac0_status.tsv').read_text())

            rows = []
            ac0.append_to_workflow_csv(report, lambda data: rows.extend(data))
            self.assertIn(['Status', 'SKIPPED_MISSING_FORMAT_FIELDS'], rows)
            self.assertIn(['Missing FORMAT', 'AD, DP, GQ'], rows)

    def test_disabled_does_not_access_inputs_or_change_matrix(self):
        original_matrix = object()
        out, report = ac0.maybe_flag_ac0(original_matrix, False, '/no/such/path')
        self.assertIs(out, original_matrix)
        self.assertEqual(report['status'], 'DISABLED')

    def test_bad_vcf_is_nonfatal_to_ac0(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / 'broken.vcf').write_text('this is not VCF\n')
            out, report = ac0.maybe_flag_ac0(object(), True, root)
            self.assertEqual(report['status'], 'SKIPPED_INPUT_AUDIT')

    def test_missing_input_directory_is_nonfatal_to_ac0(self):
        original_matrix = object()
        out, report = ac0.maybe_flag_ac0(original_matrix, True, '/no/such/path')
        self.assertIs(out, original_matrix)
        self.assertEqual(report['status'], 'SKIPPED_AC0_ERROR')


if __name__ == '__main__':
    unittest.main()
