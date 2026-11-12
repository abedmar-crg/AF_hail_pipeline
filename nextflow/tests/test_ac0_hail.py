"""Hail integration smoke tests for optional AC0 FILTER flagging.

Run this inside the rebuilt pipeline Docker image (see README.md), because
Hail/Spark is not required for the lightweight unit tests.
"""

import sys
import tempfile
from pathlib import Path

try:
    import hail as hl
except ModuleNotFoundError:
    # Allow unittest discovery on developer machines without Hail installed.
    hl = None

source_folder = Path(__file__).resolve().parents[1] / 'AF_hail_pipeline'
# In the Docker smoke test this file is mounted at /tests, while sources
# are baked into /app. In a checkout, they are at ../AF_hail_pipeline.
sys.path.insert(0, str(source_folder if source_folder.exists() else Path('/app')))
if hl is not None:
    import M4_af_annotation as ac0
else:
    ac0 = None


def write_vcf(path, missing_value=False):
    rows = [
        '##fileformat=VCFv4.2',
        '##contig=<ID=1,length=249250621>',
        '##FILTER=<ID=LowQual,Description="Low quality">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
        '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="GQ">',
        '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">',
        '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="AD">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4',
        # All ALT-containing genotypes fail, while the passing 0/0 contributes AN but no ALT AC.
        '1\t10001\t.\tA\tC\t90\tPASS\t.\tGT:GQ:DP:AD\t0/1:10:20:10,10\t0/1:40:5:2,3\t1/1:3:20:0,20\t0/0:60:20:20,0',
        # Passing ALT must NOT be flagged, original LowQual filter must be retained.
        '1\t10002\t.\tC\tT\t90\tLowQual\t.\tGT:GQ:DP:AD\t0/1:60:20:10,10\t0/0:60:20:20,0\t0/0:60:20:20,0\t0/0:60:20:20,0',
        # 1|0 is a phased het-ref and must test AD for the actual ALT allele.
        '1\t10003\t.\tG\tT\t90\tPASS\t.\tGT:GQ:DP:AD\t1|0:60:20:2,18\t0/0:60:20:20,0\t0/0:60:20:20,0\t0/0:60:20:20,0',
        # Haploid DP >= 5 qualifies.
        '1\t10004\t.\tG\tA\t90\tPASS\t.\tGT:GQ:DP:AD\t1:60:5:0,5\t0/0:60:20:20,0\t0/0:60:20:20,0\t0/0:60:20:20,0',
    ]
    if missing_value:
        rows[-1] = rows[-1].replace('1:60:5:0,5', '1:.:5:0,5')
    path.write_text('\n'.join(rows) + '\n')


def main():
    if hl is None:
        raise SystemExit('Install Hail or run this test in the Docker image.')
    hl.init(quiet=True, log='/tmp/ac0_hail_smoke.log')
    try:
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            file = folder / 'demo.vcf'
            write_vcf(file)
            mt = hl.import_vcf(str(file), reference_genome='GRCh37', min_partitions=1)
            annotated, report = ac0.maybe_flag_ac0(mt, True, str(folder))
            assert report['status'] == 'APPLIED', report
            assert report['variants_total'] == 4, report
            assert report['variants_evaluated'] == 4, report
            assert report['variants_skipped_missing_values'] == 0, report
            assert report['variants_flagged_ac0'] == 1, report
            before = mt.count()
            assert annotated.count() == before, (before, annotated.count())
            assert annotated.entry.dtype == mt.entry.dtype, 'AC0 must not mask/change FORMAT entries'

            def filt_matrix(table):
                return {
                    row.locus.position: set(row.filters)
                    for row in table.rows().select('filters').collect()
                }
            assert filt_matrix(annotated) == {
                10001: {'AC0'},
                10002: {'LowQual'},
                10003: set(),
                10004: set(),
            }, filt_matrix(annotated)

            write_vcf(file, missing_value=True)
            mt_bad = hl.import_vcf(str(file), reference_genome='GRCh37', min_partitions=1)
            partial, bad_report = ac0.maybe_flag_ac0(mt_bad, True, str(folder))
            assert bad_report['status'] == 'PARTIALLY_APPLIED_MISSING_FORMAT_VALUES', bad_report
            assert bad_report['missing_values']['missing_gq'] == 1, bad_report
            assert bad_report['variants_total'] == 4, bad_report
            assert bad_report['variants_evaluated'] == 3, bad_report
            assert bad_report['variants_skipped_missing_values'] == 1, bad_report
            assert bad_report['variants_flagged_ac0'] == 1, bad_report
            assert partial.count() == mt_bad.count(), 'No variant or sample should be removed'
            assert partial.entry.dtype == mt_bad.entry.dtype, 'FORMAT values should be preserved'
            assert filt_matrix(partial) == {
                10001: {'AC0'},
                10002: {'LowQual'},
                10003: set(),
                10004: set(),  # this row has missing GQ; leave FILTER unchanged
            }, filt_matrix(partial)
    finally:
        hl.stop()
    print('Hail AC0 integration checks PASSED')


if __name__ == '__main__':
    main()
