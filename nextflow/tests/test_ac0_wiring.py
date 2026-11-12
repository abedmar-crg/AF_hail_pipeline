"""Static integration wiring checks for the delivered Nextflow/Docker package."""

import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class WiringTests(unittest.TestCase):
    def test_nextflow_and_docker_include_ac0(self):
        nf = (ROOT / 'main.nf').read_text()
        docker = (ROOT / 'AF_hail_pipeline' / 'Dockerfile').read_text()
        stage = (ROOT / 'AF_hail_pipeline' / 'nextflow_stage.py').read_text()
        m4 = (ROOT / 'AF_hail_pipeline' / 'M4_af_annotation.py').read_text()
        self.assertIn("path 'ac0_status.tsv', emit: ac0_status", nf)
        self.assertIn("path 'ac0_status.json', emit: ac0_details", nf)
        self.assertIn('ac0_status = M4_AF_ANNOTATION.out.ac0_status', nf)
        self.assertIn("path 'qc/ac0'", nf)
        self.assertNotIn('COPY ac0_flagging.py .', docker)
        self.assertFalse((ROOT / 'AF_hail_pipeline' / 'ac0_flagging.py').exists())
        self.assertIn('COPY M4_af_annotation.py .', docker)
        self.assertIn('m4.maybe_flag_ac0(', stage)
        self.assertIn('m4.write_ac0_reports(', stage)
        self.assertIn('m4.append_to_workflow_csv(', stage)
        self.assertIn('ac0_applied=', stage)
        self.assertIn('PARTIALLY_APPLIED_MISSING_FORMAT_VALUES', stage)
        self.assertIn('def maybe_flag_ac0(', m4)
        self.assertIn('header_metadata.setdefault("filter", {})["AC0"]', m4)

    def test_main_config_has_off_switch(self):
        import yaml
        config = yaml.safe_load((ROOT / 'config.yaml').read_text())
        self.assertIs(config['ac0_flagging'], False)
        self.assertTrue(config['af_annotation'])

    def test_ac0_only_updates_filters(self):
        source = (ROOT / 'AF_hail_pipeline' / 'M4_af_annotation.py').read_text()
        source = source[source.index('def _flag_when_ac_zero('):source.index('def maybe_flag_ac0(')]
        self.assertIn('filters=hl.if_else(', source)
        self.assertNotIn('.filter_rows(', source)
        self.assertNotIn('.filter_cols(', source)
        self.assertNotIn('GT=', source)
        self.assertIn('ac0_evaluable_tmp', source)
        self.assertIn('hl.agg.all(entry_complete)', source)


if __name__ == '__main__':
    unittest.main()
