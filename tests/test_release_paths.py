"""Regression checks for orchestration and source-only re-export contracts."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

SCRIPTS=Path(__file__).resolve().parents[1]/'NHANES_Q3_Project/scripts'
def load(name):
    spec=importlib.util.spec_from_file_location(name,SCRIPTS/(name+'.py'))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

class ReleasePaths(unittest.TestCase):
    def test_new_predictions_force_new_evaluation_despite_old_summary_and_marker(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);out=root/'outputs'
            (out/'tables').mkdir(parents=True);(out/'logs').mkdir()
            summary=out/'tables/fixed_substitution_summary.json'
            summary.write_text('{"prediction_version":"old"}')
            (out/'logs/output_generation_complete.json').write_text('{}')
            (out/'checks').mkdir()
            for n in ['preflight','artificial_modeling']:(out/'checks'/f'{n}.json').write_text('{}')
            runtime=types.SimpleNamespace(limit_threads=lambda:None,ROOT=root,OUT=out,
                read_json=lambda p:{'cadence_amendment_utc':'recorded'},event=lambda *a,**k:None,
                write_json=lambda *a,**k:None,now=lambda:None)
            with patch.dict(sys.modules,{'runtime':runtime}):module=load('run_module')
            calls=[]
            def command(script):
                calls.append(script)
                if script=='run_benchmark.py':(out/'new_prediction.json').write_text('{"prediction_version":"new"}')
                if script=='evaluate_fixed_predictions.py':summary.write_bytes((out/'new_prediction.json').read_bytes())
                if script=='make_outputs.py':self.assertEqual(json.loads(summary.read_text())['prediction_version'],'new')
            with patch.object(module,'command',command),patch.object(sys,'argv',['run_module.py']):module.main()
            self.assertEqual(calls,['run_benchmark.py','evaluate_fixed_predictions.py','make_outputs.py'])
            self.assertEqual(json.loads(summary.read_text())['prediction_version'],'new')

    def test_failed_evaluation_never_renders_stale_results(self):
        runtime=types.SimpleNamespace(limit_threads=lambda:None,ROOT=Path('.'),OUT=Path('.'),
            read_json=lambda p:{'cadence_amendment_utc':'recorded'},event=lambda *a,**k:None,write_json=None,now=None)
        with patch.dict(sys.modules,{'runtime':runtime}):module=load('run_module')
        calls=[]
        def command(script):
            calls.append(script)
            if script=='evaluate_fixed_predictions.py':raise RuntimeError('synthetic evaluation failure')
        with patch.object(module,'command',command),patch.object(sys,'argv',['run_module.py']):
            with self.assertRaises(RuntimeError):module.main()
        self.assertNotIn('make_outputs.py',calls)

    def test_report_only_never_runs_training_or_evaluation(self):
        runtime=types.SimpleNamespace(limit_threads=lambda:None,ROOT=Path('.'),OUT=Path('.'),
            read_json=lambda p:self.fail('report read protocol'),event=None,write_json=None,now=None)
        seen=[]
        with patch.dict(sys.modules,{'runtime':runtime,'saved_reports':types.SimpleNamespace(generate=seen.append)}):
            module=load('run_module')
            for flag in ['--report-only','--outputs-only']:
                with patch.object(module,'command',lambda x:self.fail('report entered computation')),patch.object(sys,'argv',['run_module.py',flag,'--report-dir','saved_report']):module.main()
        self.assertEqual(seen,['saved_report','saved_report'])

    def test_reexport_does_not_require_or_redistribute_official_r_reader(self):
        module=load('export_code_only')
        with tempfile.TemporaryDirectory() as tmp:
            destination=Path(tmp)/'export'
            module.export(destination)
            self.assertFalse((destination/'sources/R_ReadInProgramAllSurveys.R').exists())
            self.assertTrue((destination/'scripts/prepare_cohort.py').is_file())
            self.assertIn('https://ftp.cdc.gov/',(destination/'sources/SOURCES.md').read_text())
            self.assertFalse((destination/'reference').exists())

if __name__=='__main__':unittest.main()
