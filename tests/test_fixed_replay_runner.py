"""Tiny synthetic bundles test the runner, not the scientific replay."""

from __future__ import annotations

import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import textwrap
import unittest
from unittest.mock import patch


RUNNER_PATH = Path(__file__).resolve().parents[1] / "tools" / "run_fixed_replay.py"
SPEC = importlib.util.spec_from_file_location("fixed_replay_runner_under_test", RUNNER_PATH)
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)
SENTINEL = "PRIVATE_REPLAY_LOG_AND_LABEL_SENTINEL"


def success_source(
    *, rows: int = 238, summary_change: dict | None = None,
    failed_row: bool = False, difference: str = "0.0", after: str = "",
    header: list[str] | None = None,
) -> str:
    summary = {
        "checks": 238, "passed": 238, "failed": 0,
        "new_random_draws": 0, "new_model_queries": 0,
        "absolute_tolerance": 2e-12,
        "extra_untrusted_field": SENTINEL,
    }
    summary.update(summary_change or {})
    return textwrap.dedent(f"""\
        import csv, json, os, sys
        from pathlib import Path
        root = Path(__file__).resolve().parent
        assert Path.cwd() != root
        assert not (root / 'output').exists()
        assert sys.flags.optimize == 0
        assert all(k not in os.environ for k in ('PYTHONPATH', 'PYTHONHOME', 'PYTHONOPTIMIZE', 'PYTHONSTARTUP'))
        assert Path(os.environ['HOME']).parent == root.parent
        assert all(os.environ[k] == '1' for k in {runner.THREAD_VARIABLES!r})
        print({SENTINEL!r})
        print({SENTINEL!r}, file=sys.stderr)
        output = root / 'output'
        output.mkdir()
        (output / 'replay_summary.json').write_text(json.dumps({summary!r}), encoding='utf-8')
        with (output / 'replay_checks.csv').open('w', newline='', encoding='utf-8') as stream:
            writer = csv.writer(stream)
            writer.writerow({header if header is not None else runner.CSV_FIELDS!r})
            for i in range({rows!r}):
                writer.writerow([{SENTINEL!r}, not ({failed_row!r} and i == 0), {difference!r}])
        """) + after + "\n"


class FixedReplayRunnerTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="fixed-runner-tests-")
        self.root = Path(self.temporary.name).resolve()
        self.bundle = self.root / "source_bundle"
        self.bundle.mkdir()
        self.output = self.root / "safe_report"

    def tearDown(self):
        self.temporary.cleanup()

    def write_entry(self, source: str):
        (self.bundle / "replay.py").write_text(source, encoding="utf-8")

    def run_entry(self, source: str, timeout: float = 5):
        self.write_entry(source)
        return runner.run_fixed_replay(self.bundle, self.output, timeout)

    def test_pass_is_sanitized_fresh_and_source_is_unchanged(self):
        old_output = self.bundle / "output"
        old_output.mkdir()
        (old_output / "replay_summary.json").write_text("stale sentinel", encoding="utf-8")
        source = success_source()
        self.write_entry(source)
        with patch.dict(os.environ, {"PYTHONPATH": "not-a-real-search-path", "PYTHONOPTIMIZE": "2"}):
            report = runner.run_fixed_replay(self.bundle, self.output, 5)
        self.assertEqual(report["status"], "entry_execution_pass")
        self.assertEqual(report["isolation_assessment"], "external runner must attest")
        self.assertEqual(report["absolute_tolerance"], 2e-12)
        self.assertEqual((self.bundle / "replay.py").read_text(encoding="utf-8"), source)
        self.assertEqual((old_output / "replay_summary.json").read_text(encoding="utf-8"), "stale sentinel")
        self.assertEqual({p.name for p in self.output.iterdir()}, {"replay_summary.json", "replay_checks_safe.csv"})
        for path in self.output.iterdir():
            self.assertNotIn(SENTINEL, path.read_text(encoding="utf-8"))
        with (self.output / "replay_checks_safe.csv").open(newline="") as stream:
            reader = csv.DictReader(stream)
            rows = list(reader)
        self.assertEqual(reader.fieldnames, ["check_id", "passed", "max_absolute_difference"])
        self.assertEqual(len(rows), 238)
        self.assertEqual([int(row["check_id"]) for row in rows], list(range(1, 239)))

    def test_stale_outputs_cannot_pass(self):
        old = self.bundle / "output"
        old.mkdir()
        (old / "replay_summary.json").write_text(json.dumps({"checks": 238, "passed": 238}))
        (old / "replay_checks.csv").write_text("stale")
        report = self.run_entry("pass\n")
        self.assertEqual(report["error_code"], "fresh_output_missing")
        self.assertFalse((self.output / "replay_checks_safe.csv").exists())

    def test_nonzero_exit_blocks_even_valid_outputs(self):
        report = self.run_entry(success_source(after="raise SystemExit(7)"))
        self.assertEqual(report["error_code"], "entry_nonzero_exit")
        self.assertEqual(report["entry_return_code"], 7)
        self.assertFalse((self.output / "replay_checks_safe.csv").exists())

    def test_failed_summary_blocks(self):
        report = self.run_entry(success_source(summary_change={"passed": 237, "failed": 1}))
        self.assertEqual(report["error_code"], "summary_acceptance_failed")

    def test_failed_csv_row_blocks_even_passing_summary(self):
        report = self.run_entry(success_source(failed_row=True))
        self.assertEqual(report["error_code"], "checks_not_all_passed")

    def test_missing_csv_column_blocks(self):
        report = self.run_entry(success_source(header=["check", "passed"]))
        self.assertEqual(report["error_code"], "checks_schema_invalid")

    def test_wrong_row_count_blocks(self):
        report = self.run_entry(success_source(rows=237))
        self.assertEqual(report["error_code"], "checks_count_invalid")

    def test_missing_summary_acceptance_field_blocks(self):
        source = success_source(after=textwrap.dedent("""\
            summary_path = output / 'replay_summary.json'
            value = json.loads(summary_path.read_text())
            del value['new_model_queries']
            summary_path.write_text(json.dumps(value))
            """))
        report = self.run_entry(source)
        self.assertEqual(report["error_code"], "summary_acceptance_failed")

    def test_tolerance_cannot_be_relaxed(self):
        report = self.run_entry(success_source(summary_change={"absolute_tolerance": 1e-9}))
        self.assertEqual(report["error_code"], "summary_tolerance_invalid")

    def test_difference_cannot_exceed_original_tolerance(self):
        report = self.run_entry(success_source(difference="2.0001e-12"))
        self.assertEqual(report["error_code"], "checks_difference_invalid")

    def test_nonfinite_difference_blocks(self):
        report = self.run_entry(success_source(difference="nan"))
        self.assertEqual(report["error_code"], "checks_difference_invalid")

    def test_timeout_blocks(self):
        report = self.run_entry("import time\ntime.sleep(60)\n", timeout=0.15)
        self.assertEqual(report["error_code"], "entry_timeout")

    def test_symlink_bundle_entry_is_rejected(self):
        self.write_entry(success_source())
        try:
            (self.bundle / "linked_entry").symlink_to(self.bundle / "replay.py")
        except (OSError, NotImplementedError):
            self.skipTest("This runtime does not allow symlink creation")
        report = runner.run_fixed_replay(self.bundle, self.output, 5)
        self.assertEqual(report["error_code"], "symlink_or_junction_rejected")

    def test_existing_report_directory_is_not_overwritten(self):
        self.output.mkdir()
        marker = self.output / "keep.txt"
        marker.write_text("unchanged")
        self.write_entry(success_source())
        with self.assertRaises(runner.ReplayError) as caught:
            runner.run_fixed_replay(self.bundle, self.output)
        self.assertEqual(caught.exception.code, "output_directory_already_exists")
        self.assertEqual(marker.read_text(), "unchanged")

    def test_report_inside_source_bundle_is_refused(self):
        self.write_entry(success_source())
        with self.assertRaises(runner.ReplayError) as caught:
            runner.run_fixed_replay(self.bundle, self.bundle / "new_output")
        self.assertEqual(caught.exception.code, "output_must_be_outside_source_bundle")
        self.assertFalse((self.bundle / "new_output").exists())

    def test_optimized_parent_cannot_disable_entry_assertions(self):
        self.write_entry(success_source(after="assert False, 'synthetic assertion failure'"))
        env = dict(os.environ, PYTHONOPTIMIZE="2")
        result = subprocess.run(
            [sys.executable, "-O", str(RUNNER_PATH), "--bundle", str(self.bundle),
             "--output-dir", str(self.output), "--timeout-seconds", "5"],
            env=env, capture_output=True, text=True, timeout=10,
        )
        self.assertEqual(result.returncode, 1)
        self.assertEqual(result.stderr, "")
        self.assertNotIn(SENTINEL, result.stdout)
        self.assertNotIn(str(self.root), result.stdout)
        report = json.loads(result.stdout)
        self.assertEqual(report["error_code"], "entry_nonzero_exit")


if __name__ == "__main__":
    unittest.main()
