#!/usr/bin/env python3
"""Run and validate one fixed replay entry; filesystem isolation is external.

The source bundle is never imported or modified. Its root ``output`` directory
is excluded from a temporary copy. Only a field-validated summary and a CSV
without original check labels leave the temporary tree. This runner is neither
a sandbox nor a privacy scanner for the rest of a release package.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import shutil
import signal
import stat
import subprocess
import sys
import tempfile
from typing import Any


EXPECTED_CHECKS = 238
ABSOLUTE_TOLERANCE = 2e-12
ISOLATION_ASSESSMENT = "external runner must attest"
CSV_FIELDS = ["check", "passed", "max_absolute_difference"]
THREAD_VARIABLES = (
    "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS",
    "NUMBA_NUM_THREADS",
)


class ReplayError(Exception):
    """Only a fixed, non-sensitive code may be published from an error."""

    def __init__(self, code: str, return_code: int | None = None):
        super().__init__(code)
        self.code = code
        self.return_code = return_code


def _absolute(path: Path | str) -> Path:
    # Do not resolve symlinks before checking for them.
    return Path(os.path.abspath(os.path.expanduser(os.fspath(path))))


def _reject_link(path: Path) -> None:
    try:
        metadata = path.lstat()
    except FileNotFoundError:
        return
    reparse_flag = getattr(stat, "FILE_ATTRIBUTE_REPARSE_POINT", 0)
    is_reparse = bool(getattr(metadata, "st_file_attributes", 0) & reparse_flag)
    if stat.S_ISLNK(metadata.st_mode) or is_reparse:
        raise ReplayError("symlink_or_junction_rejected")


def _reject_link_components(path: Path) -> None:
    for part in reversed((path, *path.parents)):
        _reject_link(part)


def _tree_entries(root: Path) -> list[tuple[Path, bool]]:
    """Check the bounded bundle tree, including links in old output."""
    entries: list[tuple[Path, bool]] = []
    def unreadable_directory(_error: OSError) -> None:
        raise ReplayError("bundle_tree_unreadable")

    for directory, dirs, files in os.walk(root, followlinks=False, onerror=unreadable_directory):
        for name in sorted(dirs + files):
            path = Path(directory) / name
            _reject_link(path)
            mode = path.lstat().st_mode
            is_dir = stat.S_ISDIR(mode)
            if not is_dir and not stat.S_ISREG(mode):
                raise ReplayError("non_regular_bundle_entry")
            entries.append((path.relative_to(root), is_dir))
    return entries


def _copy_bundle(source: Path, destination: Path) -> None:
    entries = _tree_entries(source)
    destination.mkdir()
    for relative, is_dir in entries:
        if relative.parts[0] == "output" or "__pycache__" in relative.parts:
            continue
        if not is_dir and relative.suffix in {".pyc", ".pyo"}:
            continue
        source_path = source / relative
        target = destination / relative
        _reject_link_components(source_path)
        if is_dir:
            target.mkdir(parents=True, exist_ok=True)
        else:
            if not stat.S_ISREG(source_path.lstat().st_mode):
                raise ReplayError("non_regular_bundle_entry")
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source_path, target, follow_symlinks=False)
    # Also reject a symlink introduced during copying. This is not a sandbox.
    _tree_entries(destination)
    if os.path.lexists(destination / "output"):
        raise ReplayError("old_output_was_copied")


def _child_environment(work_root: Path) -> dict[str, str]:
    # An allowlist avoids inheriting PYTHONPATH/HOME/OPTIMIZE/STARTUP, virtualenv
    # hints, application credentials, and user cache/config locations.
    env = {"PATH": os.defpath}
    for key in ("SYSTEMROOT", "WINDIR", "COMSPEC", "SYSTEMDRIVE", "PATHEXT"):
        if key in os.environ:
            env[key] = os.environ[key]
    locations = {
        "HOME": "home", "USERPROFILE": "home", "TMPDIR": "tmp",
        "TEMP": "tmp", "TMP": "tmp", "XDG_CACHE_HOME": "cache",
        "XDG_CONFIG_HOME": "config", "XDG_DATA_HOME": "user_data",
        "MPLCONFIGDIR": "matplotlib", "JOBLIB_TEMP_FOLDER": "tmp",
        "NUMBA_CACHE_DIR": "numba",
    }
    for variable, directory in locations.items():
        location = work_root / directory
        location.mkdir(exist_ok=True)
        env[variable] = str(location)
    if os.name == "nt":
        home = Path(env["HOME"])
        env["HOMEDRIVE"] = home.drive
        env["HOMEPATH"] = str(home)[len(home.drive):]
    env.update({
        "PYTHONNOUSERSITE": "1", "PYTHONDONTWRITEBYTECODE": "1",
        "PYTHONUTF8": "1", "LANG": "C.UTF-8", "LC_ALL": "C.UTF-8",
    })
    env.update({key: "1" for key in THREAD_VARIABLES})
    return env


def _kill_entry(process: subprocess.Popen) -> None:
    if os.name == "posix":
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    elif process.poll() is None:
        process.kill()
    process.wait()


def _execute(copied: Path, work_root: Path, timeout_seconds: float) -> None:
    unrelated_cwd = work_root / "unrelated_cwd"
    unrelated_cwd.mkdir()
    command = [os.path.abspath(sys.executable), str(copied / "replay.py")]
    # Temporary raw streams are never read, printed, or copied to the report.
    with tempfile.TemporaryFile(dir=work_root) as stdout_file, \
            tempfile.TemporaryFile(dir=work_root) as stderr_file:
        try:
            process = subprocess.Popen(
                command, cwd=unrelated_cwd, env=_child_environment(work_root),
                stdin=subprocess.DEVNULL, stdout=stdout_file, stderr=stderr_file,
                start_new_session=(os.name == "posix"),
            )
        except OSError:
            raise ReplayError("entry_start_failed") from None
        try:
            return_code = process.wait(timeout=timeout_seconds)
        except subprocess.TimeoutExpired:
            _kill_entry(process)
            raise ReplayError("entry_timeout") from None
        except BaseException:
            _kill_entry(process)
            raise
        if return_code != 0:
            raise ReplayError("entry_nonzero_exit", return_code)


def _json_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise ReplayError("summary_duplicate_key")
        value[key] = item
    return value


def _invalid_json_constant(_value: str) -> None:
    raise ReplayError("summary_invalid_json")


def _regular_artifact(path: Path, size_limit: int) -> None:
    _reject_link_components(path)
    if not path.is_file():
        raise ReplayError("fresh_output_missing")
    metadata = path.stat()
    if not stat.S_ISREG(metadata.st_mode) or metadata.st_size > size_limit:
        raise ReplayError("output_type_or_size_invalid")


def validate_outputs(output: Path) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Validate the observed schema, never importing replay.py or its data."""
    summary_path = output / "replay_summary.json"
    checks_path = output / "replay_checks.csv"
    _regular_artifact(summary_path, 64 * 1024)
    _regular_artifact(checks_path, 4 * 1024 * 1024)
    try:
        summary = json.loads(
            summary_path.read_text(encoding="utf-8"),
            object_pairs_hook=_json_pairs, parse_constant=_invalid_json_constant,
        )
    except (ValueError, UnicodeError):
        raise ReplayError("summary_invalid_json") from None
    if not isinstance(summary, dict):
        raise ReplayError("summary_schema_invalid")
    expected = {
        "checks": EXPECTED_CHECKS, "passed": EXPECTED_CHECKS, "failed": 0,
        "new_random_draws": 0, "new_model_queries": 0,
    }
    for key, value in expected.items():
        if type(summary.get(key)) is not int or summary[key] != value:
            raise ReplayError("summary_acceptance_failed")
    tolerance = summary.get("absolute_tolerance")
    if type(tolerance) not in (int, float) or tolerance != ABSOLUTE_TOLERANCE:
        raise ReplayError("summary_tolerance_invalid")
    safe_rows: list[dict[str, Any]] = []
    try:
        with checks_path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != CSV_FIELDS:
                raise ReplayError("checks_schema_invalid")
            for index, row in enumerate(reader, start=1):
                if index > EXPECTED_CHECKS or set(row) != set(CSV_FIELDS):
                    raise ReplayError("checks_row_schema_invalid")
                if not isinstance(row["check"], str) or not row["check"].strip():
                    raise ReplayError("checks_label_missing")
                if row["passed"] != "True":
                    raise ReplayError("checks_not_all_passed")
                try:
                    difference = float(row["max_absolute_difference"])
                except (ValueError, TypeError):
                    raise ReplayError("checks_difference_invalid") from None
                if not math.isfinite(difference) or not 0 <= difference <= ABSOLUTE_TOLERANCE:
                    raise ReplayError("checks_difference_invalid")
                # Labels need not be unique in the real bundle and may contain
                # paths. Publish only fresh positional IDs and numeric results.
                safe_rows.append({
                    "check_id": index, "passed": True,
                    "max_absolute_difference": difference,
                })
    except (csv.Error, UnicodeError):
        raise ReplayError("checks_csv_invalid") from None
    if len(safe_rows) != EXPECTED_CHECKS:
        raise ReplayError("checks_count_invalid")
    return {**expected, "absolute_tolerance": ABSOLUTE_TOLERANCE}, safe_rows


def _failure(error: ReplayError) -> dict[str, Any]:
    report: dict[str, Any] = {
        "status": "entry_execution_failed",
        "isolation_assessment": ISOLATION_ASSESSMENT,
        "error_code": error.code,
    }
    if error.return_code is not None:
        report["entry_return_code"] = error.return_code
    return report


def run_fixed_replay(
    bundle: Path | str, output_dir: Path | str, timeout_seconds: float = 300,
) -> dict[str, Any]:
    """Return a safe report; refuse an existing output directory outright."""
    if isinstance(timeout_seconds, bool) or not isinstance(timeout_seconds, (int, float)) \
            or not math.isfinite(timeout_seconds) or timeout_seconds <= 0:
        raise ReplayError("timeout_invalid")
    bundle, output_dir = _absolute(bundle), _absolute(output_dir)
    _reject_link_components(bundle)
    _reject_link_components(output_dir)
    if output_dir == bundle or bundle in output_dir.parents:
        raise ReplayError("output_must_be_outside_source_bundle")
    if os.path.lexists(output_dir):
        raise ReplayError("output_directory_already_exists")
    output_dir.mkdir(parents=True, exist_ok=False)
    safe_rows: list[dict[str, Any]] = []
    try:
        if not bundle.is_dir() or not (bundle / "replay.py").is_file():
            raise ReplayError("bundle_or_entry_missing")
        with tempfile.TemporaryDirectory(prefix="fixed-replay-") as temporary:
            work_root = Path(temporary).resolve()
            copied = work_root / "bundle"
            _copy_bundle(bundle, copied)
            _execute(copied, work_root, timeout_seconds)
            counts, safe_rows = validate_outputs(copied / "output")
            report = {
                "status": "entry_execution_pass",
                "isolation_assessment": ISOLATION_ASSESSMENT,
                "entry": "replay.py",
                "entry_return_code": 0,
                **counts,
                "safe_checks_file": "replay_checks_safe.csv",
            }
    except ReplayError as error:
        report = _failure(error)
    except (OSError, shutil.Error):
        report = _failure(ReplayError("runner_filesystem_error"))
    if report["status"] == "entry_execution_pass":
        with (output_dir / "replay_checks_safe.csv").open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=[
                "check_id", "passed", "max_absolute_difference",
            ])
            writer.writeheader()
            writer.writerows(safe_rows)
    with (output_dir / "replay_summary.json").open("x", encoding="utf-8") as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write("\n")
    return report


class _SafeArgumentParser(argparse.ArgumentParser):
    def error(self, _message: str) -> None:
        raise ReplayError("invalid_arguments")


def main(argv: list[str] | None = None) -> int:
    parser = _SafeArgumentParser(description=__doc__)
    parser.add_argument("--bundle", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--timeout-seconds", type=float, default=300)
    try:
        arguments = parser.parse_args(argv)
        report = run_fixed_replay(arguments.bundle, arguments.output_dir, arguments.timeout_seconds)
    except ReplayError as error:
        report = _failure(error)
    except Exception:
        # Never emit an exception message/traceback containing local paths.
        report = _failure(ReplayError("runner_internal_error"))
    print(json.dumps(report, sort_keys=True, allow_nan=False))
    return 0 if report["status"] == "entry_execution_pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
