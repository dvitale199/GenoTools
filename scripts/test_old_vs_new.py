#!/usr/bin/env python3
"""Compare old vs new GenoTools pipeline implementations.

This script runs both the old (genotools) and new (genotools-new) CLI
with the same input data and compares the outputs.

Usage:
    python scripts/test_old_vs_new.py

Requirements:
    - Both genotools and genotools-new must be installed
    - Synthetic test data must exist in tests/data/synthetic/
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def run_command(cmd: list[str], env: dict | None = None) -> tuple[int, str, str]:
    """Run a command and return exit code, stdout, stderr."""
    full_env = os.environ.copy()
    if env:
        full_env.update(env)

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        env=full_env,
    )
    return result.returncode, result.stdout, result.stderr


def compare_json_files(file1: Path, file2: Path) -> dict:
    """Compare two JSON files and return differences."""
    with open(file1) as f:
        data1 = json.load(f)
    with open(file2) as f:
        data2 = json.load(f)

    differences = {}

    # Compare top-level keys
    keys1 = set(data1.keys())
    keys2 = set(data2.keys())

    if keys1 != keys2:
        differences["missing_in_new"] = list(keys1 - keys2)
        differences["extra_in_new"] = list(keys2 - keys1)

    # Compare common keys
    for key in keys1 & keys2:
        if data1[key] != data2[key]:
            # For pass_fail, do deeper comparison
            if "pass_fail" in key:
                diff = compare_pass_fail(data1[key], data2[key])
                if diff:
                    differences[key] = diff
            elif key == "QC":
                diff = compare_qc_metrics(data1[key], data2[key])
                if diff:
                    differences[key] = diff
            else:
                differences[key] = "values differ"

    return differences


def compare_pass_fail(pf1: dict, pf2: dict) -> dict | None:
    """Compare pass_fail dictionaries focusing on status."""
    diffs = {}

    for step in set(pf1.keys()) | set(pf2.keys()):
        if step not in pf1:
            diffs[step] = "missing in old"
        elif step not in pf2:
            diffs[step] = "missing in new"
        else:
            status1 = pf1[step].get("status")
            status2 = pf2[step].get("status")
            if status1 != status2:
                diffs[step] = f"status differs: old={status1}, new={status2}"

    return diffs if diffs else None


def compare_qc_metrics(qc1: dict, qc2: dict) -> dict | None:
    """Compare QC metrics, allowing for minor numeric differences."""
    # QC is a dict of lists/dicts from DataFrame.to_dict()
    # Focus on step names and pruned_count

    diffs = {}

    # Get unique steps from both
    steps1 = set(qc1.get("step", {}).values()) if "step" in qc1 else set()
    steps2 = set(qc2.get("step", {}).values()) if "step" in qc2 else set()

    if steps1 != steps2:
        diffs["steps_differ"] = {"old": list(steps1), "new": list(steps2)}

    # Compare pruned counts
    counts1 = qc1.get("pruned_count", {})
    counts2 = qc2.get("pruned_count", {})

    if counts1 != counts2:
        diffs["pruned_counts_differ"] = True

    return diffs if diffs else None


def count_pfile_samples(pfile_path: Path) -> int:
    """Count samples in a pfile by reading the psam."""
    psam_path = pfile_path.with_suffix(".psam")
    if not psam_path.exists():
        return -1

    with open(psam_path) as f:
        # Skip header, count lines
        lines = [l for l in f if not l.startswith("#")]
        # First line might be header without #
        if lines and lines[0].startswith("FID") or lines[0].startswith("#"):
            return len(lines) - 1
        return len(lines)


def main():
    # Paths
    project_root = Path(__file__).parent.parent
    test_data = project_root / "tests" / "data" / "synthetic" / "genotools_test"

    # Verify test data exists
    if not test_data.with_suffix(".pgen").exists():
        print(f"ERROR: Test data not found at {test_data}")
        sys.exit(1)

    print("=" * 60)
    print("GenoTools Old vs New Pipeline Comparison")
    print("=" * 60)
    print(f"Test data: {test_data}")
    print()

    # Create temporary directory for outputs
    with tempfile.TemporaryDirectory(prefix="genotools_test_") as tmpdir:
        tmpdir = Path(tmpdir)
        old_out = tmpdir / "old" / "output"
        new_out = tmpdir / "new" / "output"
        new_modules_out = tmpdir / "new_modules" / "output"

        # Create output directories
        old_out.parent.mkdir(parents=True)
        new_out.parent.mkdir(parents=True)
        new_modules_out.parent.mkdir(parents=True)

        # Test 1: Old pipeline (genotools)
        print("Running: OLD pipeline (genotools --all_sample)")
        print("-" * 60)
        cmd_old = [
            "genotools",
            "--pfile", str(test_data),
            "--out", str(old_out),
            "--all_sample",
            "--full_output", "True",  # Use full output to match new pipeline
        ]
        print(f"Command: {' '.join(cmd_old)}")
        exit_old, stdout_old, stderr_old = run_command(cmd_old)
        print(f"Exit code: {exit_old}")
        if exit_old != 0:
            print(f"STDERR:\n{stderr_old[:2000]}")
        print()

        # Test 2: New pipeline with legacy modules (genotools-new)
        print("Running: NEW CLI with LEGACY modules (genotools-new --all-sample)")
        print("-" * 60)
        cmd_new = [
            "genotools-new",
            "--pfile", str(test_data),
            "--out", str(new_out),
            "--all-sample",
            "--full-output",  # Use full output for proper path handling
        ]
        print(f"Command: {' '.join(cmd_new)}")
        exit_new, stdout_new, stderr_new = run_command(cmd_new)
        print(f"Exit code: {exit_new}")
        if exit_new != 0:
            print(f"STDERR:\n{stderr_new[:2000]}")
        print()

        # Test 3: New pipeline with new modules (GENOTOOLS_USE_NEW_MODULES=1)
        print("Running: NEW CLI with NEW modules (GENOTOOLS_USE_NEW_MODULES=1)")
        print("-" * 60)
        cmd_new_modules = [
            "genotools-new",
            "--pfile", str(test_data),
            "--out", str(new_modules_out),
            "--all-sample",
            "--full-output",  # Use full output for proper path handling
        ]
        print(f"Command: GENOTOOLS_USE_NEW_MODULES=1 {' '.join(cmd_new_modules)}")
        exit_new_modules, stdout_new_modules, stderr_new_modules = run_command(
            cmd_new_modules,
            env={"GENOTOOLS_USE_NEW_MODULES": "1"}
        )
        print(f"Exit code: {exit_new_modules}")
        if exit_new_modules != 0:
            print(f"STDERR:\n{stderr_new_modules[:2000]}")
        print()

        # Results summary
        print("=" * 60)
        print("RESULTS SUMMARY")
        print("=" * 60)

        # Check output files exist
        old_json = old_out.with_suffix(".json")
        new_json = new_out.with_suffix(".json")
        new_modules_json = new_modules_out.with_suffix(".json")

        old_pgen = old_out.with_suffix(".pgen")
        new_pgen = new_out.with_suffix(".pgen")
        new_modules_pgen = new_modules_out.with_suffix(".pgen")

        print("\nOutput files:")
        print(f"  Old JSON exists:         {old_json.exists()}")
        print(f"  New JSON exists:         {new_json.exists()}")
        print(f"  New modules JSON exists: {new_modules_json.exists()}")
        print(f"  Old pgen exists:         {old_pgen.exists()}")
        print(f"  New pgen exists:         {new_pgen.exists()}")
        print(f"  New modules pgen exists: {new_modules_pgen.exists()}")

        # Count samples in output
        print("\nSample counts:")
        old_samples = count_pfile_samples(old_out)
        new_samples = count_pfile_samples(new_out)
        new_modules_samples = count_pfile_samples(new_modules_out)
        print(f"  Old:         {old_samples}")
        print(f"  New:         {new_samples}")
        print(f"  New modules: {new_modules_samples}")

        # Compare JSON outputs
        print("\nJSON Comparison (Old vs New CLI w/ legacy modules):")
        if old_json.exists() and new_json.exists():
            diffs = compare_json_files(old_json, new_json)
            if diffs:
                print("  DIFFERENCES FOUND:")
                for key, val in diffs.items():
                    print(f"    {key}: {val}")
            else:
                print("  IDENTICAL (or equivalent)")
        else:
            print("  Cannot compare - missing files")

        print("\nJSON Comparison (Old vs New CLI w/ new modules):")
        if old_json.exists() and new_modules_json.exists():
            diffs = compare_json_files(old_json, new_modules_json)
            if diffs:
                print("  DIFFERENCES FOUND:")
                for key, val in diffs.items():
                    print(f"    {key}: {val}")
            else:
                print("  IDENTICAL (or equivalent)")
        else:
            print("  Cannot compare - missing files")

        # Final verdict
        print("\n" + "=" * 60)
        print("VERDICT")
        print("=" * 60)

        all_passed = True

        if exit_old != 0:
            print("FAIL: Old pipeline exited with error")
            all_passed = False

        if exit_new != 0:
            print("FAIL: New CLI (legacy modules) exited with error")
            all_passed = False

        if exit_new_modules != 0:
            print("FAIL: New CLI (new modules) exited with error")
            all_passed = False

        if old_samples != new_samples:
            print(f"WARN: Sample counts differ between old ({old_samples}) and new ({new_samples})")

        if old_samples != new_modules_samples:
            print(f"WARN: Sample counts differ between old ({old_samples}) and new_modules ({new_modules_samples})")

        if all_passed and old_samples == new_samples == new_modules_samples:
            print("PASS: All pipelines completed successfully with matching sample counts")
            return 0
        elif all_passed:
            print("PARTIAL: Pipelines completed but outputs may differ")
            return 1
        else:
            print("FAIL: One or more pipelines failed")
            return 1


if __name__ == "__main__":
    sys.exit(main())
