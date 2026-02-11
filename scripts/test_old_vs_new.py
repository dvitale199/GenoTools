#!/usr/bin/env python3
"""Compare old vs new GenoTools pipeline implementations.

This script runs both the old (genotools) and new (genotools-new) CLI
with the same input data and compares the outputs across multiple test
scenarios: sample QC, sample+variant QC, and full pipeline with ancestry.

QC scenarios run in parallel. Ancestry scenarios run sequentially to avoid
race conditions on shared ref panel intermediate files.

Usage:
    # QC only
    python scripts/test_old_vs_new.py --geno /path/to/data --out /path/to/output

    # QC + ancestry (default if ref panel exists)
    python scripts/test_old_vs_new.py

Requirements:
    - Both genotools and genotools-new must be installed
"""

import argparse
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
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
    diffs = {}

    # Get unique steps from both
    steps1 = set(qc1.get("step", {}).values()) if "step" in qc1 else set()
    steps2 = set(qc2.get("step", {}).values()) if "step" in qc2 else set()

    if steps1 != steps2:
        diffs["steps_differ"] = {"old": list(steps1), "new": list(steps2)}

    # Compare counts (old uses "pruned_count", new uses "count")
    counts1 = qc1.get("pruned_count", qc1.get("count", {}))
    counts2 = qc2.get("pruned_count", qc2.get("count", {}))

    if list(counts1.values()) != list(counts2.values()):
        diffs["counts_differ"] = True

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


def run_single_pipeline(
    job_id: str,
    cmd: list[str],
    out_path: Path,
    env: dict | None = None,
) -> dict:
    """Run a single pipeline invocation. Designed for use with ThreadPoolExecutor."""
    out_path.parent.mkdir(parents=True, exist_ok=True)

    start = time.time()
    exit_code, stdout, stderr = run_command(cmd, env=env)
    elapsed = time.time() - start

    return {
        "job_id": job_id,
        "exit": exit_code,
        "out": out_path,
        "stderr": stderr,
        "elapsed": elapsed,
        "cmd": cmd,
    }


def report_scenario(label: str, results: dict) -> bool:
    """Print comparison results for a scenario. Returns True if all passed."""
    print(f"\n  --- {label} ---")

    old_out = results["old"]["out"]
    new_out = results["new"]["out"]
    nm_out = results["new_modules"]["out"]

    # Timing
    for key in ("old", "new", "new_modules"):
        r = results[key]
        status = "OK" if r["exit"] == 0 else f"FAILED (exit {r['exit']})"
        print(f"  {key:>12s}: {status}  ({r['elapsed']:.1f}s)")

    # Check output files
    old_json = old_out.with_suffix(".json")
    new_json = new_out.with_suffix(".json")
    nm_json = nm_out.with_suffix(".json")

    print(f"  JSON files: old={old_json.exists()}, new={new_json.exists()}, new_modules={nm_json.exists()}")
    print(f"  Pgen files: old={old_out.with_suffix('.pgen').exists()}, new={new_out.with_suffix('.pgen').exists()}, new_modules={nm_out.with_suffix('.pgen').exists()}")

    # Sample counts
    old_samples = count_pfile_samples(old_out)
    new_samples = count_pfile_samples(new_out)
    nm_samples = count_pfile_samples(nm_out)
    print(f"  Sample counts: old={old_samples}, new={new_samples}, new_modules={nm_samples}")

    # JSON comparison
    if old_json.exists() and new_json.exists():
        diffs = compare_json_files(old_json, new_json)
        if diffs:
            print(f"  JSON diff (old vs new-legacy): {diffs}")
        else:
            print(f"  JSON diff (old vs new-legacy): IDENTICAL")

    if old_json.exists() and nm_json.exists():
        diffs = compare_json_files(old_json, nm_json)
        if diffs:
            print(f"  JSON diff (old vs new-modules): {diffs}")
        else:
            print(f"  JSON diff (old vs new-modules): IDENTICAL")

    # Verdict for this scenario
    all_ok = True
    if results["old"]["exit"] != 0:
        print(f"  FAIL: Old pipeline exited with error")
        all_ok = False
    if results["new"]["exit"] != 0:
        print(f"  FAIL: New CLI (legacy modules) exited with error")
        all_ok = False
    if results["new_modules"]["exit"] != 0:
        print(f"  FAIL: New CLI (new modules) exited with error")
        all_ok = False
    if old_samples != new_samples:
        print(f"  WARN: Sample counts differ: old={old_samples} vs new={new_samples}")
    if old_samples != nm_samples:
        print(f"  WARN: Sample counts differ: old={old_samples} vs new_modules={nm_samples}")

    if all_ok and old_samples == new_samples == nm_samples:
        print(f"  PASS")
    elif all_ok:
        print(f"  PARTIAL: Pipelines completed but outputs may differ")
        all_ok = False
    else:
        print(f"  FAIL")

    return all_ok


def parse_args():
    project_root = Path(__file__).resolve().parent.parent
    home = Path.home()

    parser = argparse.ArgumentParser(
        description="Compare old vs new GenoTools pipeline implementations."
    )
    parser.add_argument(
        "--geno", type=Path,
        default=project_root / "tests" / "data" / "synthetic" / "genotools_test",
        help="Input genotype file prefix (without extension)"
    )
    parser.add_argument(
        "--out", type=Path,
        default=Path.cwd() / "test_old_vs_new_results",
        help="Output directory for comparison results"
    )
    parser.add_argument(
        "--ref-panel", type=Path,
        default=home / ".genotools" / "ref" / "ref_panel" / "ref_panel_gp2_prune_rm_underperform_pos_update",
        help="Reference panel file prefix for ancestry testing (without extension)"
    )
    parser.add_argument(
        "--ref-labels", type=Path,
        default=home / ".genotools" / "ref" / "ref_panel" / "ref_panel_ancestry_updated.txt",
        help="Reference panel ancestry labels file (FID IID label, no header)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    test_data = args.geno
    out_dir = args.out

    # Verify test data exists
    if not test_data.with_suffix(".pgen").exists() and not test_data.with_suffix(".bed").exists():
        print(f"ERROR: Test data not found at {test_data}")
        sys.exit(1)

    # Determine which scenarios to run
    ref_panel_exists = args.ref_panel and (
        args.ref_panel.with_suffix(".bed").exists() or args.ref_panel.with_suffix(".pgen").exists()
    )
    ref_labels_exists = args.ref_labels and args.ref_labels.exists()
    run_ancestry = ref_panel_exists and ref_labels_exists
    if not run_ancestry:
        print("NOTE: Skipping ancestry scenario (ref panel/labels not found)")
        print(f"  ref_panel: {args.ref_panel} (found={ref_panel_exists})")
        print(f"  ref_labels: {args.ref_labels} (found={ref_labels_exists})")

    print("=" * 60)
    print("GenoTools Old vs New Pipeline Comparison")
    print("=" * 60)
    print(f"Test data:  {test_data}")
    print(f"Output dir: {out_dir}")
    if run_ancestry:
        print(f"Ref panel:  {args.ref_panel}")
        print(f"Ref labels: {args.ref_labels}")
    print()

    out_dir.mkdir(parents=True, exist_ok=True)

    # Build jobs: QC jobs run in parallel, ancestry jobs run sequentially
    # (ancestry writes intermediate files to shared ref panel directory)
    # Each job: (job_id, scenario_name, pipeline_variant, cmd, out_path, env)
    qc_jobs = []
    ancestry_jobs = []

    # --- Scenario 1: Sample QC ---
    s1 = out_dir / "sample_qc"
    qc_jobs.append(("sample_qc/old", "sample_qc", "old", [
        "genotools", "--pfile", str(test_data), "--out", str(s1 / "old" / "output"),
        "--full_output", "True", "--all_sample",
    ], s1 / "old" / "output", None))
    qc_jobs.append(("sample_qc/new", "sample_qc", "new", [
        "genotools-new", "--pfile", str(test_data), "--out", str(s1 / "new" / "output"),
        "--full-output", "--all-sample",
    ], s1 / "new" / "output", None))
    qc_jobs.append(("sample_qc/new_modules", "sample_qc", "new_modules", [
        "genotools-new", "--pfile", str(test_data), "--out", str(s1 / "new_modules" / "output"),
        "--full-output", "--all-sample",
    ], s1 / "new_modules" / "output", {"GENOTOOLS_USE_NEW_MODULES": "1"}))

    # --- Scenario 2: Sample + Variant QC ---
    s2 = out_dir / "full_qc"
    qc_jobs.append(("full_qc/old", "full_qc", "old", [
        "genotools", "--pfile", str(test_data), "--out", str(s2 / "old" / "output"),
        "--full_output", "True", "--all_sample", "--all_variant",
    ], s2 / "old" / "output", None))
    qc_jobs.append(("full_qc/new", "full_qc", "new", [
        "genotools-new", "--pfile", str(test_data), "--out", str(s2 / "new" / "output"),
        "--full-output", "--all-sample", "--all-variant",
    ], s2 / "new" / "output", None))
    qc_jobs.append(("full_qc/new_modules", "full_qc", "new_modules", [
        "genotools-new", "--pfile", str(test_data), "--out", str(s2 / "new_modules" / "output"),
        "--full-output", "--all-sample", "--all-variant",
    ], s2 / "new_modules" / "output", {"GENOTOOLS_USE_NEW_MODULES": "1"}))

    # --- Scenario 3: Full QC + Ancestry (sequential to avoid ref panel race) ---
    if run_ancestry:
        s3 = out_dir / "ancestry"
        ref_p = str(args.ref_panel)
        ref_l = str(args.ref_labels)
        ancestry_jobs.append(("ancestry/old", "ancestry", "old", [
            "genotools", "--pfile", str(test_data), "--out", str(s3 / "old" / "output"),
            "--full_output", "True", "--all_sample", "--all_variant",
            "--ancestry", "True", "--ref_panel", ref_p, "--ref_labels", ref_l,
        ], s3 / "old" / "output", None))
        ancestry_jobs.append(("ancestry/new", "ancestry", "new", [
            "genotools-new", "--pfile", str(test_data), "--out", str(s3 / "new" / "output"),
            "--full-output", "--all-sample", "--all-variant",
            "--ancestry", "--ref-panel", ref_p, "--ref-labels", ref_l,
        ], s3 / "new" / "output", None))
        ancestry_jobs.append(("ancestry/new_modules", "ancestry", "new_modules", [
            "genotools-new", "--pfile", str(test_data), "--out", str(s3 / "new_modules" / "output"),
            "--full-output", "--all-sample", "--all-variant",
            "--ancestry", "--ref-panel", ref_p, "--ref-labels", ref_l,
        ], s3 / "new_modules" / "output", {"GENOTOOLS_USE_NEW_MODULES": "1"}))

    all_jobs = qc_jobs + ancestry_jobs
    n_jobs = len(all_jobs)
    print(f"Running {len(qc_jobs)} QC jobs in parallel, {len(ancestry_jobs)} ancestry jobs sequentially...")
    print("-" * 60)
    for job_id, _, _, cmd, _, env in all_jobs:
        env_prefix = "GENOTOOLS_USE_NEW_MODULES=1 " if env else ""
        print(f"  {job_id}: {env_prefix}{' '.join(cmd)}")
    print("-" * 60)
    print()

    total_start = time.time()
    scenario_results: dict[str, dict] = {}

    # Phase 1: Run QC jobs in parallel
    if qc_jobs:
        print(f"Phase 1: Running {len(qc_jobs)} QC jobs in parallel...")
        with ThreadPoolExecutor(max_workers=len(qc_jobs)) as executor:
            futures = {}
            for job_id, scenario, variant, cmd, out_path, env in qc_jobs:
                future = executor.submit(run_single_pipeline, job_id, cmd, out_path, env)
                futures[future] = (scenario, variant)

            for future in as_completed(futures):
                scenario, variant = futures[future]
                result = future.result()
                status = "OK" if result["exit"] == 0 else "FAILED"
                print(f"  Completed: {result['job_id']} [{status}] ({result['elapsed']:.1f}s)")

                if scenario not in scenario_results:
                    scenario_results[scenario] = {}
                scenario_results[scenario][variant] = result

    # Phase 2: Run ancestry jobs sequentially (shared ref panel intermediates)
    if ancestry_jobs:
        print(f"\nPhase 2: Running {len(ancestry_jobs)} ancestry jobs sequentially...")
        for job_id, scenario, variant, cmd, out_path, env in ancestry_jobs:
            print(f"  Running: {job_id}...")
            result = run_single_pipeline(job_id, cmd, out_path, env)
            status = "OK" if result["exit"] == 0 else "FAILED"
            print(f"  Completed: {result['job_id']} [{status}] ({result['elapsed']:.1f}s)")

            if scenario not in scenario_results:
                scenario_results[scenario] = {}
            scenario_results[scenario][variant] = result

    total_elapsed = time.time() - total_start
    print(f"\nAll {n_jobs} runs completed in {total_elapsed:.1f}s")

    # Summary
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)

    all_passed = True
    for name in ["sample_qc", "full_qc", "ancestry"]:
        if name in scenario_results:
            passed = report_scenario(name, scenario_results[name])
            if not passed:
                all_passed = False

    # Print stderr for any failures
    any_failures = False
    for name, results in scenario_results.items():
        for variant, r in results.items():
            if r["exit"] != 0:
                if not any_failures:
                    print("\n" + "=" * 60)
                    print("FAILURE DETAILS")
                    print("=" * 60)
                    any_failures = True
                print(f"\n  {r['job_id']} (exit {r['exit']}):")
                print(f"  {r['stderr'][:3000]}")

    # Final verdict
    print("\n" + "=" * 60)
    print("FINAL VERDICT")
    print("=" * 60)

    if all_passed:
        print(f"PASS: All scenarios completed successfully with matching outputs ({total_elapsed:.1f}s)")
        return 0
    else:
        print(f"FAIL: One or more scenarios had differences or errors ({total_elapsed:.1f}s)")
        return 1


if __name__ == "__main__":
    sys.exit(main())
