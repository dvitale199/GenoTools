#!/usr/bin/env python3
"""A/B test: legacy ancestry.py vs new ancestry/ module.

Runs both genotools (legacy ancestry) and genotools-new (new AncestryModel)
with identical inputs and compares the predicted ancestry labels.

Usage:
    python tests/scripts/test_ancestry_ab.py \
        --pfile /path/to/data \
        --ref-panel /path/to/ref \
        --ref-labels /path/to/labels.txt \
        --out-dir /path/to/results

Both pipelines use the same PLINK preprocessing (get_raw_files). Only the
ML pipeline differs: legacy uses Ancestry class methods, new uses AncestryModel.
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path

import pandas as pd


def run_pipeline(cmd: list, label: str) -> int:
    """Run a pipeline command and return exit code."""
    print(f"\n{'=' * 60}")
    print(f"Running {label}")
    print(f"Command: {' '.join(cmd)}")
    print(f"{'=' * 60}\n")

    result = subprocess.run(cmd)

    print(f"\n{label} exit code: {result.returncode}")
    return result.returncode


def find_labels_file(out_prefix: str) -> Path:
    """Find the predicted labels file for an output prefix."""
    candidates = [
        Path(f"{out_prefix}_ancestry_umap_linearsvc_predicted_labels.txt"),
        Path(f"{out_prefix}_umap_linearsvc_predicted_labels.txt"),
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(
        f"Labels file not found. Tried:\n"
        + "\n".join(f"  {p}" for p in candidates)
    )


def compare_predictions(legacy_path: Path, new_path: Path) -> dict:
    """Compare ancestry predictions from legacy and new pipelines."""
    legacy_df = pd.read_csv(legacy_path, sep="\t")
    new_df = pd.read_csv(new_path, sep="\t")

    # Normalize column names
    if "label" in legacy_df.columns and "predicted_ancestry" not in legacy_df.columns:
        legacy_df = legacy_df.rename(columns={"label": "predicted_ancestry"})
    if "label" in new_df.columns and "predicted_ancestry" not in new_df.columns:
        new_df = new_df.rename(columns={"label": "predicted_ancestry"})

    # Merge on IID
    merged = legacy_df.merge(new_df, on="IID", suffixes=("_legacy", "_new"))
    total = len(merged)

    if total == 0:
        return {
            "total_samples": 0,
            "error": "No overlapping samples between legacy and new",
        }

    # Per-sample agreement
    agreement = (
        merged["predicted_ancestry_legacy"] == merged["predicted_ancestry_new"]
    ).sum()
    agreement_pct = agreement / total * 100

    # Disagreements
    disagreements = merged[
        merged["predicted_ancestry_legacy"] != merged["predicted_ancestry_new"]
    ]

    # Ancestry counts comparison
    legacy_counts = legacy_df["predicted_ancestry"].value_counts().to_dict()
    new_counts = new_df["predicted_ancestry"].value_counts().to_dict()

    all_labels = sorted(set(list(legacy_counts.keys()) + list(new_counts.keys())))
    count_comparison = {}
    for label in all_labels:
        count_comparison[label] = {
            "legacy": legacy_counts.get(label, 0),
            "new": new_counts.get(label, 0),
            "diff": new_counts.get(label, 0) - legacy_counts.get(label, 0),
        }

    # Build disagreement details
    if len(disagreements) <= 100:
        disagreement_details = disagreements[
            ["IID", "predicted_ancestry_legacy", "predicted_ancestry_new"]
        ].to_dict("records")
    else:
        disagreement_details = (
            f"{len(disagreements)} disagreements (showing first 100)"
        )
        disagreement_details = disagreements.head(100)[
            ["IID", "predicted_ancestry_legacy", "predicted_ancestry_new"]
        ].to_dict("records")

    return {
        "total_samples": total,
        "agreement_count": int(agreement),
        "agreement_pct": round(float(agreement_pct), 2),
        "disagreement_count": int(total - agreement),
        "count_comparison": count_comparison,
        "disagreement_details": disagreement_details,
    }


def compare_split_pfiles(legacy_prefix: str, new_prefix: str, labels: list) -> dict:
    """Compare split pfiles by sample count per ancestry."""
    results = {}
    for label in labels:
        legacy_psam = Path(f"{legacy_prefix}_{label}.psam")
        new_psam = Path(f"{new_prefix}_{label}.psam")

        legacy_n = 0
        new_n = 0

        if legacy_psam.exists():
            legacy_n = len(pd.read_csv(legacy_psam, sep=r"\s+"))
        if new_psam.exists():
            new_n = len(pd.read_csv(new_psam, sep=r"\s+"))

        results[label] = {
            "legacy_samples": legacy_n,
            "new_samples": new_n,
            "match": legacy_n == new_n,
        }

    return results


def print_report(results: dict) -> None:
    """Print comparison report."""
    print(f"\n{'=' * 60}")
    print("A/B TEST RESULTS: Legacy vs New Ancestry")
    print(f"{'=' * 60}")
    print(f"Total samples:    {results['total_samples']}")
    print(f"Agreement:        {results['agreement_count']} ({results['agreement_pct']}%)")
    print(f"Disagreements:    {results['disagreement_count']}")
    print()
    print("Ancestry counts:")
    print(f"  {'Label':<10} {'Legacy':>8} {'New':>8} {'Diff':>8}")
    print(f"  {'-' * 36}")
    for label, counts in sorted(results["count_comparison"].items()):
        diff_str = f"{counts['diff']:+d}" if counts["diff"] != 0 else "0"
        print(f"  {label:<10} {counts['legacy']:>8} {counts['new']:>8} {diff_str:>8}")

    if results["disagreement_count"] > 0 and isinstance(
        results["disagreement_details"], list
    ):
        print(f"\nDisagreement details (first {min(20, len(results['disagreement_details']))}):")
        for d in results["disagreement_details"][:20]:
            print(
                f"  {d['IID']}: {d['predicted_ancestry_legacy']} -> {d['predicted_ancestry_new']}"
            )


def main():
    parser = argparse.ArgumentParser(
        description="A/B test legacy vs new ancestry prediction"
    )
    parser.add_argument("--pfile", type=str, required=True, help="Input pfile prefix")
    parser.add_argument("--ref-panel", type=str, required=True, help="Reference panel prefix")
    parser.add_argument("--ref-labels", type=str, required=True, help="Reference labels file")
    parser.add_argument("--out-dir", type=str, required=True, help="Output directory")
    parser.add_argument(
        "--extra-args",
        nargs="*",
        default=[],
        help="Extra args to pass to both commands",
    )
    parser.add_argument(
        "--skip-legacy",
        action="store_true",
        help="Skip running legacy (if output already exists)",
    )
    parser.add_argument(
        "--skip-new",
        action="store_true",
        help="Skip running new (if output already exists)",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    legacy_out = str(out_dir / "legacy")
    new_out = str(out_dir / "new")

    # Common args for both pipelines
    common = [
        "--pfile", args.pfile,
        "--ref-panel", args.ref_panel,
        "--ref-labels", args.ref_labels,
        "--ancestry",
        "--full-output",
    ] + args.extra_args

    # Run legacy pipeline
    rc1 = 0
    if not args.skip_legacy:
        rc1 = run_pipeline(
            ["genotools", "--out", legacy_out] + common, "LEGACY (ancestry.py)"
        )
    else:
        print("Skipping legacy run (--skip-legacy)")

    # Run new pipeline
    rc2 = 0
    if not args.skip_new:
        rc2 = run_pipeline(
            ["genotools-new", "--out", new_out] + common, "NEW (AncestryModel)"
        )
    else:
        print("Skipping new run (--skip-new)")

    # Check exit codes
    if rc1 != 0:
        print(f"\nLEGACY pipeline FAILED (exit code {rc1})")
    if rc2 != 0:
        print(f"\nNEW pipeline FAILED (exit code {rc2})")
    if rc1 != 0 or rc2 != 0:
        sys.exit(1)

    # Find and compare prediction files
    try:
        legacy_labels = find_labels_file(legacy_out)
        new_labels = find_labels_file(new_out)
    except FileNotFoundError as e:
        print(f"\nERROR: {e}")
        sys.exit(1)

    print(f"\nLegacy labels: {legacy_labels}")
    print(f"New labels:    {new_labels}")

    # Compare predictions
    results = compare_predictions(legacy_labels, new_labels)

    # Print report
    print_report(results)

    # Save full results
    report_path = out_dir / "ab_test_report.json"
    with open(report_path, "w") as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nFull report saved to: {report_path}")

    # Exit status based on agreement
    if results.get("error"):
        print(f"\nERROR: {results['error']}")
        sys.exit(1)
    elif results["agreement_pct"] < 80:
        print("\nFAIL: Agreement below 80% -- investigate differences")
        sys.exit(2)
    elif results["agreement_pct"] < 95:
        print("\nWARN: Agreement below 95% -- minor differences detected")
        sys.exit(0)
    else:
        print("\nPASS: High agreement between legacy and new ancestry")
        sys.exit(0)


if __name__ == "__main__":
    main()
