#!/usr/bin/env python
# Copyright 2023 The GenoTools Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================

"""Quantify what unpinning ``umap_learn`` does to ancestry calls.

``umap_learn==0.5.3`` was pinned for reproducibility, and unpinning it changes
the embedding that feeds the classifier — so it can move ancestry labels. This
runs the same cohort through two Python environments that differ only in their
installed libraries (the GenoTools source is identical, since both install it
editable from this working tree) and reports how many calls moved.

Two questions, which are not the same question
----------------------------------------------
``--mode retrain`` (default)
    Each environment trains its own model from the reference panel, then
    predicts. This is the real upgrade path: what a user gets after upgrading
    and retraining. It is the number that decides whether to unpin.

``--mode reuse-model``
    Both environments load the *same* pre-trained model and only predict (the
    reference panel is still required — prediction re-derives the reference PCA
    from it — but the classifier is not refitted). This
    isolates inference drift from training drift, and quantifies the silent
    failure that model version-recording now warns about: an old model run
    under a new umap. Much faster, since nothing is trained.

Running both tells you whether a difference came from the fit or the transform.

Usage
-----
    python tests/scripts/compare_umap_versions.py \\
        --baseline-python .venv/bin/python \\
        --candidate-python .venv-next/bin/python \\
        --geno ~/parity_data/GP2_r12_subset10k \\
        --ref-panel ~/.genotools/ref/ref_panel/ref_panel_gp2_prune_rm_underperform_pos_update \\
        --ref-labels ~/.genotools/ref/ref_panel/ref_panel_ancestry_updated.txt \\
        --work /tmp/umap-compare

A full ``--ancestry --all-sample`` run on the 10k subset takes roughly 25
minutes per environment, so ``--mode retrain`` is an hour-plus job. Completed
runs are reused on re-invocation unless ``--force`` is passed, so an
interrupted comparison can be resumed.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional

REPO_ROOT = Path(__file__).resolve().parents[2]
COMPARATOR = REPO_ROOT / "tests" / "scripts" / "compare_ancestry_run.py"

def environment_versions(python: Path) -> Dict[str, str]:
    """Return the library versions installed for one interpreter."""
    code = (
        "import json, platform\n"
        "import importlib.metadata as md\n"
        "out = {}\n"
        "for p in ['umap-learn','scikit-learn','xgboost','numpy','pandas','scipy','numba']:\n"
        "    try: out[p] = md.version(p)\n"
        "    except Exception: out[p] = 'absent'\n"
        "out['python'] = platform.python_version()\n"
        "print(json.dumps(out))\n"
    )
    result = subprocess.run(
        [str(python), "-c", code], capture_output=True, text=True, check=True
    )
    return json.loads(result.stdout)


def describe_difference(base: Dict[str, str], cand: Dict[str, str]) -> List[str]:
    """Render the environment delta, so attribution is not guesswork."""
    return [
        f"  {name:15s} {base[name]:>12s}  ->  {cand[name]}"
        for name in base
        if base.get(name) != cand.get(name)
    ]


def run_pipeline(
    python: Path,
    out_prefix: Path,
    geno: Path,
    ref_panel: Optional[Path],
    ref_labels: Optional[Path],
    model: Optional[Path],
    force: bool,
) -> None:
    """Run one ancestry pipeline, reusing a completed run unless forced.

    Reuse is keyed on the JSON report, which the runner writes last — a partial
    run leaves no report and is redone rather than silently compared.
    """
    report = out_prefix.with_suffix(".json")
    if report.exists() and not force:
        print(f"  [reuse] {report} already exists")
        return

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    command = [
        str(python), "-m", "genotools",
        "--pfile", str(geno),
        "--out", str(out_prefix),
        "--ancestry",
        "--all-sample",
        "--full-output",
    ]
    # --ref-panel/--ref-labels are required even with --model: inference
    # re-derives the reference PCA from the panel (round 12). --model is
    # additive, not an alternative.
    command += ["--ref-panel", str(ref_panel), "--ref-labels", str(ref_labels)]
    if model is not None:
        command += ["--model", str(model)]

    print(f"  [run]   {' '.join(command)}")
    # cwd=REPO_ROOT so both environments import the same working tree: the
    # libraries are the variable, the GenoTools source is not.
    result = subprocess.run(command, cwd=REPO_ROOT)
    if result.returncode != 0:
        raise SystemExit(f"pipeline failed under {python} (exit {result.returncode})")


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--baseline-python", required=True, type=Path,
                        help="Interpreter with the pinned umap (e.g. .venv/bin/python)")
    parser.add_argument("--candidate-python", required=True, type=Path,
                        help="Interpreter with the new umap (e.g. .venv-next/bin/python)")
    parser.add_argument("--geno", required=True, type=Path, help="Cohort pfile prefix")
    parser.add_argument("--ref-panel", required=True, type=Path,
                        help="Reference panel bfile prefix (needed in both modes)")
    parser.add_argument("--ref-labels", required=True, type=Path,
                        help="Reference panel labels TSV (needed in both modes)")
    parser.add_argument("--model", type=Path, default=None,
                        help="Pre-trained model, required by --mode reuse-model")
    parser.add_argument("--work", required=True, type=Path,
                        help="Directory for both runs' outputs")
    parser.add_argument("--mode", choices=("retrain", "reuse-model"), default="retrain")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even where a completed run exists")
    parser.add_argument("--skip-genotypes", action="store_true",
                        help="Passed through to the comparator (much faster)")
    args = parser.parse_args()

    # Both modes need the panel: prediction re-derives the reference PCA from
    # it even when the classifier itself is loaded from --model.
    if not (args.ref_panel and args.ref_labels):
        parser.error("--ref-panel and --ref-labels are required in both modes")
    if args.mode == "reuse-model" and not args.model:
        parser.error("--mode reuse-model needs --model")

    base_versions = environment_versions(args.baseline_python)
    cand_versions = environment_versions(args.candidate_python)

    print("=" * 72)
    print(f"umap comparison — mode: {args.mode}")
    print("=" * 72)
    delta = describe_difference(base_versions, cand_versions)
    if delta:
        print("Environment differences (baseline -> candidate):")
        print("\n".join(delta))
    else:
        print("WARNING: the two environments are identical; this compares nothing.")
    print()

    baseline_out = args.work / "baseline" / "out"
    candidate_out = args.work / "candidate" / "out"

    model = args.model if args.mode == "reuse-model" else None
    print("Baseline run:")
    run_pipeline(args.baseline_python, baseline_out, args.geno,
                 args.ref_panel, args.ref_labels, model, args.force)
    print("Candidate run:")
    run_pipeline(args.candidate_python, candidate_out, args.geno,
                 args.ref_panel, args.ref_labels, model, args.force)

    # Reuse the existing report comparator rather than reimplementing it: both
    # runs emit the same JSON schema, so "old vs new CLI" and "old vs new umap"
    # are the same comparison.
    command = [
        sys.executable, str(COMPARATOR),
        "--old", str(baseline_out),
        "--new", str(candidate_out),
    ]
    if args.skip_genotypes:
        command.append("--skip-genotypes")
    print("\n" + "=" * 72)
    print("Comparing reports")
    print("=" * 72)
    comparison = subprocess.run(command, cwd=REPO_ROOT)

    versions_path = args.work / "environments.json"
    versions_path.parent.mkdir(parents=True, exist_ok=True)
    versions_path.write_text(
        json.dumps({"baseline": base_versions, "candidate": cand_versions}, indent=2)
    )
    print(f"\nEnvironments recorded in {versions_path}")

    # A non-zero comparator exit means the calls moved. For this script that is
    # a finding to read, not a failure to fix, so it is reported rather than
    # propagated as a crash.
    if comparison.returncode != 0:
        print(
            "\nThe two environments do NOT agree. That is the number this "
            "script exists to produce - read the per-check detail above and "
            "record it in MIGRATION_2.0.md."
        )
    return comparison.returncode


if __name__ == "__main__":
    raise SystemExit(main())
