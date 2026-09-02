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

"""Quantify what fixing ``get_common_snps`` does to ancestry calls.

Sibling of ``compare_umap_versions.py``, asking the same question of a
different variable: that script holds the source constant and varies the
installed libraries, this one holds the libraries constant and varies the
source. Both delegate to ``compare_ancestry_run.py``, so the movement table
means the same thing in both.

The change under measurement (REFACTOR.md item 13) has two halves:

palindromic sites are excluded
    Unreachable on a reference panel built per ``docs/prep_reference_panel.md``,
    which already drops them — the GP2 panel in use has 0 of 209,517, so no
    palindromic site can survive the intersection and this half is expected to
    move nothing. Run it against a panel that kept them and it will.

position ties are broken by call rate
    The live half. A cohort carries the same site under several probe IDs, and
    the previous ``drop_duplicates`` kept whichever row the merge emitted
    first, which is arbitrary and sensitive to input variant order.

Why ``--mode reuse-model`` is the right measurement here
-------------------------------------------------------
``get_common_snps`` is called twice, with its arguments **reversed**. The first
call extracts from the reference panel and decides the model's common-SNP list;
the second extracts from the cohort and decides which probe represents each
position in the feature matrix. Only the second call's output changes (the
first is byte-identical, since its candidate rows all name the same reference
variant), so the fitted model and the reference PCA are untouched and the whole
effect lands on prediction.

That makes ``--mode reuse-model`` exact rather than merely cheaper: both arms
load the same model, so every moved call is attributable to the changed probe
selection. ``--mode retrain`` would additionally fold in UMAP's own fit
stochasticity, which is noise for this question — round 16 measured that noise
separately and it is larger than the signal here.

Usage
-----
    python tests/scripts/compare_snp_matching.py \\
        --python .venv/bin/python \\
        --baseline-repo /path/to/worktree/at/main \\
        --geno ~/parity_data/GP2_r12_subset10k \\
        --ref-panel ~/.genotools/ref/ref_panel/ref_panel_gp2_prune_rm_underperform_pos_update \\
        --ref-labels ~/.genotools/ref/ref_panel/ref_panel_ancestry_updated.txt \\
        --model ~/parity_data/models/new_ancestry_ancestry_model \\
        --work /tmp/snp-matching-compare

Create the baseline tree with ``git worktree add <dir> <commit>``; it must be a
checkout of the code *before* the fix, and the candidate defaults to this
working tree. Completed runs are reused unless ``--force`` is passed.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from compare_umap_versions import COMPARATOR, REPO_ROOT, run_pipeline

PALINDROMIC = {"AT", "TA", "CG", "GC"}


def describe_panel(ref_panel: Path) -> dict:
    """Report the panel's exposure to each half of the fix.

    Printed up front because it decides what the run can possibly show: a
    palindrome-free, position-unique panel cannot exercise either half on the
    reference side, and reading a 0.00% movement without that context invites
    the wrong conclusion ("the fix does nothing") instead of the right one
    ("this panel was never exposed").
    """
    import pandas as pd

    bim = pd.read_csv(
        f"{ref_panel}.bim", sep="\t", header=None, dtype={0: str}, low_memory=False
    )
    bim.columns = ["chr", "rsid", "kb", "pos", "a1", "a2"]
    return {
        "variants": len(bim),
        "palindromic": int((bim["a1"] + bim["a2"]).isin(PALINDROMIC).sum()),
        "duplicate_positions": int(bim.duplicated(subset=["chr", "pos"]).sum()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--python", required=True, type=Path,
                        help="One interpreter, used for both arms")
    parser.add_argument("--baseline-repo", required=True, type=Path,
                        help="Source tree WITHOUT the fix (git worktree at the prior commit)")
    parser.add_argument("--candidate-repo", type=Path, default=REPO_ROOT,
                        help="Source tree WITH the fix (default: this working tree)")
    parser.add_argument("--geno", required=True, type=Path, help="Cohort pfile prefix")
    parser.add_argument("--ref-panel", required=True, type=Path)
    parser.add_argument("--ref-labels", required=True, type=Path)
    parser.add_argument("--model", type=Path, default=None,
                        help="Pre-trained model, required by --mode reuse-model")
    parser.add_argument("--work", required=True, type=Path,
                        help="Directory for both runs' outputs")
    parser.add_argument("--mode", choices=("reuse-model", "retrain"), default="reuse-model")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even where a completed run exists")
    parser.add_argument("--skip-genotypes", action="store_true",
                        help="Passed through to the comparator (much faster)")
    args = parser.parse_args()

    if args.mode == "reuse-model" and not args.model:
        parser.error("--mode reuse-model needs --model")
    for repo in (args.baseline_repo, args.candidate_repo):
        if not (repo / "genotools" / "ancestry" / "preprocessing.py").is_file():
            parser.error(f"{repo} is not a GenoTools source tree")
    if args.baseline_repo.resolve() == args.candidate_repo.resolve():
        parser.error("baseline and candidate repos are the same tree; this compares nothing")

    print("=" * 72)
    print(f"SNP-matching comparison — mode: {args.mode}")
    print("=" * 72)
    print(f"  baseline source : {args.baseline_repo}")
    print(f"  candidate source: {args.candidate_repo}")
    print(f"  interpreter     : {args.python}")

    panel = describe_panel(args.ref_panel)
    print(f"\nReference panel exposure ({panel['variants']} variants):")
    print(f"  palindromic sites   : {panel['palindromic']}")
    print(f"  duplicate positions : {panel['duplicate_positions']}")
    if not panel["palindromic"]:
        print("  -> panel is palindrome-free, so the palindrome half cannot move a call here")
    print()

    baseline_out = args.work / "baseline" / "out"
    candidate_out = args.work / "candidate" / "out"
    model = args.model if args.mode == "reuse-model" else None
    keep = not args.skip_genotypes

    print("Baseline run (without the fix):")
    run_pipeline(args.python, baseline_out, args.geno, args.ref_panel,
                 args.ref_labels, model, args.force, keep, cwd=args.baseline_repo)
    print("Candidate run (with the fix):")
    run_pipeline(args.python, candidate_out, args.geno, args.ref_panel,
                 args.ref_labels, model, args.force, keep, cwd=args.candidate_repo)

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

    provenance = args.work / "sources.json"
    provenance.parent.mkdir(parents=True, exist_ok=True)

    def head(repo: Path) -> str:
        result = subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                                capture_output=True, text=True)
        return result.stdout.strip() or "unknown"

    provenance.write_text(json.dumps({
        "mode": args.mode,
        "interpreter": str(args.python),
        "baseline": {"repo": str(args.baseline_repo), "head": head(args.baseline_repo)},
        "candidate": {"repo": str(args.candidate_repo), "head": head(args.candidate_repo)},
        "ref_panel": {"path": str(args.ref_panel), **panel},
    }, indent=2))
    print(f"\nSources recorded in {provenance}")

    # A non-zero comparator exit means calls moved, which is the number this
    # script exists to produce — reported, not propagated as a crash.
    if comparison.returncode != 0:
        print(
            "\nThe two source trees do NOT agree. Read the per-label movement "
            "table above and record it in REFACTOR.md item 13."
        )
    return comparison.returncode


if __name__ == "__main__":
    raise SystemExit(main())
