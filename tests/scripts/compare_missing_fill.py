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

"""Quantify what changing the missing-SNP fill does to ancestry calls.

Third in the series after ``compare_umap_versions.py`` (vary the libraries) and
``compare_snp_matching.py`` (vary the source tree). This one needs neither: both
fill strategies are reachable from one tree through
``--ancestry-missing-fill``, so it holds *everything* constant -- one process,
one loaded model, one reference PCA -- and varies only the value written into
the columns the cohort could not supply. Every moved call is therefore
attributable to the fill and to nothing else.

It also doubles as the first thing to run against a cohort whose prediction
looks wrong. The header it prints -- how much of the model's SNP list the
cohort carried -- is the number that decides whether there is anything else
worth investigating:

    matched 43,173 / 43,173 (100.00%)   the fill cannot matter; look elsewhere
    matched   412 / 43,173 (  0.95%)    nothing else matters until this is fixed

Usage
-----
    python tests/scripts/compare_missing_fill.py \\
        --geno ~/parity_data/GP2_r12_subset10k \\
        --ref-panel ~/.genotools/ref/ref_panel/ref_panel_gp2_prune_rm_underperform_pos_update \\
        --ref-labels ~/.genotools/ref/ref_panel/ref_panel_ancestry_updated.txt \\
        --model ~/parity_data/models/new_ancestry_ancestry_model \\
        --work /tmp/missing-fill-compare

Requires a pre-trained model: the point is to hold the model fixed. Run from the
repo root, which puts the working tree first on ``sys.path``.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Any, Dict

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]

STRATEGIES = ("constant", "ref-mean")


def _predict(
    strategy: str, args: argparse.Namespace, model: Any
) -> Dict[str, Any]:
    """Build the feature matrix under one fill strategy and predict from it."""
    from genotools.ancestry.diagnostics import AncestryDiagnostics, format_pc_drift
    from genotools.ancestry.preprocessing import get_raw_files

    work = args.work / strategy
    work.mkdir(parents=True, exist_ok=True)

    snps_file = work / "model.common_snps"
    snps_file.write_text("\n".join(model.common_snps) + "\n")

    diagnostics = AncestryDiagnostics()
    raw = get_raw_files(
        geno_path=str(args.geno),
        ref_panel=str(args.ref_panel),
        ref_labels=str(args.ref_labels),
        out_path=str(work / "out"),
        train=False,
        common_snps_file=str(snps_file),
        fill_strategy=strategy,
        # The comparison is the point; a refusal here would prevent it. The
        # overlap is printed either way, which is what a refusal is for.
        max_missing_fraction=1.0,
        diagnostics=diagnostics,
    )
    geno = raw["raw_geno"]
    predictions = model.predict(
        geno.drop(columns=["label"]),
        geno[["FID", "IID"]],
        out_path=work / "predict",
        diagnostics=diagnostics,
        diagnostics_prefix=str(work / "predict"),
    )

    print(f"\n--- fill = {strategy} " + "-" * (56 - len(strategy)))
    print(f"  {diagnostics.snp_overlap.format_summary()}")
    print(format_pc_drift(diagnostics.pc_drift))
    counts = predictions.predictions["predicted_ancestry"].value_counts()
    print(f"  labels: {counts.to_dict()}")
    if diagnostics.admixture:
        print(
            f"  CAH: {diagnostics.admixture['n_cah']} of "
            f"{diagnostics.admixture['n_samples']}"
        )
    return {
        "diagnostics": diagnostics,
        "labels": predictions.predictions.set_index("IID")["predicted_ancestry"],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--geno", required=True, type=Path, help="Cohort pfile prefix")
    parser.add_argument("--ref-panel", required=True, type=Path)
    parser.add_argument("--ref-labels", required=True, type=Path)
    parser.add_argument("--model", required=True, type=Path,
                        help="Pre-trained model directory (held constant across arms)")
    parser.add_argument("--work", required=True, type=Path,
                        help="Directory for both arms' intermediates")
    args = parser.parse_args()
    args.geno = args.geno.expanduser()
    args.ref_panel = args.ref_panel.expanduser()
    args.ref_labels = args.ref_labels.expanduser()
    args.model = args.model.expanduser()
    args.work = args.work.expanduser()

    from genotools.ancestry import AncestryModel
    from genotools.core.logging import setup_logging

    setup_logging(level="INFO")
    model = AncestryModel.load(args.model)
    print(f"model: {args.model} ({len(model.common_snps)} common SNPs)")
    print(f"cohort: {args.geno}")

    arms = {strategy: _predict(strategy, args, model) for strategy in STRATEGIES}

    baseline, candidate = (arms[s]["labels"] for s in STRATEGIES)
    shared = baseline.index.intersection(candidate.index)
    baseline, candidate = baseline.loc[shared], candidate.loc[shared]
    moved = baseline != candidate

    print("\n=== fill: constant -> ref-mean " + "=" * 45)
    print(f"  {int(moved.sum())} of {len(shared)} calls moved "
          f"({100 * moved.mean():.3f}%)")
    if moved.any():
        transitions = (
            pd.DataFrame({"old": baseline[moved], "new": candidate[moved]})
            .value_counts()
            .sort_values(ascending=False)
        )
        for (old, new), count in transitions.items():
            print(f"    {old:>6} -> {new:<6} {count}")

    overlaps = {s: arms[s]["diagnostics"].snp_overlap.n_filled for s in STRATEGIES}
    assert len(set(overlaps.values())) == 1, (
        f"the arms disagree about how many SNPs were filled ({overlaps}), so "
        f"they are not comparing the same match"
    )
    if overlaps["ref-mean"] == 0:
        print("\n  Nothing was filled, so this cohort cannot be affected by the")
        print("  change. That is the expected result for a cohort that carries")
        print("  the model's SNP list in full.")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(REPO_ROOT))
    raise SystemExit(main())
