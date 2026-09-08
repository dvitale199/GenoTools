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

"""Refit an ancestry model on the reference panel, for a SNP list already fixed.

A normal `--ancestry` training run derives its common-SNP list from the
intersection of the panel and the cohort, so retraining means re-reading the
cohort -- at GP2 release scale, a dense 8-byte matrix peaking near 187 GiB
(REFACTOR item 40). But the *model* is fitted on the panel's labels; the cohort
only ever enters through that intersection. So when the SNP list is already
known and is not what you want to change, the cohort is not needed at all: the
panel restricted to those SNPs is the entire training input.

That is what this does, and the case it exists for is a library upgrade. A model
records the versions it was fitted under, and ancestry calls move with them
(~1.2% for the umap unpin alone, round 16). A model fitted under versions the
package no longer permits is a model nobody can reproduce -- which is what
happened to the round-19 GP2 model, fitted under umap-learn 0.5.3 while
`setup.py` requires >=0.5.5, because the development venv had drifted from
`requirements-lock.txt`.

What this deliberately does NOT do:

  - re-derive the SNP list. It reuses the given one, so the new model is
    comparable to the old, variant for variant. Re-deriving needs the cohort.
  - predict anything. No cohort, so no predictions and no per-cohort
    diagnostics. Validate the result against a cohort separately.

The new model is not bit-identical to the old one and is not meant to be: the
grid's top is a plateau (48 of 216 candidates within one fold-std, REFACTOR
item 41), so selection among near-equals moves with any perturbation. Expect
the substantive parameters to reproduce and the UMAP shape pair to wander.

Usage:

    python tests/scripts/retrain_reference_model.py \\
        --ref-panel ~/.genotools/ref/ref_panel/<panel prefix> \\
        --ref-labels ~/.genotools/ref/ref_panel/ref_panel_ancestry_updated.txt \\
        --snplist <old model>/common_snps.txt \\
        --out ~/retrained/nba_gp2_r12

Compare against the model it replaces with:

    python tests/scripts/check_model_health.py <out dir>
"""

import argparse
import json
import logging
import sys
import time
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from genotools.ancestry import AncestryModel  # noqa: E402
from genotools.ancestry.config import AncestryConfig  # noqa: E402
from genotools.core.executors import run_command  # noqa: E402
from genotools.dependencies import check_plink2  # noqa: E402

logger = logging.getLogger("genotools.retrain")


def build_reference_matrix(
    ref_panel: Path, ref_labels: Path, snplist: Path, workdir: Path
) -> pd.DataFrame:
    """The panel, restricted to `snplist`, as a labeled feature frame.

    This mirrors the reference half of `ancestry/preprocessing.get_raw_files`
    (the `--extract`, the `--recode A`, the six-column drop, the `_`-suffix
    strip, and the fam/label merge). It is a copy rather than a call because
    that function derives the SNP list from a cohort, which is the one thing
    this script exists to avoid. Keep the two in step.

    Returns:
        Frame of `FID`, `IID`, one column per SNP, and `label`.
    """
    plink2 = check_plink2()
    prefix = workdir / "ref_common"

    run_command(
        f"{plink2} --bfile {ref_panel} --extract {snplist} "
        f"--make-bed --out {prefix}",
        tool_name="plink2",
    )
    run_command(f"{plink2} --bfile {prefix} --recode A --out {prefix}", tool_name="plink2")

    raw = pd.read_csv(f"{prefix}.raw", sep=r"\s+")
    ids = raw[["FID", "IID"]]
    snps = raw.drop(columns=["FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"])
    snp_cols = snps.columns.str.extract("(.*)_")[0]
    snps.columns = snp_cols
    ref_raw = pd.concat([ids, snps], axis=1)
    ref_raw.columns = ["FID", "IID"] + list(snp_cols)

    ancestry = pd.read_csv(
        ref_labels, sep="\t", header=None, names=["FID", "IID", "label"]
    )
    ref_fam = pd.read_csv(f"{ref_panel}.fam", sep=r"\s+", header=None)
    ref_labeled = ref_fam.merge(
        ancestry, how="left", left_on=[0, 1], right_on=["FID", "IID"]
    )
    labeled = ref_raw.merge(ref_labeled, how="left", on=["FID", "IID"])
    labeled.drop(columns=[0, 1, 2, 3, 4, 5], inplace=True)
    return labeled


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ref-panel", type=Path, required=True, help="Panel bfile prefix")
    parser.add_argument("--ref-labels", type=Path, required=True, help="TSV: FID IID label")
    parser.add_argument("--snplist", type=Path, required=True, help="One rsID per line")
    parser.add_argument("--out", type=Path, required=True, help="Model directory to write")
    parser.add_argument(
        "--workdir", type=Path, default=None, help="Intermediates (default: <out>_work)"
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s"
    )

    workdir = args.workdir or Path(f"{args.out}_work")
    workdir.mkdir(parents=True, exist_ok=True)

    started = time.time()
    labeled = build_reference_matrix(
        args.ref_panel, args.ref_labels, args.snplist, workdir
    )

    labels = pd.Series(
        labeled["label"].values, index=labeled["IID"].values, name="label"
    )
    ref_data = labeled.drop(columns=["label"])
    snp_columns = [c for c in ref_data.columns if c not in ("FID", "IID")]

    logger.info(
        f"Training matrix: {len(ref_data)} panel samples x {len(snp_columns)} SNPs"
    )
    logger.info(f"Label counts:\n{labels.value_counts()}")
    unlabeled = int(labels.isna().sum())
    if unlabeled:
        logger.warning(f"{unlabeled} panel sample(s) carry no label")

    # Exactly the runner's config for a default run: both training flags
    # default to None/3, which is what TrainingConfig already holds.
    model = AncestryModel(config=AncestryConfig())
    model.fit(ref_data, labels, out_path=args.out.parent)
    model.common_snps = list(snp_columns)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    model.save(args.out)

    cv_results = getattr(model, "_cv_results", None)
    if cv_results is not None:
        grid_path = Path(f"{args.out}_grid_search.txt")
        cv_results.to_csv(grid_path, sep="\t", index=False)
        logger.info(f"Grid written to {grid_path}")

    elapsed = time.time() - started
    metadata = json.loads((args.out / "metadata.json").read_text())
    logger.info(f"Saved model to {args.out} in {elapsed / 60:.1f} min")
    logger.info(f"best_params: {json.dumps(metadata['best_params'])}")
    logger.info(
        f"train/test accuracy: "
        f"{metadata['training_metrics']['train_accuracy']:.6f} / "
        f"{metadata['training_metrics']['test_accuracy']:.6f}"
    )
    logger.info(f"versions: {json.dumps(metadata['versions'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
