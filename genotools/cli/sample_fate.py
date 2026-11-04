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

"""The per-sample fate table, ``{out}_sample_fate.tsv``.

One row per input sample, saying what happened to it: which step removed it
(or ``kept``), the ancestry group it was QC'd in, and the per-sample statistics
the run kept. The assembly is pure pandas over frames the report already holds,
so it can be tested against hand-built ones without PLINK.

``removed_by`` is read from the report's own ``pruned_samples`` table, never
recomputed, so the two cannot disagree.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from ..core.logging import get_logger

logger = get_logger(__name__)

KEPT = "kept"
NO_ANCESTRY = "all"
STAT_COLUMNS = ["F_MISS", "sex_F", "het_F", "max_kinship"]
FATE_COLUMNS = ["FID", "IID", "ancestry", "removed_by"] + STAT_COLUMNS

# Subdirectories of {out}_plot_data/ that are not ancestry groups.
_NOT_A_GROUP = {"ancestry"}


def _ids(frame: pd.DataFrame, fid: str = "FID", iid: str = "IID") -> pd.DataFrame:
    """Copy of `frame` with PLINK's ``#FID``/``#IID`` headers and string IDs.

    A numeric ID column reads back as integers from one file and strings from
    another; comparing them as strings is what makes the joins line up.
    """
    frame = frame.rename(columns={"#FID": "FID", "#IID": "IID"})
    if fid != "FID" and fid in frame.columns:
        frame = frame.rename(columns={fid: "FID"})
    if iid != "IID" and iid in frame.columns:
        frame = frame.rename(columns={iid: "IID"})
    frame = frame.copy()
    for col in ("FID", "IID"):
        if col in frame.columns:
            frame[col] = frame[col].astype(str)
    return frame


def _join_keys(left: pd.DataFrame, right: pd.DataFrame) -> list:
    """FID+IID when both sides have a FID (a .psam may carry only IID)."""
    return ["FID", "IID"] if "FID" in left and "FID" in right else ["IID"]


def kinship_max(kinship: pd.DataFrame) -> pd.DataFrame:
    """Each sample's largest kinship over every pair it appears in.

    Args:
        kinship: KING pair table (``#FID1``, ``IID1``, ``FID2``, ``IID2``,
            ``KINSHIP``).

    Returns:
        ``FID``, ``IID``, ``max_kinship``, one row per sample in any pair.
    """
    empty = pd.DataFrame(columns=["FID", "IID", "max_kinship"])
    needed = {"IID1", "IID2", "KINSHIP"}
    kinship = kinship.rename(columns={"#IID1": "IID1"})
    if not needed.issubset(kinship.columns):
        return empty
    kinship = kinship.rename(columns={"#FID1": "FID1"})
    has_fid = {"FID1", "FID2"}.issubset(kinship.columns)
    sides = []
    for n in ("1", "2"):
        side = pd.DataFrame(
            {"IID": kinship[f"IID{n}"], "max_kinship": kinship["KINSHIP"]}
        )
        if has_fid:
            side.insert(0, "FID", kinship[f"FID{n}"])
        sides.append(side)
    stacked = _ids(pd.concat(sides, ignore_index=True))
    keys = ["FID", "IID"] if has_fid else ["IID"]
    out = stacked.groupby(keys, as_index=False, sort=False)["max_kinship"].max()
    return out if not out.empty else empty


def collect_sample_stats(tables: list[dict]) -> pd.DataFrame:
    """Per-sample statistics from the plot-data tables of every group.

    Args:
        tables: One dict per ancestry group, mapping a table's file stem
            (``callrate_sample_missingness``, ``sex_check``, ``het_f``,
            ``related_kinship``) to its frame. Missing tables are fine.

    Returns:
        ``FID``, ``IID`` and whichever of `STAT_COLUMNS` could be filled, one
        row per sample. A sample the step never saw has no value.
    """
    parts = {c: [] for c in STAT_COLUMNS}
    for group in tables:
        miss = group.get("callrate_sample_missingness")
        if miss is not None and "F_MISS" in miss.columns:
            parts["F_MISS"].append(_ids(miss)[["FID", "IID", "F_MISS"]])
        sex = group.get("sex_check")
        if sex is not None and {"F", "check"}.issubset(sex.columns):
            # The all-X pass: the filtered-X pass is a second opinion on the
            # same sample, not a second value.
            x = _ids(sex[sex["check"] == "all_x"])
            parts["sex_F"].append(x[["FID", "IID", "F"]].rename(columns={"F": "sex_F"}))
        het = group.get("het_f")
        if het is not None and "F" in het.columns:
            parts["het_F"].append(
                _ids(het)[["FID", "IID", "F"]].rename(columns={"F": "het_F"})
            )
        kin = group.get("related_kinship")
        if kin is not None:
            parts["max_kinship"].append(kinship_max(kin))

    merged: Optional[pd.DataFrame] = None
    for col in STAT_COLUMNS:
        if not parts[col]:
            continue
        stat = pd.concat(parts[col], ignore_index=True)
        keys = ["FID", "IID"] if "FID" in stat.columns else ["IID"]
        stat = stat.drop_duplicates(subset=keys, keep="first")
        if merged is None:
            merged = stat
        else:
            on = _join_keys(merged, stat)
            merged = merged.merge(stat, on=on, how="outer")
    if merged is None:
        return pd.DataFrame(columns=["FID", "IID"] + STAT_COLUMNS)
    for col in STAT_COLUMNS:
        if col not in merged.columns:
            merged[col] = float("nan")
    return merged


def assemble_sample_fate(
    input_samples: pd.DataFrame,
    pruned_samples: Optional[pd.DataFrame] = None,
    ancestry_labels: Optional[pd.DataFrame] = None,
    stats: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """One row per input sample: its group, the step that removed it, its stats.

    Args:
        input_samples: The input ``.psam`` table (``#FID``/``FID``, ``IID``).
        pruned_samples: The report's ``pruned_samples`` (``FID``, ``IID``,
            ``step``). A sample listed more than once reports its first step.
        ancestry_labels: The report's ``ancestry_labels`` (``FID``, ``IID``,
            ``label``), or None for a run without ``--ancestry``, where every
            sample's ancestry is ``all``.
        stats: Output of `collect_sample_stats`, or None.

    Returns:
        `FATE_COLUMNS`, in input order; every input sample exactly once
        (an ID repeated in the input keeps its first row).
    """
    fate = _ids(input_samples)
    if "IID" not in fate.columns:
        raise ValueError("input samples have no IID column")
    fate = fate.drop_duplicates(
        subset=["FID", "IID"] if "FID" in fate.columns else ["IID"], ignore_index=True
    )
    if "FID" not in fate.columns:
        fate.insert(0, "FID", "0")
    fate = fate[["FID", "IID"]]

    def attach(frame: pd.DataFrame, cols: list) -> pd.DataFrame:
        frame = _ids(frame)
        keys = _join_keys(fate, frame)
        frame = frame.drop_duplicates(subset=keys, keep="first")
        return fate.merge(frame[keys + cols], on=keys, how="left")

    if ancestry_labels is not None and len(ancestry_labels) > 0:
        fate["ancestry"] = attach(ancestry_labels, ["label"])["label"].astype(object)
    else:
        fate["ancestry"] = NO_ANCESTRY

    removed = pd.Series([None] * len(fate), dtype=object)
    if pruned_samples is not None and len(pruned_samples) > 0:
        removed = attach(pruned_samples, ["step"])["step"].astype(object)
        # A dropped sample that no prediction covers still has a group.
        if "label" in pruned_samples.columns and ancestry_labels is not None:
            lab = attach(pruned_samples, ["label"])["label"]
            fate["ancestry"] = fate["ancestry"].where(fate["ancestry"].notna(), lab)
    fate["removed_by"] = removed.where(removed.notna(), KEPT)

    if stats is not None and len(stats) > 0:
        fate = attach_stats(fate, stats)
    for col in STAT_COLUMNS:
        if col not in fate.columns:
            fate[col] = float("nan")
    return fate[FATE_COLUMNS]


def attach_stats(fate: pd.DataFrame, stats: pd.DataFrame) -> pd.DataFrame:
    """Left-join the statistics columns onto the fate table."""
    stats = _ids(stats)
    keys = _join_keys(fate, stats)
    cols = [c for c in STAT_COLUMNS if c in stats.columns]
    stats = stats.drop_duplicates(subset=keys, keep="first")
    return fate.merge(stats[keys + cols], on=keys, how="left")


def read_plot_tables(plot_data_dir: Optional[Path]) -> list[dict]:
    """Load the per-group tables the fate table reads from ``{out}_plot_data/``.

    A missing directory, group or file yields fewer tables, never an error:
    the statistics are a convenience and the fate table is the point.
    """
    if plot_data_dir is None or not Path(plot_data_dir).is_dir():
        return []
    stems = ("callrate_sample_missingness", "sex_check", "het_f", "related_kinship")
    groups = []
    for group_dir in sorted(Path(plot_data_dir).iterdir()):
        if not group_dir.is_dir() or group_dir.name in _NOT_A_GROUP:
            continue
        tables = {}
        for stem in stems:
            path = group_dir / f"{stem}.tsv"
            if not path.is_file():
                continue
            try:
                tables[stem] = pd.read_csv(path, sep="\t", dtype=str)
            except (OSError, ValueError, pd.errors.ParserError) as e:
                logger.warning(f"Sample fate: could not read {path}: {e}")
                continue
            for col in ("F_MISS", "F", "KINSHIP"):
                if col in tables[stem].columns:
                    tables[stem][col] = pd.to_numeric(tables[stem][col], errors="coerce")
        if tables:
            groups.append(tables)
    return groups


def write_sample_fate(
    out_path: Path,
    input_samples: Optional[pd.DataFrame],
    pruned_samples: Optional[pd.DataFrame],
    ancestry_labels: Optional[pd.DataFrame],
    plot_data_dir: Optional[Path],
) -> Optional[Path]:
    """Assemble and write ``{out_path}_sample_fate.tsv``.

    Returns:
        The path written, or None when there is nothing to write (no input
        sample table) or it could not be written. Never raises: the table is a
        report, and the run's real outputs are already on disk.
    """
    if input_samples is None or len(input_samples) == 0:
        return None
    path = Path(f"{out_path}_sample_fate.tsv")
    try:
        stats = collect_sample_stats(read_plot_tables(plot_data_dir))
        fate = assemble_sample_fate(
            input_samples, pruned_samples, ancestry_labels, stats
        )
        fate.to_csv(path, sep="\t", index=False, na_rep="NA")
    except (OSError, ValueError, KeyError) as e:
        logger.warning(f"Sample fate table not written: {e}")
        return None
    return path
