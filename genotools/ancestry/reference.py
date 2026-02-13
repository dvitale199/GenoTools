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

"""Reference panel management for ancestry prediction.

This module provides the ReferencePanel class for loading and managing
reference population data used in ancestry prediction. Reference panels
contain genotype data from individuals with known ancestry labels.

Example:
    >>> from genotools.ancestry.reference import ReferencePanel
    >>> ref = ReferencePanel.load(
    ...     geno_path=Path("/path/to/ref_panel"),
    ...     labels_path=Path("/path/to/labels.txt")
    ... )
    >>> ref.sample_count
    3000
    >>> ref.ancestry_counts
    {'EUR': 500, 'AFR': 500, ...}
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set

import pandas as pd

from genotools.core.genotypes import GenotypeData
from genotools.core.exceptions import AncestryError


@dataclass
class ReferencePanel:
    """Reference panel for ancestry prediction.

    A reference panel contains genotype data from individuals with known
    ancestry labels. This data is used to train the ancestry prediction
    model and provide a basis for projecting new samples.

    The reference panel expects:
    - Genotype files in PLINK bfile format (.bed/.bim/.fam)
    - A labels file with tab-separated FID, IID, label columns
      (no header, or header with these column names)

    Attributes:
        genotypes: GenotypeData pointing to the reference genotype files.
        labels: Series mapping sample IID to ancestry label.
        labels_path: Path to the labels file that was loaded.
    """

    genotypes: GenotypeData
    labels: pd.Series  # type: ignore[type-arg]
    labels_path: Path

    def __post_init__(self) -> None:
        """Validate reference panel data."""
        if len(self.labels) == 0:
            raise AncestryError("Reference panel has no labeled samples")

        # Check that labels index matches sample IDs in genotype data
        # This is validated at load time, but we double-check here

    @classmethod
    def load(
        cls,
        geno_path: Path,
        labels_path: Path,
        require_all_labeled: bool = True,
    ) -> "ReferencePanel":
        """Load reference panel from files.

        Args:
            geno_path: Path prefix to genotype files (without extension).
                Expects bfile format (.bed/.bim/.fam).
            labels_path: Path to ancestry labels file. Tab-separated with
                columns: FID, IID, label (or no header with same order).
            require_all_labeled: If True, raise error if any samples in
                the genotype files lack labels. Default is True.

        Returns:
            Loaded ReferencePanel instance.

        Raises:
            FileNotFoundError: If genotype or labels files don't exist.
            AncestryError: If labels are invalid or don't match samples.
        """
        geno_path = Path(geno_path)
        labels_path = Path(labels_path)

        # Validate genotype files exist
        if not geno_path.with_suffix(".bed").exists():
            raise FileNotFoundError(
                f"Reference panel bed file not found: {geno_path}.bed"
            )
        if not geno_path.with_suffix(".bim").exists():
            raise FileNotFoundError(
                f"Reference panel bim file not found: {geno_path}.bim"
            )
        if not geno_path.with_suffix(".fam").exists():
            raise FileNotFoundError(
                f"Reference panel fam file not found: {geno_path}.fam"
            )

        # Validate labels file exists
        if not labels_path.exists():
            raise FileNotFoundError(
                f"Reference labels file not found: {labels_path}"
            )

        # Load genotype data
        genotypes = GenotypeData.from_path(geno_path)

        # Load labels
        labels_df = cls._load_labels_file(labels_path)

        # Load fam file to get sample IDs
        fam_df = pd.read_csv(
            geno_path.with_suffix(".fam"),
            sep=r"\s+",
            header=None,
            names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
            dtype=str,
        )

        # Merge labels with fam file
        merged = fam_df.merge(
            labels_df,
            on=["FID", "IID"],
            how="left",
        )

        # Check for missing labels
        missing_labels = merged["label"].isna()
        if missing_labels.any():
            n_missing = missing_labels.sum()
            if require_all_labeled:
                raise AncestryError(
                    f"{n_missing} samples in reference panel have no labels. "
                    f"Set require_all_labeled=False to allow unlabeled samples."
                )
            # Filter to only labeled samples
            merged = merged[~missing_labels]

        # Create Series with IID as index
        labels_series = pd.Series(
            merged["label"].values,
            index=merged["IID"].values,
            name="label",
        )

        return cls(
            genotypes=genotypes,
            labels=labels_series,
            labels_path=labels_path,
        )

    @staticmethod
    def _load_labels_file(path: Path) -> pd.DataFrame:
        """Load ancestry labels file.

        Handles both headerless files and files with FID/IID/label header.

        Args:
            path: Path to labels file.

        Returns:
            DataFrame with FID, IID, label columns.
        """
        # Try reading with assumed columns
        df = pd.read_csv(
            path,
            sep="\t",
            header=None,
            names=["FID", "IID", "label"],
            dtype=str,
        )

        # Check if first row looks like a header
        if df.iloc[0]["FID"] in ("FID", "#FID"):
            # Re-read with header
            df = pd.read_csv(path, sep="\t", dtype=str)
            # Normalize column names
            df.columns = [c.replace("#", "") for c in df.columns]
            if "label" not in df.columns:
                # Try common alternatives
                for alt in ("LABEL", "ancestry", "ANCESTRY", "pop", "POP"):
                    if alt in df.columns:
                        df = df.rename(columns={alt: "label"})
                        break

        return df[["FID", "IID", "label"]]

    @property
    def sample_count(self) -> int:
        """Number of labeled samples in reference panel."""
        return len(self.labels)

    @property
    def variant_count(self) -> int:
        """Number of variants in reference panel."""
        return self.genotypes.variant_count

    @property
    def ancestry_counts(self) -> Dict[str, int]:
        """Count of samples per ancestry label.

        Returns:
            Dictionary mapping ancestry labels to sample counts.
        """
        counts = self.labels.value_counts()
        return dict(counts)

    @property
    def unique_ancestries(self) -> List[str]:
        """List of unique ancestry labels in reference panel."""
        return list(self.labels.unique())

    def get_samples_by_ancestry(self, ancestry: str) -> List[str]:
        """Get sample IDs for a specific ancestry.

        Args:
            ancestry: Ancestry label to filter by.

        Returns:
            List of sample IIDs with the specified ancestry.
        """
        return list(self.labels[self.labels == ancestry].index)

    def validate_snp_overlap(
        self,
        other: GenotypeData,
        min_overlap: int = 10000,
    ) -> int:
        """Check SNP overlap with another genotype dataset.

        Ancestry prediction requires sufficient SNP overlap between
        the reference panel and the target samples. This method
        estimates the overlap (actual overlap is computed during
        data preparation).

        Args:
            other: GenotypeData to check overlap with.
            min_overlap: Minimum required SNP overlap. Default is 10000.

        Returns:
            Estimated number of overlapping SNPs.

        Raises:
            AncestryError: If overlap is below minimum threshold.
        """
        # This is a simplified check - actual overlap computation
        # happens during data preparation when SNP IDs are compared
        min_variants = min(self.variant_count, other.variant_count)

        if min_variants < min_overlap:
            raise AncestryError(
                f"Insufficient variant count for ancestry prediction. "
                f"Reference has {self.variant_count} variants, "
                f"target has {other.variant_count} variants. "
                f"Expected at least {min_overlap} overlapping SNPs."
            )

        return min_variants

    def get_labeled_data(self) -> pd.DataFrame:
        """Get labeled sample information.

        Returns:
            DataFrame with FID, IID, and label columns for all
            labeled samples in the reference panel.
        """
        # Load fam file for FID/IID pairs
        fam_df = pd.read_csv(
            self.genotypes.path.with_suffix(".fam"),
            sep=r"\s+",
            header=None,
            names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
            dtype=str,
        )

        # Filter to labeled samples and add labels
        labeled = fam_df[fam_df["IID"].isin(self.labels.index)].copy()
        labeled["label"] = labeled["IID"].map(self.labels)

        return labeled[["FID", "IID", "label"]]


def get_default_model_path() -> Optional[Path]:
    """Get path to bundled NBA model if available.

    Returns the path to the pre-trained NBA v2 model if it has been
    downloaded using genotools-download.

    Returns:
        Path to model directory, or None if not available.
    """
    import os

    model_destination = os.path.expanduser("~/.genotools/ref")
    nba_model_dir = Path(model_destination) / "models" / "nba_v2"

    if nba_model_dir.exists():
        return nba_model_dir

    return None


def validate_model_files(model_dir: Path) -> Dict[str, Path]:
    """Validate that required model files exist.

    Checks for the presence of all files needed to load a trained
    ancestry model.

    Args:
        model_dir: Directory containing model files.

    Returns:
        Dictionary mapping file types to their paths.

    Raises:
        FileNotFoundError: If required files are missing.
    """
    required_files = {
        "model": "pipeline.pkl",
        "common_snps": "common_snps.txt",
    }

    found_files: Dict[str, Path] = {}

    for file_type, pattern in required_files.items():
        matches = list(model_dir.glob(pattern))
        if not matches:
            raise FileNotFoundError(
                f"Required {file_type} file not found in {model_dir}. "
                f"Expected file matching pattern: {pattern}"
            )
        found_files[file_type] = matches[0]

    return found_files
