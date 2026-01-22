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

"""QC step configuration dataclasses.

Each QC step has a corresponding frozen configuration dataclass that
validates parameters on construction. All configs inherit validation
helpers from ThresholdConfig.
"""

from dataclasses import dataclass

from genotools.core.config import ThresholdConfig


@dataclass(frozen=True)
class CallrateConfig(ThresholdConfig):
    """Configuration for sample call rate filtering.

    Removes samples with missing genotype rate exceeding the threshold.
    Uses PLINK2 --mind flag.

    Attributes:
        mind: Maximum allowable proportion of missing genotypes per sample.
            Default is 0.05 (5% missing allowed).
    """

    mind: float = 0.05

    def __post_init__(self) -> None:
        self._validate_proportion(self.mind, "mind")


@dataclass(frozen=True)
class SexConfig(ThresholdConfig):
    """Configuration for sex check filtering.

    Removes samples with sex discrepancies based on X chromosome
    F-statistic values. Uses PLINK 1.9 --check-sex flag.

    Attributes:
        female_max_f: F-statistic upper bound for females.
            Samples with reported sex female and F > this value are flagged.
            Default is 0.25.
        male_min_f: F-statistic lower bound for males.
            Samples with reported sex male and F < this value are flagged.
            Default is 0.75.
    """

    female_max_f: float = 0.25
    male_min_f: float = 0.75

    def __post_init__(self) -> None:
        self._validate_proportion(self.female_max_f, "female_max_f")
        self._validate_proportion(self.male_min_f, "male_min_f")
        if self.female_max_f >= self.male_min_f:
            raise ValueError(
                f"female_max_f ({self.female_max_f}) must be less than "
                f"male_min_f ({self.male_min_f})"
            )


@dataclass(frozen=True)
class HetConfig(ThresholdConfig):
    """Configuration for heterozygosity filtering.

    Removes samples with extreme heterozygosity rates, either using
    fixed F-statistic bounds or automatic detection based on standard
    deviations from the mean.

    Attributes:
        f_lower: Lower F-statistic bound. Samples with F < this are removed.
            Set to -1.0 along with f_upper to enable auto_detect mode.
            Default is -0.15.
        f_upper: Upper F-statistic bound. Samples with F > this are removed.
            Set to -1.0 along with f_lower to enable auto_detect mode.
            Default is 0.15.
        auto_detect: If True (or if both f_lower and f_upper are -1.0),
            use 3 standard deviations from the mean heterozygosity rate
            to determine outliers. Default is False.
    """

    f_lower: float = -0.15
    f_upper: float = 0.15
    auto_detect: bool = False

    def __post_init__(self) -> None:
        # Special case: [-1, -1] means auto-detect mode
        if self.f_lower == -1.0 and self.f_upper == -1.0:
            return
        if not self.auto_detect:
            if self.f_lower >= self.f_upper:
                raise ValueError(
                    f"f_lower ({self.f_lower}) must be less than "
                    f"f_upper ({self.f_upper})"
                )


@dataclass(frozen=True)
class RelatedConfig(ThresholdConfig):
    """Configuration for relatedness filtering.

    Removes related and/or duplicated samples based on KING kinship
    coefficients. Uses PLINK2 --make-king-table and --king-cutoff.

    Kinship coefficient reference (Manichaikul et al. 2010):
        - >= 0.354: MZ twin/duplicate
        - 0.177-0.354: 1st degree (parent-child, full sibling)
        - 0.0884-0.177: 2nd degree (half-sibling, avuncular, grandparent)
        - 0.0442-0.0884: 3rd degree (first cousin)

    Attributes:
        related_cutoff: Kinship coefficient threshold for relatedness.
            Pairs above this are considered related. Default is 0.0884
            (2nd degree relatives).
        duplicated_cutoff: Kinship coefficient threshold for duplicates.
            Pairs above this are considered duplicates. Default is 0.354.
        prune_related: Whether to remove related samples. Default is True.
        prune_duplicated: Whether to remove duplicate samples. Default is True.
    """

    related_cutoff: float = 0.0884
    duplicated_cutoff: float = 0.354
    prune_related: bool = True
    prune_duplicated: bool = True

    def __post_init__(self) -> None:
        self._validate_proportion(self.related_cutoff, "related_cutoff", 0.0, 0.5)
        self._validate_proportion(self.duplicated_cutoff, "duplicated_cutoff", 0.0, 0.5)
        if self.prune_related and not self.prune_duplicated:
            raise ValueError(
                "Cannot prune related samples without also pruning duplicates. "
                "Set prune_duplicated=True or prune_related=False."
            )


@dataclass(frozen=True)
class GenoConfig(ThresholdConfig):
    """Configuration for variant missingness filtering.

    Removes variants with high missing genotype rates.
    Uses PLINK2 --geno flag.

    Attributes:
        geno: Maximum allowable proportion of missing genotypes per variant.
            Default is 0.05 (5% missing allowed).
    """

    geno: float = 0.05

    def __post_init__(self) -> None:
        self._validate_proportion(self.geno, "geno")


@dataclass(frozen=True)
class CaseControlConfig(ThresholdConfig):
    """Configuration for case-control differential missingness filtering.

    Removes variants with significantly different missing rates between
    cases and controls. Uses PLINK 1.9 --test-missing.

    Attributes:
        p_threshold: P-value threshold for differential missingness test.
            Variants with P <= threshold are removed. Default is 1e-4.
    """

    p_threshold: float = 1e-4

    def __post_init__(self) -> None:
        self._validate_proportion(self.p_threshold, "p_threshold")


@dataclass(frozen=True)
class HaplotypeConfig(ThresholdConfig):
    """Configuration for haplotype missingness filtering.

    Removes variants with non-random patterns of missing data that
    suggest genotyping problems. Uses PLINK 1.9 --test-mishap.

    Attributes:
        p_threshold: P-value threshold for haplotype missingness test.
            Variants with P <= threshold are removed. Default is 1e-4.
        adjust_for_sample_size: If True and sample size > 10,000, use
            Bonferroni-corrected threshold (0.05/n). Default is True.
    """

    p_threshold: float = 1e-4
    adjust_for_sample_size: bool = True

    def __post_init__(self) -> None:
        self._validate_proportion(self.p_threshold, "p_threshold")


@dataclass(frozen=True)
class HWEConfig(ThresholdConfig):
    """Configuration for Hardy-Weinberg equilibrium filtering.

    Removes variants violating Hardy-Weinberg equilibrium, which may
    indicate genotyping errors. Uses PLINK 1.9 --hwe.

    Attributes:
        hwe_threshold: P-value threshold for HWE test.
            Variants with P <= threshold are removed. Default is 1e-4.
        filter_controls: If True, only test HWE in control samples.
            Useful for case-control studies where cases may have
            true deviations from HWE. Default is False.
    """

    hwe_threshold: float = 1e-4
    filter_controls: bool = False

    def __post_init__(self) -> None:
        self._validate_proportion(self.hwe_threshold, "hwe_threshold")


@dataclass(frozen=True)
class LDConfig(ThresholdConfig):
    """Configuration for linkage disequilibrium pruning.

    Removes variants in high LD to create a set of approximately
    independent variants. Uses PLINK2 --indep-pairwise.

    Attributes:
        window_size: Size of sliding window in variant count.
            Default is 50.
        step_size: Step size for sliding window. Default is 5.
        r2_threshold: r² threshold. Variant pairs with r² > threshold
            have one variant removed. Default is 0.5.
    """

    window_size: int = 50
    step_size: int = 5
    r2_threshold: float = 0.5

    def __post_init__(self) -> None:
        self._validate_positive(float(self.window_size), "window_size")
        self._validate_positive(float(self.step_size), "step_size")
        self._validate_proportion(self.r2_threshold, "r2_threshold")
