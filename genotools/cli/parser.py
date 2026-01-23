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

"""Command-line argument parsing for GenoTools pipeline.

This module provides typed argument parsing with validation, replacing the
legacy string-based boolean parsing with proper argparse actions.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Tuple, Sequence

from ..qc.config import (
    CallrateConfig,
    SexConfig,
    HetConfig,
    RelatedConfig,
    GenoConfig,
    CaseControlConfig,
    HaplotypeConfig,
    HWEConfig,
    LDConfig,
)
from ..ancestry.config import InferenceMode


@dataclass
class InputArgs:
    """Input file arguments."""

    bfile: Optional[Path] = None
    pfile: Optional[Path] = None
    vcf: Optional[Path] = None

    def __post_init__(self) -> None:
        """Validate that at least one input is provided."""
        if self.bfile is None and self.pfile is None and self.vcf is None:
            raise ValueError(
                "At least one input must be provided: --bfile, --pfile, or --vcf"
            )

    @property
    def geno_path(self) -> Path:
        """Get the resolved genotype path (after conversion to pfiles)."""
        if self.pfile is not None:
            return self.pfile
        elif self.bfile is not None:
            return self.bfile
        elif self.vcf is not None:
            # VCF path without .vcf extension
            vcf_str = str(self.vcf)
            if ".vcf" in vcf_str:
                return Path(vcf_str.split(".vcf")[0])
            return self.vcf
        raise ValueError("No input file specified")

    @property
    def input_format(self) -> str:
        """Get the input format type."""
        if self.pfile is not None:
            return "pfile"
        elif self.bfile is not None:
            return "bfile"
        elif self.vcf is not None:
            return "vcf"
        raise ValueError("No input file specified")


@dataclass
class SampleQCArgs:
    """Sample-level QC arguments."""

    # Callrate
    run_callrate: bool = False
    callrate_threshold: float = 0.02

    # Sex check
    run_sex: bool = False
    female_max_f: float = 0.25
    male_min_f: float = 0.75

    # Heterozygosity
    run_het: bool = False
    het_lower: float = -0.15
    het_upper: float = 0.15
    amr_het: bool = False

    # Relatedness
    run_related: bool = False
    related_cutoff: float = 0.0884
    duplicated_cutoff: float = 0.354
    prune_related: bool = False
    prune_duplicated: bool = True

    # Kinship check
    run_kinship_check: bool = False

    def to_callrate_config(self) -> CallrateConfig:
        """Convert to CallrateConfig."""
        return CallrateConfig(mind=self.callrate_threshold)

    def to_sex_config(self) -> SexConfig:
        """Convert to SexConfig."""
        return SexConfig(
            female_max_f=self.female_max_f,
            male_min_f=self.male_min_f,
        )

    def to_het_config(self) -> HetConfig:
        """Convert to HetConfig."""
        return HetConfig(
            f_lower=self.het_lower,
            f_upper=self.het_upper,
        )

    def to_related_config(self) -> RelatedConfig:
        """Convert to RelatedConfig."""
        return RelatedConfig(
            related_cutoff=self.related_cutoff,
            duplicated_cutoff=self.duplicated_cutoff,
            prune_related=self.prune_related,
            prune_duplicated=self.prune_duplicated,
        )


@dataclass
class VariantQCArgs:
    """Variant-level QC arguments."""

    # Genotype missingness
    run_geno: bool = False
    geno_threshold: float = 0.05

    # Case-control
    run_case_control: bool = False
    case_control_threshold: float = 1e-4

    # Haplotype
    run_haplotype: bool = False
    haplotype_threshold: float = 1e-4

    # HWE
    run_hwe: bool = False
    hwe_threshold: float = 1e-4
    filter_controls: bool = False

    # LD pruning
    run_ld: bool = False
    ld_window_size: int = 50
    ld_step_size: int = 5
    ld_r2_threshold: float = 0.5

    def to_geno_config(self) -> GenoConfig:
        """Convert to GenoConfig."""
        return GenoConfig(geno=self.geno_threshold)

    def to_case_control_config(self) -> CaseControlConfig:
        """Convert to CaseControlConfig."""
        return CaseControlConfig(p_threshold=self.case_control_threshold)

    def to_haplotype_config(self) -> HaplotypeConfig:
        """Convert to HaplotypeConfig."""
        return HaplotypeConfig(p_threshold=self.haplotype_threshold)

    def to_hwe_config(self) -> HWEConfig:
        """Convert to HWEConfig."""
        return HWEConfig(
            hwe_threshold=self.hwe_threshold,
            filter_controls=self.filter_controls,
        )

    def to_ld_config(self) -> LDConfig:
        """Convert to LDConfig."""
        return LDConfig(
            window_size=self.ld_window_size,
            step_size=self.ld_step_size,
            r2_threshold=self.ld_r2_threshold,
        )


@dataclass
class AncestryArgs:
    """Ancestry prediction arguments."""

    run_ancestry: bool = False
    ref_panel: Optional[Path] = None
    ref_labels: Optional[Path] = None
    model_path: Optional[Path] = None
    subset_ancestry: Optional[List[str]] = None
    min_samples: int = 0

    # Inference mode
    use_container: bool = False
    use_singularity: bool = False
    use_cloud: bool = False

    @property
    def inference_mode(self) -> InferenceMode:
        """Get the inference mode."""
        if self.use_cloud:
            return InferenceMode.CLOUD
        elif self.use_singularity:
            return InferenceMode.SINGULARITY
        elif self.use_container:
            return InferenceMode.CONTAINER
        return InferenceMode.LOCAL

    def __post_init__(self) -> None:
        """Validate ancestry arguments."""
        if self.use_container and self.model_path is not None:
            raise ValueError(
                "Cannot use both --model and --container. "
                "Container mode uses its own pre-trained model."
            )


@dataclass
class GWASArgs:
    """GWAS and PCA arguments."""

    run_pca: bool = False
    n_pcs: int = 10
    build: str = "hg38"

    run_gwas: bool = False
    covar_path: Optional[Path] = None
    covar_names: Optional[str] = None
    maf_lambdas: bool = False


@dataclass
class OutputArgs:
    """Output configuration arguments."""

    out_path: Path = field(default_factory=lambda: Path("genotools_output"))
    full_output: bool = False
    skip_fails: bool = False
    warn_only: bool = True  # Continue on step failure


@dataclass
class PipelineArgs:
    """Complete validated pipeline arguments.

    This dataclass aggregates all argument groups and provides
    convenience methods for pipeline configuration.
    """

    input: InputArgs
    output: OutputArgs
    sample_qc: SampleQCArgs = field(default_factory=SampleQCArgs)
    variant_qc: VariantQCArgs = field(default_factory=VariantQCArgs)
    ancestry: AncestryArgs = field(default_factory=AncestryArgs)
    gwas: GWASArgs = field(default_factory=GWASArgs)

    @property
    def geno_path(self) -> Path:
        """Get the input genotype path."""
        return self.input.geno_path

    @property
    def out_path(self) -> Path:
        """Get the output path."""
        return self.output.out_path

    @property
    def warn_only(self) -> bool:
        """Get warn_only flag."""
        return self.output.warn_only

    @property
    def full_output(self) -> bool:
        """Get full_output flag."""
        return self.output.full_output

    def get_enabled_sample_steps(self) -> List[str]:
        """Get list of enabled sample QC steps."""
        steps = []
        if self.sample_qc.run_callrate:
            steps.append("callrate")
        if self.sample_qc.run_sex:
            steps.append("sex")
        if self.sample_qc.run_het:
            steps.append("het")
        if self.sample_qc.run_related:
            steps.append("related")
        if self.sample_qc.run_kinship_check:
            steps.append("kinship_check")
        return steps

    def get_enabled_variant_steps(self) -> List[str]:
        """Get list of enabled variant QC steps."""
        steps = []
        if self.variant_qc.run_case_control:
            steps.append("case_control")
        if self.variant_qc.run_haplotype:
            steps.append("haplotype")
        if self.variant_qc.run_hwe:
            steps.append("hwe")
        if self.variant_qc.run_geno:
            steps.append("geno")
        if self.variant_qc.run_ld:
            steps.append("ld")
        return steps

    def get_all_enabled_steps(self) -> List[str]:
        """Get all enabled QC steps in order."""
        steps = self.get_enabled_sample_steps() + self.get_enabled_variant_steps()
        if self.gwas.run_pca or self.gwas.run_gwas:
            steps.append("assoc")
        return steps

    def has_any_qc_steps(self) -> bool:
        """Check if any QC steps are enabled."""
        return len(self.get_enabled_sample_steps()) > 0 or len(
            self.get_enabled_variant_steps()
        ) > 0

    def to_legacy_dict(self) -> dict:
        """Convert to legacy args_dict format for backward compatibility.

        This allows gradual migration by supporting the existing
        execute_pipeline and execute_ancestry_predictions functions.
        """
        return {
            # Input/output
            "bfile": str(self.input.bfile) if self.input.bfile else None,
            "pfile": str(self.input.pfile) if self.input.pfile else None,
            "vcf": str(self.input.vcf) if self.input.vcf else None,
            "geno_path": str(self.geno_path),
            "out": str(self.out_path),
            "full_output": self.full_output,
            "skip_fails": self.output.skip_fails,
            "warn": self.warn_only,
            # Sample QC
            "callrate": self.sample_qc.callrate_threshold
            if self.sample_qc.run_callrate
            else None,
            "sex": [self.sample_qc.female_max_f, self.sample_qc.male_min_f]
            if self.sample_qc.run_sex
            else None,
            "het": [self.sample_qc.het_lower, self.sample_qc.het_upper]
            if self.sample_qc.run_het
            else None,
            "amr_het": self.sample_qc.amr_het,
            "related": self.sample_qc.run_related,
            "related_cutoff": self.sample_qc.related_cutoff,
            "duplicated_cutoff": self.sample_qc.duplicated_cutoff,
            "prune_related": self.sample_qc.prune_related,
            "prune_duplicated": self.sample_qc.prune_duplicated,
            "kinship_check": self.sample_qc.run_kinship_check,
            "all_sample": False,  # Handled by enabling individual steps
            # Variant QC
            "geno": self.variant_qc.geno_threshold
            if self.variant_qc.run_geno
            else None,
            "case_control": self.variant_qc.case_control_threshold
            if self.variant_qc.run_case_control
            else None,
            "haplotype": self.variant_qc.haplotype_threshold
            if self.variant_qc.run_haplotype
            else None,
            "hwe": self.variant_qc.hwe_threshold if self.variant_qc.run_hwe else None,
            "filter_controls": self.variant_qc.filter_controls,
            "ld": [
                self.variant_qc.ld_window_size,
                self.variant_qc.ld_step_size,
                self.variant_qc.ld_r2_threshold,
            ]
            if self.variant_qc.run_ld
            else None,
            "all_variant": False,  # Handled by enabling individual steps
            # Ancestry
            "ancestry": self.ancestry.run_ancestry,
            "ref_panel": str(self.ancestry.ref_panel)
            if self.ancestry.ref_panel
            else None,
            "ref_labels": str(self.ancestry.ref_labels)
            if self.ancestry.ref_labels
            else None,
            "model": str(self.ancestry.model_path)
            if self.ancestry.model_path
            else None,
            "container": self.ancestry.use_container,
            "singularity": self.ancestry.use_singularity,
            "subset_ancestry": self.ancestry.subset_ancestry,
            "min_samples": self.ancestry.min_samples,
            # GWAS
            "pca": self.gwas.n_pcs if self.gwas.run_pca else None,
            "build": self.gwas.build,
            "gwas": self.gwas.run_gwas,
            "covars": str(self.gwas.covar_path) if self.gwas.covar_path else None,
            "covar_names": self.gwas.covar_names,
            "maf_lambdas": self.gwas.maf_lambdas,
        }


def create_parser() -> argparse.ArgumentParser:
    """Create argument parser with proper types.

    Returns:
        Configured ArgumentParser instance.
    """
    parser = argparse.ArgumentParser(
        prog="genotools",
        description="Genotype QC and ancestry prediction pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all sample QC
  genotools --pfile data/mydata --out results/output --all-sample

  # Run ancestry prediction only
  genotools --pfile data/mydata --out results/output --ancestry \\
            --ref-panel refs/panel --ref-labels refs/labels.txt

  # Run full pipeline with ancestry split
  genotools --pfile data/mydata --out results/output --ancestry \\
            --ref-panel refs/panel --ref-labels refs/labels.txt \\
            --all-sample --all-variant
        """,
    )

    # Input file group
    input_group = parser.add_argument_group("Input files")
    input_group.add_argument(
        "--bfile",
        type=Path,
        default=None,
        help="PLINK 1.9 format genotype file path (without .bed/.bim/.fam extension)",
    )
    input_group.add_argument(
        "--pfile",
        type=Path,
        default=None,
        help="PLINK 2 format genotype file path (without .pgen/.pvar/.psam extension)",
    )
    input_group.add_argument(
        "--vcf",
        type=Path,
        default=None,
        help="VCF format genotype file path",
    )

    # Output group
    output_group = parser.add_argument_group("Output configuration")
    output_group.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output prefix (including path)",
    )
    output_group.add_argument(
        "--full-output",
        action="store_true",
        help="Output all intermediate files",
    )
    output_group.add_argument(
        "--skip-fails",
        action="store_true",
        help="Skip up-front validation checks",
    )
    output_group.add_argument(
        "--no-warn",
        action="store_true",
        help="Stop pipeline on first error (default: continue with warnings)",
    )

    # Sample QC group
    sample_group = parser.add_argument_group("Sample-level QC")
    sample_group.add_argument(
        "--callrate",
        type=float,
        nargs="?",
        const=0.02,
        default=None,
        metavar="THRESHOLD",
        help="Sample callrate threshold (default: 0.02 if flag used)",
    )
    sample_group.add_argument(
        "--sex",
        type=float,
        nargs="*",
        default=None,
        metavar="VALUE",
        help="Sex check thresholds [female_max_f male_min_f] (default: 0.25 0.75)",
    )
    sample_group.add_argument(
        "--het",
        type=float,
        nargs="*",
        default=None,
        metavar="VALUE",
        help="Heterozygosity thresholds [lower upper] (default: -0.15 0.15)",
    )
    sample_group.add_argument(
        "--amr-het",
        action="store_true",
        help="Use auto-detect heterozygosity for AMR samples",
    )
    sample_group.add_argument(
        "--related",
        action="store_true",
        help="Run relatedness analysis",
    )
    sample_group.add_argument(
        "--related-cutoff",
        type=float,
        default=0.0884,
        metavar="CUTOFF",
        help="Relatedness cutoff (default: 0.0884)",
    )
    sample_group.add_argument(
        "--duplicated-cutoff",
        type=float,
        default=0.354,
        metavar="CUTOFF",
        help="Duplicate cutoff (default: 0.354)",
    )
    sample_group.add_argument(
        "--prune-related",
        action="store_true",
        help="Prune related samples",
    )
    sample_group.add_argument(
        "--no-prune-duplicated",
        action="store_true",
        help="Do not prune duplicated samples",
    )
    sample_group.add_argument(
        "--kinship-check",
        action="store_true",
        help="Run kinship confirmation (Linux only)",
    )
    sample_group.add_argument(
        "--all-sample",
        action="store_true",
        help="Run all sample-level QC with default thresholds",
    )

    # Variant QC group
    variant_group = parser.add_argument_group("Variant-level QC")
    variant_group.add_argument(
        "--geno",
        type=float,
        nargs="?",
        const=0.05,
        default=None,
        metavar="THRESHOLD",
        help="Variant missingness threshold (default: 0.05 if flag used)",
    )
    variant_group.add_argument(
        "--case-control",
        type=float,
        nargs="?",
        const=1e-4,
        default=None,
        metavar="P_THRESHOLD",
        help="Case-control missingness p-value threshold (default: 1e-4)",
    )
    variant_group.add_argument(
        "--haplotype",
        type=float,
        nargs="?",
        const=1e-4,
        default=None,
        metavar="P_THRESHOLD",
        help="Haplotype test p-value threshold (default: 1e-4)",
    )
    variant_group.add_argument(
        "--hwe",
        type=float,
        nargs="?",
        const=1e-4,
        default=None,
        metavar="P_THRESHOLD",
        help="HWE p-value threshold (default: 1e-4)",
    )
    variant_group.add_argument(
        "--filter-controls",
        action="store_true",
        help="Apply HWE filter to controls only",
    )
    variant_group.add_argument(
        "--ld",
        type=float,
        nargs="*",
        default=None,
        metavar="VALUE",
        help="LD pruning [window_size step_size r2_threshold] (default: 50 5 0.5)",
    )
    variant_group.add_argument(
        "--all-variant",
        action="store_true",
        help="Run all variant-level QC with default thresholds (except LD)",
    )

    # Ancestry group
    ancestry_group = parser.add_argument_group("Ancestry prediction")
    ancestry_group.add_argument(
        "--ancestry",
        action="store_true",
        help="Run ancestry prediction",
    )
    ancestry_group.add_argument(
        "--ref-panel",
        type=Path,
        default=None,
        metavar="PATH",
        help="Reference panel genotype file path",
    )
    ancestry_group.add_argument(
        "--ref-labels",
        type=Path,
        default=None,
        metavar="PATH",
        help="Reference panel ancestry labels file (FID IID label, no header)",
    )
    ancestry_group.add_argument(
        "--model",
        type=Path,
        default=None,
        metavar="PATH",
        help="Pre-trained ancestry model pickle file",
    )
    ancestry_group.add_argument(
        "--container",
        action="store_true",
        help="Run ancestry prediction in Docker container",
    )
    ancestry_group.add_argument(
        "--singularity",
        action="store_true",
        help="Run ancestry prediction in Singularity container",
    )
    ancestry_group.add_argument(
        "--cloud",
        action="store_true",
        help="Run ancestry prediction on Google Cloud AI Platform",
    )
    ancestry_group.add_argument(
        "--subset-ancestry",
        type=str,
        nargs="*",
        default=None,
        metavar="LABEL",
        help="Subset to specific ancestry labels for downstream analysis",
    )
    ancestry_group.add_argument(
        "--min-samples",
        type=int,
        default=0,
        metavar="N",
        help="Minimum samples per ancestry for downstream analysis (default: 0)",
    )

    # GWAS group
    gwas_group = parser.add_argument_group("GWAS and PCA")
    gwas_group.add_argument(
        "--pca",
        type=int,
        nargs="?",
        const=10,
        default=None,
        metavar="N_PCS",
        help="Run PCA with N principal components (default: 10)",
    )
    gwas_group.add_argument(
        "--build",
        type=str,
        default="hg38",
        choices=["hg19", "hg38"],
        help="Genome build for PCA (default: hg38)",
    )
    gwas_group.add_argument(
        "--gwas",
        action="store_true",
        help="Run GWAS",
    )
    gwas_group.add_argument(
        "--covars",
        type=Path,
        default=None,
        metavar="PATH",
        help="External covariates file path",
    )
    gwas_group.add_argument(
        "--covar-names",
        type=str,
        default=None,
        metavar="NAMES",
        help="Covariate names to use from external file",
    )
    gwas_group.add_argument(
        "--maf-lambdas",
        action="store_true",
        help="MAF prune before lambda calculations",
    )

    return parser


def _parse_sex_args(
    sex_values: Optional[List[float]],
) -> Tuple[bool, float, float]:
    """Parse sex check arguments.

    Args:
        sex_values: List of threshold values or None.

    Returns:
        Tuple of (run_sex, female_max_f, male_min_f).
    """
    if sex_values is None:
        return False, 0.25, 0.75

    if len(sex_values) == 0:
        # Flag used without values - use defaults
        return True, 0.25, 0.75
    elif len(sex_values) == 2:
        return True, sex_values[0], sex_values[1]
    else:
        raise ValueError(
            f"--sex requires 0 or 2 values, got {len(sex_values)}: {sex_values}"
        )


def _parse_het_args(
    het_values: Optional[List[float]],
) -> Tuple[bool, float, float]:
    """Parse heterozygosity arguments.

    Args:
        het_values: List of threshold values or None.

    Returns:
        Tuple of (run_het, het_lower, het_upper).
    """
    if het_values is None:
        return False, -0.15, 0.15

    if len(het_values) == 0:
        # Flag used without values - use defaults
        return True, -0.15, 0.15
    elif len(het_values) == 2:
        return True, het_values[0], het_values[1]
    else:
        raise ValueError(
            f"--het requires 0 or 2 values, got {len(het_values)}: {het_values}"
        )


def _parse_ld_args(
    ld_values: Optional[List[float]],
) -> Tuple[bool, int, int, float]:
    """Parse LD pruning arguments.

    Args:
        ld_values: List of LD parameters or None.

    Returns:
        Tuple of (run_ld, window_size, step_size, r2_threshold).
    """
    if ld_values is None:
        return False, 50, 5, 0.5

    if len(ld_values) == 0:
        # Flag used without values - use defaults
        return True, 50, 5, 0.5
    elif len(ld_values) == 3:
        return True, int(ld_values[0]), int(ld_values[1]), ld_values[2]
    else:
        raise ValueError(
            f"--ld requires 0 or 3 values, got {len(ld_values)}: {ld_values}"
        )


def parse_args(args: Optional[Sequence[str]] = None) -> PipelineArgs:
    """Parse and validate command-line arguments.

    Args:
        args: Command-line arguments (defaults to sys.argv).

    Returns:
        Validated PipelineArgs instance.

    Raises:
        SystemExit: If argument parsing fails.
        ValueError: If argument validation fails.
    """
    parser = create_parser()
    ns = parser.parse_args(args)

    # Validate input files
    if ns.bfile is None and ns.pfile is None and ns.vcf is None:
        parser.error("At least one input required: --bfile, --pfile, or --vcf")

    # Parse complex arguments
    run_sex, female_max_f, male_min_f = _parse_sex_args(ns.sex)
    run_het, het_lower, het_upper = _parse_het_args(ns.het)
    run_ld, ld_window, ld_step, ld_r2 = _parse_ld_args(ns.ld)

    # Handle --all-sample
    run_callrate = ns.callrate is not None
    callrate_threshold = ns.callrate if ns.callrate is not None else 0.02
    run_related = ns.related

    if ns.all_sample:
        run_callrate = True
        callrate_threshold = 0.05  # Different default for all_sample
        run_sex = True
        run_het = True
        run_related = True

    # Handle --all-variant
    run_geno = ns.geno is not None
    geno_threshold = ns.geno if ns.geno is not None else 0.05
    run_case_control = ns.case_control is not None
    case_control_threshold = ns.case_control if ns.case_control is not None else 1e-4
    run_haplotype = ns.haplotype is not None
    haplotype_threshold = ns.haplotype if ns.haplotype is not None else 1e-4
    run_hwe = ns.hwe is not None
    hwe_threshold = ns.hwe if ns.hwe is not None else 1e-4

    if ns.all_variant:
        run_geno = True
        geno_threshold = 0.05
        run_case_control = True
        case_control_threshold = 1e-4
        run_haplotype = True
        haplotype_threshold = 1e-4
        run_hwe = True
        hwe_threshold = 1e-4
        # Note: LD is NOT included in all_variant (matches legacy behavior)

    # Build argument groups
    input_args = InputArgs(
        bfile=ns.bfile,
        pfile=ns.pfile,
        vcf=ns.vcf,
    )

    output_args = OutputArgs(
        out_path=ns.out,
        full_output=ns.full_output,
        skip_fails=ns.skip_fails,
        warn_only=not ns.no_warn,
    )

    sample_qc_args = SampleQCArgs(
        run_callrate=run_callrate,
        callrate_threshold=callrate_threshold,
        run_sex=run_sex,
        female_max_f=female_max_f,
        male_min_f=male_min_f,
        run_het=run_het,
        het_lower=het_lower,
        het_upper=het_upper,
        amr_het=ns.amr_het,
        run_related=run_related,
        related_cutoff=ns.related_cutoff,
        duplicated_cutoff=ns.duplicated_cutoff,
        prune_related=ns.prune_related,
        prune_duplicated=not ns.no_prune_duplicated,
        run_kinship_check=ns.kinship_check,
    )

    variant_qc_args = VariantQCArgs(
        run_geno=run_geno,
        geno_threshold=geno_threshold,
        run_case_control=run_case_control,
        case_control_threshold=case_control_threshold,
        run_haplotype=run_haplotype,
        haplotype_threshold=haplotype_threshold,
        run_hwe=run_hwe,
        hwe_threshold=hwe_threshold,
        filter_controls=ns.filter_controls,
        run_ld=run_ld,
        ld_window_size=ld_window,
        ld_step_size=ld_step,
        ld_r2_threshold=ld_r2,
    )

    ancestry_args = AncestryArgs(
        run_ancestry=ns.ancestry,
        ref_panel=ns.ref_panel,
        ref_labels=ns.ref_labels,
        model_path=ns.model,
        subset_ancestry=ns.subset_ancestry,
        min_samples=ns.min_samples,
        use_container=ns.container,
        use_singularity=ns.singularity,
        use_cloud=ns.cloud,
    )

    gwas_args = GWASArgs(
        run_pca=ns.pca is not None,
        n_pcs=ns.pca if ns.pca is not None else 10,
        build=ns.build,
        run_gwas=ns.gwas,
        covar_path=ns.covars,
        covar_names=ns.covar_names,
        maf_lambdas=ns.maf_lambdas,
    )

    return PipelineArgs(
        input=input_args,
        output=output_args,
        sample_qc=sample_qc_args,
        variant_qc=variant_qc_args,
        ancestry=ancestry_args,
        gwas=gwas_args,
    )
