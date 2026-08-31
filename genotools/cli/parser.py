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
from typing import Dict, List, Optional, Tuple, Sequence

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


@dataclass(frozen=True)
class HetSpec:
    """One resolved heterozygosity specification.

    Both ``--het`` and ``--het-ancestry LABEL`` accept the same four-row
    grammar, so both resolve to one of these:

    ==================  ====================================================
    Tokens              Meaning
    ==================  ====================================================
    *(none)*            fixed ``-0.15 0.15`` (the defaults)
    ``LOWER UPPER``     fixed bounds on PLINK's F
    ``sd``              ``mean +/- 3 * sd`` of this dataset's F
    ``sd N``            ``mean +/- N * sd``
    ==================  ====================================================

    Sharing one spec type is what keeps the two flags from drifting into
    different rules: :func:`_parse_het_spec` is the only thing that builds one.
    """

    auto: bool = False
    sd: float = 3.0
    lower: float = -0.15
    upper: float = 0.15

    def to_het_config(self) -> HetConfig:
        """Convert to the step-level HetConfig.

        ``lower``/``upper`` ride along unused when ``auto`` is set. They are
        deliberately *not* forced to the legacy ``[-1, -1]`` sentinel, which
        would select HetConfig's rate-based path instead of the F-based one.
        """
        return HetConfig(
            f_lower=self.lower,
            f_upper=self.upper,
            auto_detect=self.auto,
            auto_sd=self.sd,
        )


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

    # Heterozygosity. het_lower/het_upper (or het_auto/het_auto_sd) are the
    # base, applying to every dataset; het_by_ancestry overrides it for named
    # ancestry groups. Resolve the pair through het_config_for -- never read
    # them apart, or the two runner paths can disagree again.
    run_het: bool = False
    het_lower: float = -0.15
    het_upper: float = 0.15
    het_auto: bool = False
    het_auto_sd: float = 3.0
    het_by_ancestry: Dict[str, HetSpec] = field(default_factory=dict)

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

    def het_base_spec(self) -> HetSpec:
        """The base spec, applying wherever no per-ancestry override does."""
        return HetSpec(
            auto=self.het_auto,
            sd=self.het_auto_sd,
            lower=self.het_lower,
            upper=self.het_upper,
        )

    def het_config_for(self, label: Optional[str] = None) -> HetConfig:
        """Resolve the heterozygosity config for one dataset.

        The single point where base and per-ancestry override are combined.
        Both runner paths call it -- the flat run with ``label=None``, the
        ancestry loop with each group's label -- so neither can decide het
        config on its own. ``--amr-het`` was silently inert outside
        ``--ancestry`` precisely because the two paths *did* decide separately
        and only one of them knew about the flag.

        Args:
            label: Ancestry group being configured, or None for a flat run
                (which has no labels, so only the base can apply).

        Returns:
            HetConfig for this dataset: the override for ``label`` if one was
            given, otherwise the base.
        """
        spec = self.het_by_ancestry.get(label) if label is not None else None
        if spec is None:
            spec = self.het_base_spec()
        return spec.to_het_config()

    def to_het_config(self) -> HetConfig:
        """Convert to HetConfig (the base, with no ancestry override)."""
        return self.het_config_for(None)

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


#: Remote-execution flags that 2.0 does not implement, mapped to
#: ``(AncestryArgs attribute, explanatory detail)``.
#:
#: 1.x ran ancestry prediction inside a Docker/Singularity image whose
#: ``run.py`` unpickles a 1.x ``umap_linearsvc`` model. 2.0's ``AncestryModel``
#: cannot load that format, so honouring ``--container`` would need a rebuilt,
#: republished image carrying a 2.0-format model. Until that exists these flags
#: fail loudly rather than silently running local prediction.
_UNSUPPORTED_INFERENCE_FLAGS: Dict[str, Tuple[str, str]] = {
    "--container": (
        "use_container",
        "1.x ran prediction in a Docker image built around a 1.x model, which "
        "2.0 cannot load. Drop the flag to predict locally, or pin "
        "'genotools<2.0' to keep the 1.x container. See MIGRATION_2.0.md.",
    ),
    "--singularity": (
        "use_singularity",
        "1.x ran prediction in a Singularity image built around a 1.x model, "
        "which 2.0 cannot load. Drop the flag to predict locally, or pin "
        "'genotools<2.0' to keep the 1.x container. See MIGRATION_2.0.md.",
    ),
    "--cloud": (
        "use_cloud",
        "Cloud prediction has never been implemented -- not in 2.0 and not in "
        "1.x. Drop the flag to predict locally. See MIGRATION_2.0.md.",
    ),
}


@dataclass
class AncestryArgs:
    """Ancestry prediction arguments."""

    run_ancestry: bool = False
    ref_panel: Optional[Path] = None
    ref_labels: Optional[Path] = None
    model_path: Optional[Path] = None
    subset_ancestry: Optional[List[str]] = None
    min_samples: int = 0

    # Remote-execution flags. Accepted by the parser so that a 1.x command line
    # gets a targeted error instead of argparse's bare "unrecognized arguments",
    # but rejected in __post_init__ -- see _UNSUPPORTED_INFERENCE_FLAGS.
    use_container: bool = False
    use_singularity: bool = False
    use_cloud: bool = False

    def __post_init__(self) -> None:
        """Validate ancestry arguments."""
        for flag, (attr, detail) in _UNSUPPORTED_INFERENCE_FLAGS.items():
            if getattr(self, attr):
                raise ValueError(f"{flag} is not supported in GenoTools 2.0. {detail}")


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
    quiet: bool = False  # Suppress the console stream (log files still written)
    debug: bool = False  # Widen both streams to DEBUG level


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

    def __post_init__(self) -> None:
        """Validate arguments that span more than one group."""
        # A flat run has no ancestry labels, so a per-label override cannot be
        # honoured. Erroring is the point of this feature: --amr-het's defining
        # bug was accepting exactly this and quietly doing nothing.
        if self.sample_qc.het_by_ancestry and not self.ancestry.run_ancestry:
            labels = ", ".join(sorted(self.sample_qc.het_by_ancestry))
            raise ValueError(
                f"--het-ancestry ({labels}) requires --ancestry: a run without "
                f"ancestry prediction has no ancestry labels to match. To set "
                f"bounds for the whole input, use --het instead."
            )

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
        if self.variant_qc.run_geno:
            steps.append("geno")
        if self.variant_qc.run_case_control:
            steps.append("case_control")
        if self.variant_qc.run_haplotype:
            steps.append("haplotype")
        if self.variant_qc.run_hwe:
            steps.append("hwe")
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
            "het_auto": self.sample_qc.het_auto,
            "het_auto_sd": self.sample_qc.het_auto_sd,
            "het_by_ancestry": dict(self.sample_qc.het_by_ancestry),
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
        "--full_output",  # deprecated underscore spelling
        action="store_true",
        help="Output all intermediate files",
    )
    output_group.add_argument(
        "--skip-fails",
        "--skip_fails",  # deprecated underscore spelling
        action="store_true",
        help="Skip up-front validation checks",
    )
    output_group.add_argument(
        "--no-warn",
        action="store_true",
        help="Stop pipeline on first error (default: continue with warnings)",
    )
    # Deprecated presence-only alias: 1.x defaulted warn on, which is now the
    # default, so passing it changes nothing. Use --no-warn to disable.
    output_group.add_argument(
        "--warn",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    output_group.add_argument(
        "--quiet",
        action="store_true",
        help="Suppress console progress output (log files are still written)",
    )
    output_group.add_argument(
        "--debug",
        action="store_true",
        help="Log at DEBUG level (console and consolidated log)",
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
    # type=str, not float: the spec grammar mixes numbers with the "sd"
    # keyword, so coercion moves into _parse_het_spec. argparse still parses
    # leading negatives here -- its negative-number heuristic keys off whether
    # the parser has number-like option strings, not off the argument type,
    # and GenoTools has none (nor any positionals for greedy nargs to eat).
    sample_group.add_argument(
        "--het",
        type=str,
        nargs="*",
        default=None,
        metavar="SPEC",
        help=(
            "Heterozygosity bounds for every dataset: [LOWER UPPER] fixed "
            "bounds on F, or 'sd [N]' for bounds at N standard deviations of "
            "this dataset's own F (default: -0.15 0.15, or 3 sd for bare 'sd')"
        ),
    )
    sample_group.add_argument(
        # No underscore alias: underscores exist only as deprecated 1.x
        # spellings, and this flag is new. Same as every other 2.0 addition.
        "--het-ancestry",
        type=str,
        nargs="+",
        action="append",
        default=None,
        metavar="SPEC",
        help=(
            "Override --het for one ancestry group: LABEL followed by the "
            "same spec ('sd [N]' or LOWER UPPER). Repeatable, one group per "
            "use. Requires --ancestry"
        ),
    )
    sample_group.add_argument(
        "--related",
        action="store_true",
        help="Run relatedness analysis",
    )
    sample_group.add_argument(
        "--related-cutoff",
        "--related_cutoff",  # deprecated underscore spelling
        type=float,
        default=0.0884,
        metavar="CUTOFF",
        help="Relatedness cutoff (default: 0.0884)",
    )
    sample_group.add_argument(
        "--duplicated-cutoff",
        "--duplicated_cutoff",  # deprecated underscore spelling
        type=float,
        default=0.354,
        metavar="CUTOFF",
        help="Duplicate cutoff (default: 0.354)",
    )
    sample_group.add_argument(
        "--prune-related",
        "--prune_related",  # deprecated underscore spelling
        action="store_true",
        help="Prune related samples",
    )
    # Deprecated presence-only alias: 1.x defaulted duplicate pruning on, which
    # is now the default. Use --no-prune-duplicated to disable.
    sample_group.add_argument(
        "--prune_duplicated",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    sample_group.add_argument(
        "--no-prune-duplicated",
        action="store_true",
        help="Do not prune duplicated samples",
    )
    sample_group.add_argument(
        "--kinship-check",
        "--kinship_check",  # deprecated underscore spelling
        action="store_true",
        help="Run kinship confirmation (Linux only)",
    )
    sample_group.add_argument(
        "--all-sample",
        "--all_sample",  # deprecated underscore spelling
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
        "--case_control",  # deprecated underscore spelling
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
        "--filter_controls",  # deprecated underscore spelling
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
        "--all_variant",  # deprecated underscore spelling
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
        "--ref_panel",  # deprecated underscore spelling
        type=Path,
        default=None,
        metavar="PATH",
        help="Reference panel genotype file path",
    )
    ancestry_group.add_argument(
        "--ref-labels",
        "--ref_labels",  # deprecated underscore spelling
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
        help="Pre-trained ancestry model: a directory, or a single .pkl "
             "written by GenoTools 2.0 (1.x models cannot be loaded)",
    )
    ancestry_group.add_argument(
        "--container",
        action="store_true",
        help="Not supported in 2.0 (1.x-only Docker prediction); errors if passed",
    )
    ancestry_group.add_argument(
        "--singularity",
        action="store_true",
        help="Not supported in 2.0 (1.x-only Singularity prediction); errors if passed",
    )
    ancestry_group.add_argument(
        "--cloud",
        action="store_true",
        help="Not supported in 2.0 (never implemented); errors if passed",
    )
    ancestry_group.add_argument(
        "--subset-ancestry",
        "--subset_ancestry",  # deprecated underscore spelling
        type=str,
        nargs="*",
        default=None,
        metavar="LABEL",
        help="Subset to specific ancestry labels for downstream analysis",
    )
    ancestry_group.add_argument(
        "--min-samples",
        "--min_samples",  # deprecated underscore spelling
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
        "--covar_names",  # deprecated underscore spelling
        type=str,
        default=None,
        metavar="NAMES",
        help="Covariate names to use from external file",
    )
    gwas_group.add_argument(
        "--maf-lambdas",
        "--maf_lambdas",  # deprecated underscore spelling
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


#: Spelling of the data-derived mode on the command line. The config layer
#: calls it ``auto_detect``, which predates the multiplier being user-visible;
#: "auto 3" would read as "automatic... three?". Only the location and scale
#: are derived from the data - the multiplier is a human choice - so the CLI
#: says ``sd``, and ``--het sd 2`` says exactly what it does.
_HET_SD_TOKEN = "sd"


def _het_float(token: str, flag: str) -> float:
    """Coerce one spec token to a float with an actionable error.

    ``--het`` is ``type=str`` so the grammar can carry the ``sd`` keyword,
    which means argparse no longer supplies the type error. This does.
    """
    try:
        return float(token)
    except ValueError:
        raise ValueError(
            f"{flag} expects numbers or '{_HET_SD_TOKEN}', got {token!r}"
        ) from None


def _parse_het_spec(tokens: Sequence[str], flag: str) -> HetSpec:
    """Resolve one heterozygosity spec - the single grammar for both flags.

    ``--het`` and ``--het-ancestry LABEL`` share this so they cannot drift into
    different rules or different error wording. The forms never collide:
    ``sd`` is a keyword, so a two-token spec is fixed bounds when it starts
    with a number and a multiplier when it starts with ``sd``.

    Args:
        tokens: Spec tokens, with any ancestry label already stripped.
        flag: Flag name, for error messages.

    Returns:
        The resolved HetSpec.

    Raises:
        ValueError: On unparseable tokens, a bad arity, a non-positive
            multiplier, or bounds that are not ordered.
    """
    if len(tokens) == 0:
        return HetSpec()

    if tokens[0].lower() == _HET_SD_TOKEN:
        if len(tokens) == 1:
            return HetSpec(auto=True)
        if len(tokens) == 2:
            sd = _het_float(tokens[1], flag)
            if sd <= 0:
                raise ValueError(
                    f"{flag} '{_HET_SD_TOKEN}' multiplier must be greater "
                    f"than 0, got {sd}"
                )
            return HetSpec(auto=True, sd=sd)
        raise ValueError(
            f"{flag} '{_HET_SD_TOKEN}' takes at most one multiplier, got "
            f"{len(tokens) - 1}: {list(tokens[1:])}"
        )

    if len(tokens) == 2:
        lower = _het_float(tokens[0], flag)
        upper = _het_float(tokens[1], flag)
        # The 1.x sentinel for "derive the bounds", still reachable on old
        # command lines. It has a real keyword now, so accept but redirect:
        # a silent second spelling is worse than a deprecation.
        if lower == -1.0 and upper == -1.0:
            from ..core.logging import get_logger

            get_logger(__name__).warning(
                f"{flag} -1 -1 is deprecated; use "
                f"'{flag} {_HET_SD_TOKEN}' instead."
            )
            return HetSpec(auto=True)
        if lower >= upper:
            raise ValueError(
                f"{flag} lower bound ({lower}) must be less than the upper "
                f"bound ({upper})"
            )
        return HetSpec(lower=lower, upper=upper)

    # Report an unparseable token as a type error before falling back to the
    # arity message: "--het abc" wanted 'abc' explained, not counted.
    for token in tokens:
        _het_float(token, flag)
    raise ValueError(
        f"{flag} requires 0 or 2 values, or '{_HET_SD_TOKEN}' [N], got "
        f"{len(tokens)}: {list(tokens)}"
    )


def _parse_het_args(
    het_values: Optional[List[str]],
) -> Tuple[bool, HetSpec]:
    """Parse the base --het argument.

    Args:
        het_values: Raw spec tokens, or None when the flag was not passed.

    Returns:
        Tuple of (run_het, spec).
    """
    if het_values is None:
        return False, HetSpec()
    return True, _parse_het_spec(het_values, "--het")


def _parse_het_ancestry_args(
    occurrences: Optional[List[List[str]]],
) -> Dict[str, HetSpec]:
    """Parse repeated --het-ancestry occurrences into per-label specs.

    Args:
        occurrences: One token list per use of the flag, each beginning with
            an ancestry label, or None when the flag was not passed.

    Returns:
        Mapping of ancestry label to its spec. Empty when the flag is unused.

    Raises:
        ValueError: On a missing spec, a label that looks like a spec token,
            or the same label overridden twice.
    """
    specs: Dict[str, HetSpec] = {}
    for tokens in occurrences or []:
        label, rest = tokens[0], tokens[1:]
        # Guard the easy transposition: --het-ancestry sd 2 names no group.
        if label.lower() == _HET_SD_TOKEN or _is_number(label):
            raise ValueError(
                f"--het-ancestry expects an ancestry label first, got "
                f"{label!r}. Usage: --het-ancestry LABEL "
                f"['{_HET_SD_TOKEN}' [N] | LOWER UPPER]"
            )
        # A bare label would silently mean "fixed defaults for this group",
        # which is the class of quiet no-op this flag exists to remove.
        if not rest:
            raise ValueError(
                f"--het-ancestry {label} needs a spec: "
                f"'{_HET_SD_TOKEN}' [N], or LOWER UPPER"
            )
        if label in specs:
            raise ValueError(
                f"--het-ancestry {label} given more than once; pass one "
                f"spec per ancestry group"
            )
        specs[label] = _parse_het_spec(rest, f"--het-ancestry {label}")
    return specs


def _is_number(token: str) -> bool:
    """Whether a token parses as a float."""
    try:
        float(token)
    except ValueError:
        return False
    return True


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


# Flags whose spelling changed in 2.0. The old spelling is still accepted so
# existing scripts keep working; using one emits a single deprecation warning.
_DEPRECATED_SPELLINGS = {
    "--all_sample": "--all-sample",
    "--all_variant": "--all-variant",
    "--case_control": "--case-control",
    "--covar_names": "--covar-names",
    "--duplicated_cutoff": "--duplicated-cutoff",
    "--filter_controls": "--filter-controls",
    "--full_output": "--full-output",
    "--kinship_check": "--kinship-check",
    "--maf_lambdas": "--maf-lambdas",
    "--min_samples": "--min-samples",
    "--prune_related": "--prune-related",
    "--ref_labels": "--ref-labels",
    "--ref_panel": "--ref-panel",
    "--related_cutoff": "--related-cutoff",
    "--skip_fails": "--skip-fails",
    "--subset_ancestry": "--subset-ancestry",
}

#: Flags removed outright, mapped to what to use instead. Declared here rather
#: than in the parser so a 1.x command line gets a targeted error naming the
#: replacement instead of argparse's bare "unrecognized arguments".
#:
#: ``--amr-het`` switched the het step to data-derived bounds for one
#: hardcoded label. It was read in exactly one place, inside the --ancestry
#: branch, so it did nothing at all in a flat run - which is how the
#: per-ancestry production workflow runs QC. ``--het-ancestry AMR sd`` is the
#: same rule, on any label, in either run shape.
_REMOVED_FLAGS: Dict[str, str] = {
    "--amr-het": (
        "Use '--het-ancestry AMR sd' instead, which also works in a flat run "
        "(--amr-het silently did nothing without --ancestry). See "
        "MIGRATION_2.0.md."
    ),
    "--amr_het": (
        "Use '--het-ancestry AMR sd' instead, which also works in a flat run "
        "(--amr_het silently did nothing without --ancestry). See "
        "MIGRATION_2.0.md."
    ),
}


def _reject_removed_flags(args: Optional[Sequence[str]] = None) -> None:
    """Give an actionable error for a flag that 2.1 removed.

    Runs before argparse so the message survives - and before
    ``_reject_boolean_values``, so ``--amr-het False`` reports the removal
    rather than the (now irrelevant) presence-flag rule.
    """
    import sys

    argv = list(args) if args is not None else sys.argv[1:]
    used = {tok.split("=", 1)[0] for tok in argv if tok.startswith("--")}
    for flag in sorted(used & set(_REMOVED_FLAGS)):
        raise ValueError(
            f"{flag} was removed in GenoTools 2.0.1. {_REMOVED_FLAGS[flag]}"
        )


# Flags replaced by a --no- inverse. 1.x defaulted both to on, which is now the
# default, so passing them changes nothing.
_DEPRECATED_BOOLEANS = {
    "--warn": (
        "warn-and-continue is now the default; pass --no-warn to stop on the "
        "first error"
    ),
    "--prune_duplicated": (
        "duplicate pruning is now the default; pass --no-prune-duplicated to "
        "disable it"
    ),
}

# 1.x declared every boolean flag as `type=str, nargs='?', const='True'`,
# copying the shape of the threshold flags (--callrate and friends) where an
# optional value is genuinely wanted. A presence flag has no value to override,
# so that pattern only made `--all_sample False` parse and then quietly mean
# "off". 2.0 declares them action="store_true", which rejects a value - but
# argparse reports it as a bare "unrecognized arguments: False", which does not
# tell anyone migrating what to do instead.
_BOOLEAN_VALUE_TOKENS = frozenset({"True", "False", "true", "false"})


def _reject_boolean_values(
    parser: argparse.ArgumentParser, args: Optional[Sequence[str]] = None
) -> None:
    """Give an actionable error for the 1.x ``--flag True/False`` form.

    Derives the presence-only flags from the parser itself (an action that
    consumes no argument has ``nargs == 0``) so any flag added later is covered
    without touching this list.
    """
    import sys

    presence_only = {
        opt
        for action in parser._actions
        if action.nargs == 0
        for opt in action.option_strings
    }

    argv = list(args) if args is not None else sys.argv[1:]
    for flag, value in zip(argv, argv[1:]):
        if flag not in presence_only or value not in _BOOLEAN_VALUE_TOKENS:
            continue
        hint = _DEPRECATED_BOOLEANS.get(flag)
        if hint is None:
            hint = (
                "pass the flag on its own to enable it, or omit it to disable "
                "it"
            )
        parser.error(
            f"{flag} takes no value (got {value!r}). In 1.x every boolean flag "
            f"accepted True/False; they are now presence-only. {hint}."
        )


def _warn_deprecated_flags(args: Optional[Sequence[str]] = None) -> None:
    """Log one warning naming every deprecated flag spelling that was used."""
    import sys

    argv = list(args) if args is not None else sys.argv[1:]
    # Only the flag token matters; strip any "=value" form.
    used = {tok.split("=", 1)[0] for tok in argv if tok.startswith("--")}

    renamed = sorted(used & set(_DEPRECATED_SPELLINGS))
    booleans = sorted(used & set(_DEPRECATED_BOOLEANS))
    if not renamed and not booleans:
        return

    from ..core.logging import get_logger

    logger = get_logger(__name__)
    for flag in renamed:
        logger.warning(
            f"{flag} is deprecated and will be removed in a future release; "
            f"use {_DEPRECATED_SPELLINGS[flag]} instead."
        )
    for flag in booleans:
        logger.warning(
            f"{flag} is deprecated: {_DEPRECATED_BOOLEANS[flag]}."
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
    _reject_removed_flags(args)
    _reject_boolean_values(parser, args)
    ns = parser.parse_args(args)
    _warn_deprecated_flags(args)

    # Validate input files
    if ns.bfile is None and ns.pfile is None and ns.vcf is None:
        parser.error("At least one input required: --bfile, --pfile, or --vcf")

    # Parse complex arguments
    run_sex, female_max_f, male_min_f = _parse_sex_args(ns.sex)
    run_het, het_spec = _parse_het_args(ns.het)
    het_by_ancestry = _parse_het_ancestry_args(ns.het_ancestry)
    run_ld, ld_window, ld_step, ld_r2 = _parse_ld_args(ns.ld)

    # An override implies the system it overrides is on. (Unlike the round-10
    # treatment of --container, which fails loudly: those flags could not work,
    # this one can.)
    if het_by_ancestry:
        run_het = True

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
        quiet=ns.quiet,
        debug=ns.debug,
    )

    sample_qc_args = SampleQCArgs(
        run_callrate=run_callrate,
        callrate_threshold=callrate_threshold,
        run_sex=run_sex,
        female_max_f=female_max_f,
        male_min_f=male_min_f,
        run_het=run_het,
        het_lower=het_spec.lower,
        het_upper=het_spec.upper,
        het_auto=het_spec.auto,
        het_auto_sd=het_spec.sd,
        het_by_ancestry=het_by_ancestry,
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

    # Ancestry prediction needs the reference panel in *both* modes. Training
    # builds the reference PCA from it; inference (--model) still re-derives
    # that PCA, subsetting the panel to the model's common SNPs
    # (ancestry/preprocessing.py) and reading the labels TSV. Without this
    # check the missing path is stringified into a PLINK command and surfaces
    # as "Failed to open None.bed" / "No such file or directory: 'None.bim'" -
    # naming neither the flag nor the reason. See REFACTOR.md item 23.
    #
    # Deliberately after AncestryArgs is constructed: its __post_init__ rejects
    # --container/--singularity/--cloud, and a flag that can never work is more
    # useful to report than a panel the user could supply.
    if ancestry_args.run_ancestry:
        missing = [
            flag
            for flag, value in (
                ("--ref-panel", ancestry_args.ref_panel),
                ("--ref-labels", ancestry_args.ref_labels),
            )
            if value is None
        ]
        if missing:
            detail = (
                "--model still needs it: inference re-derives the reference "
                "PCA from the panel."
                if ancestry_args.model_path
                else "Training builds the reference PCA from it."
            )
            raise ValueError(
                f"--ancestry requires {' and '.join(missing)}. {detail}"
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
