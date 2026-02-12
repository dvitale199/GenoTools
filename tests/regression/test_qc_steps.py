"""
Regression tests for QC steps.

These tests compare the new implementation against golden reference outputs
to detect any regressions during refactoring.

Run after generating golden files:
    pytest tests/regression/test_qc_steps.py -v
"""

import pytest
from pathlib import Path
import pandas as pd

from .compare import (
    compare_pfiles,
    compare_metrics,
    load_golden_metrics,
)


# =============================================================================
# Sample QC Regression Tests
# =============================================================================


class TestCallrateRegression:
    """Regression tests for callrate filtering."""

    def test_callrate_output_matches_golden(
        self, test_geno_path: Path, golden_callrate: Path, tmp_output_dir: Path
    ):
        """New implementation produces same sample/variant counts as golden."""
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.callrate import filter_callrate
        from genotools.qc.config import CallrateConfig

        # Load input data
        data = GenotypeData.from_path(test_geno_path)

        # Run current implementation
        output_prefix = tmp_output_dir / "callrate_output"
        config = CallrateConfig(mind=0.05)  # Matches --all_sample default
        result = filter_callrate(data, config, output_prefix)

        # Compare to golden
        golden_output = golden_callrate / "output"
        comparison = compare_pfiles(golden_output, result.output.path)

        assert comparison.equal, (
            f"Output differs from golden:\n"
            f"  Sample diff: {comparison.sample_diff}\n"
            f"  Variant diff: {comparison.variant_diff}\n"
            f"  {comparison.message}"
        )

    def test_callrate_metrics_match(
        self, test_geno_path: Path, golden_callrate: Path, tmp_output_dir: Path
    ):
        """Outlier counts match golden."""
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.callrate import filter_callrate
        from genotools.qc.config import CallrateConfig

        # Load golden metrics
        golden_metrics = load_golden_metrics(golden_callrate.parent, "callrate")
        assert golden_metrics is not None, "Golden metrics not found"

        # Load input data
        data = GenotypeData.from_path(test_geno_path)

        # Run current implementation
        output_prefix = tmp_output_dir / "callrate_output"
        config = CallrateConfig(mind=0.05)  # Matches --all_sample default
        result = filter_callrate(data, config, output_prefix)

        # Compare metrics
        expected_outliers = golden_metrics.get("metrics", {}).get("outlier_count", 0)
        actual_outliers = result.metrics.get("outlier_count", 0)

        assert actual_outliers == expected_outliers, (
            f"Outlier count mismatch: expected {expected_outliers}, got {actual_outliers}"
        )

    def test_callrate_outliers_match_manifest(
        self, test_geno_path: Path, outlier_manifest: Path, tmp_output_dir: Path
    ):
        """Detected outliers match the injected callrate outliers."""
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.callrate import filter_callrate
        from genotools.qc.config import CallrateConfig

        # Load expected outliers from manifest
        manifest_df = pd.read_csv(outlier_manifest)
        expected_outliers = set(
            manifest_df[manifest_df["outlier_type"] == "callrate_samples"]["sample_1"]
        )

        # Load input data
        data = GenotypeData.from_path(test_geno_path)

        # Run callrate filtering
        output_prefix = tmp_output_dir / "callrate_output"
        config = CallrateConfig(mind=0.05)  # Matches --all_sample default
        result = filter_callrate(data, config, output_prefix)

        # Get detected outliers
        if result.pruned_samples_file and result.pruned_samples_file.exists():
            outliers_df = pd.read_csv(result.pruned_samples_file, sep="\t")
            detected_outliers = set(outliers_df["IID"])
        else:
            detected_outliers = set()

        # All expected outliers should be detected
        missing = expected_outliers - detected_outliers

        assert len(missing) == 0, f"Expected outliers not detected: {missing}"
        # Note: extra outliers might be acceptable if the data has other issues


class TestSexRegression:
    """Regression tests for sex check."""

    def test_sex_output_matches_golden(
        self, golden_callrate: Path, golden_sex: Path, tmp_output_dir: Path
    ):
        """Sex check produces same output as golden.

        Note: Sex check runs on the callrate-filtered output, not raw data.
        """
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.sex import filter_sex
        from genotools.qc.config import SexConfig

        # Input is the callrate-filtered output (mimics the pipeline)
        callrate_output = golden_callrate / "output"
        data = GenotypeData.from_path(callrate_output)

        # Run current implementation
        output_prefix = tmp_output_dir / "sex_output"
        config = SexConfig()
        result = filter_sex(data, config, output_prefix)

        # Compare to golden
        golden_output = golden_sex / "output"
        comparison = compare_pfiles(golden_output, result.output.path)

        assert comparison.equal, (
            f"Output differs from golden:\n"
            f"  Sample diff: {comparison.sample_diff}\n"
            f"  {comparison.message}"
        )


class TestHetRegression:
    """Regression tests for heterozygosity filtering."""

    def test_het_output_matches_golden(
        self, golden_sex: Path, golden_het: Path, tmp_output_dir: Path
    ):
        """Het filtering produces same output as golden.

        Note: Het check runs on the sex-filtered output.
        """
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.heterozygosity import filter_heterozygosity
        from genotools.qc.config import HetConfig

        # Input is the sex-filtered output (mimics the pipeline)
        sex_output = golden_sex / "output"
        data = GenotypeData.from_path(sex_output)

        # Run current implementation
        output_prefix = tmp_output_dir / "het_output"
        config = HetConfig()
        result = filter_heterozygosity(data, config, output_prefix)

        # Compare to golden
        golden_output = golden_het / "output"
        comparison = compare_pfiles(golden_output, result.output.path)

        assert comparison.equal, (
            f"Output differs from golden:\n"
            f"  Sample diff: {comparison.sample_diff}\n"
            f"  {comparison.message}"
        )


class TestRelatedRegression:
    """Regression tests for relatedness filtering."""

    def test_related_output_matches_golden(
        self, golden_het: Path, golden_related: Path, tmp_output_dir: Path
    ):
        """Relatedness filtering produces same output as golden.

        Note: Relatedness check runs on the het-filtered output, not raw data.
        """
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.relatedness import filter_relatedness
        from genotools.qc.config import RelatedConfig

        # Input is the het-filtered output (mimics the pipeline)
        het_output = golden_het / "output"
        data = GenotypeData.from_path(het_output)

        # Run current implementation
        output_prefix = tmp_output_dir / "related_output"
        config = RelatedConfig(related_cutoff=0.0884)
        result = filter_relatedness(data, config, output_prefix)

        # Compare to golden
        golden_output = golden_related / "output"
        comparison = compare_pfiles(golden_output, result.output.path)

        assert comparison.equal, (
            f"Output differs from golden:\n"
            f"  Sample diff: {comparison.sample_diff}\n"
            f"  {comparison.message}"
        )


class TestDirectOutlierValidation:
    """Tests that validate outlier detection against the known manifest."""

    def test_all_callrate_outliers_detected(
        self, test_geno_path: Path, outlier_manifest: Path, tmp_output_dir: Path
    ):
        """All injected callrate outliers should be detected."""
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.callrate import filter_callrate
        from genotools.qc.config import CallrateConfig

        manifest_df = pd.read_csv(outlier_manifest)
        expected = set(
            manifest_df[manifest_df["outlier_type"] == "callrate_samples"]["sample_1"]
        )

        data = GenotypeData.from_path(test_geno_path)
        output_prefix = tmp_output_dir / "callrate_output"
        config = CallrateConfig(mind=0.05)  # Matches --all_sample default
        result = filter_callrate(data, config, output_prefix)

        if result.pruned_samples_file and result.pruned_samples_file.exists():
            detected = set(pd.read_csv(result.pruned_samples_file, sep="\t")["IID"])
        else:
            detected = set()

        missing = expected - detected
        assert len(missing) == 0, f"Missed callrate outliers: {missing}"

    def test_sex_discordance_detection(
        self, test_geno_path: Path, outlier_manifest: Path, tmp_output_dir: Path
    ):
        """Sex discordant samples should be detected."""
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.sex import filter_sex
        from genotools.qc.config import SexConfig

        manifest_df = pd.read_csv(outlier_manifest)
        expected = set(
            manifest_df[manifest_df["outlier_type"] == "sex_discordance_samples"]["sample_1"]
        )

        data = GenotypeData.from_path(test_geno_path)
        output_prefix = tmp_output_dir / "sex_output"
        config = SexConfig()
        result = filter_sex(data, config, output_prefix)

        if result.pruned_samples_file and result.pruned_samples_file.exists():
            detected = set(pd.read_csv(result.pruned_samples_file, sep="\t")["IID"])
        else:
            detected = set()

        # Check how many expected outliers were detected
        detected_expected = expected & detected
        detection_rate = len(detected_expected) / len(expected) if expected else 0

        # We expect most (>80%) of the injected outliers to be detected
        assert detection_rate >= 0.8, (
            f"Low detection rate for sex outliers: {detection_rate:.1%}. "
            f"Expected {len(expected)}, detected {len(detected_expected)}"
        )

    def test_het_outliers_detection(
        self, test_geno_path: Path, outlier_manifest: Path, tmp_output_dir: Path
    ):
        """Heterozygosity outliers should be detected."""
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.heterozygosity import filter_heterozygosity
        from genotools.qc.config import HetConfig

        manifest_df = pd.read_csv(outlier_manifest)
        expected_high = set(
            manifest_df[manifest_df["outlier_type"] == "het_high_samples"]["sample_1"]
        )
        expected_low = set(
            manifest_df[manifest_df["outlier_type"] == "het_low_samples"]["sample_1"]
        )
        expected = expected_high | expected_low

        data = GenotypeData.from_path(test_geno_path)
        output_prefix = tmp_output_dir / "het_output"
        config = HetConfig()
        result = filter_heterozygosity(data, config, output_prefix)

        if result.pruned_samples_file and result.pruned_samples_file.exists():
            detected = set(pd.read_csv(result.pruned_samples_file, sep="\t")["IID"])
        else:
            detected = set()

        detected_expected = expected & detected
        detection_rate = len(detected_expected) / len(expected) if expected else 0

        # We expect most (>70%) of the injected outliers to be detected
        # Het detection can be more variable due to population statistics
        assert detection_rate >= 0.7, (
            f"Low detection rate for het outliers: {detection_rate:.1%}. "
            f"Expected {len(expected)}, detected {len(detected_expected)}"
        )


# =============================================================================
# Variant QC Regression Tests
# =============================================================================


class TestGenoRegression:
    """Regression tests for variant missingness (geno) filtering."""

    def test_geno_output_matches_golden(
        self, golden_related: Path, golden_geno: Path, tmp_output_dir: Path
    ):
        """Geno filtering produces same output as golden.

        Note: Variant QC runs after sample QC, so input is related-filtered output.
        """
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.variant_missingness import filter_variant_missingness
        from genotools.qc.config import GenoConfig

        # Input is the sample-QC-filtered output
        related_output = golden_related / "output"
        data = GenotypeData.from_path(related_output)

        # Run current implementation
        output_prefix = tmp_output_dir / "geno_output"
        config = GenoConfig(geno=0.05)
        result = filter_variant_missingness(data, config, output_prefix)

        # Compare to golden
        golden_output = golden_geno / "output"
        comparison = compare_pfiles(golden_output, result.output.path)

        assert comparison.equal, (
            f"Output differs from golden:\n"
            f"  Variant diff: {comparison.variant_diff}\n"
            f"  {comparison.message}"
        )


class TestHWERegression:
    """Regression tests for Hardy-Weinberg equilibrium filtering."""

    def test_hwe_output_matches_golden(
        self, golden_geno: Path, golden_hwe: Path, tmp_output_dir: Path
    ):
        """HWE filtering produces same output as golden.

        Note: HWE runs after geno filtering.
        """
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.hwe import filter_hwe
        from genotools.qc.config import HWEConfig

        # Input is the geno-filtered output
        geno_output = golden_geno / "output"
        data = GenotypeData.from_path(geno_output)

        # Run current implementation
        output_prefix = tmp_output_dir / "hwe_output"
        config = HWEConfig(hwe_threshold=1e-4)
        result = filter_hwe(data, config, output_prefix)

        # Compare to golden
        golden_output = golden_hwe / "output"
        comparison = compare_pfiles(golden_output, result.output.path)

        assert comparison.equal, (
            f"Output differs from golden:\n"
            f"  Variant diff: {comparison.variant_diff}\n"
            f"  {comparison.message}"
        )


class TestCaseControlRegression:
    """Regression tests for case-control differential missingness filtering."""

    def test_case_control_output_matches_golden(
        self, golden_hwe: Path, golden_case_control: Path, tmp_output_dir: Path
    ):
        """Case-control filtering produces same output as golden.

        Note: Case-control runs after HWE filtering.
        """
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.case_control import filter_case_control
        from genotools.qc.config import CaseControlConfig

        # Input is the hwe-filtered output
        hwe_output = golden_hwe / "output"
        data = GenotypeData.from_path(hwe_output)

        # Run current implementation
        output_prefix = tmp_output_dir / "case_control_output"
        config = CaseControlConfig(p_threshold=1e-4)
        result = filter_case_control(data, config, output_prefix)

        # Compare to golden
        golden_output = golden_case_control / "output"
        comparison = compare_pfiles(golden_output, result.output.path)

        assert comparison.equal, (
            f"Output differs from golden:\n"
            f"  Variant diff: {comparison.variant_diff}\n"
            f"  {comparison.message}"
        )


class TestVariantOutlierValidation:
    """Tests that validate variant outlier detection against the known manifest."""

    def test_geno_miss_variants_detected(
        self, golden_related: Path, outlier_manifest: Path, tmp_output_dir: Path
    ):
        """Variants with high missingness should be filtered."""
        from genotools.core.genotypes import GenotypeData
        from genotools.qc.steps.variant_missingness import filter_variant_missingness
        from genotools.qc.config import GenoConfig

        manifest_df = pd.read_csv(outlier_manifest)
        # Count expected variant outliers (we injected 150 geno_miss_variants)
        expected_count = len(
            manifest_df[manifest_df["outlier_type"] == "geno_miss_variants"]
        )

        # Run geno filtering on related output
        related_output = golden_related / "output"
        data = GenotypeData.from_path(related_output)
        output_prefix = tmp_output_dir / "geno_output"
        config = GenoConfig(geno=0.05)
        result = filter_variant_missingness(data, config, output_prefix)

        # Check that we filtered some variants
        pruned_count = result.metrics.get("geno_removed_count", 0)

        # We should detect at least 50% of the injected outliers
        # (threshold interactions may cause some to pass)
        assert pruned_count >= expected_count * 0.5, (
            f"Low detection: expected ~{expected_count} bad variants, "
            f"only filtered {pruned_count}"
        )


# =============================================================================
# QCPipeline API Regression Tests (New Implementation)
# =============================================================================


class TestQCPipelineSampleQC:
    """Regression tests for QCPipeline.sample_qc() against golden output.

    These tests verify the NEW QCPipeline API produces the same results
    as the old SampleQC class-based implementation.
    """

    def test_sample_qc_pipeline_matches_golden(
        self, test_geno_path: Path, golden_all_sample: Path, tmp_output_dir: Path
    ):
        """QCPipeline.sample_qc() produces same final output as golden.

        This is the key regression test for Phase 2 completion.
        """
        from genotools.core.genotypes import GenotypeData
        from genotools.qc import (
            QCPipeline,
            CallrateConfig,
            SexConfig,
            HetConfig,
            RelatedConfig,
        )

        # Load input data
        data = GenotypeData.from_path(test_geno_path)

        # Create pipeline with same settings as --all_sample
        pipeline = QCPipeline.sample_qc(
            callrate_config=CallrateConfig(mind=0.05),
            sex_config=SexConfig(),
            het_config=HetConfig(),
            related_config=RelatedConfig(related_cutoff=0.0884),
        )

        # Run the pipeline
        output_dir = tmp_output_dir / "sample_qc_pipeline"
        result = pipeline.run(data, output_dir)

        # Compare final output to golden
        golden_output = golden_all_sample / "output"
        comparison = compare_pfiles(golden_output, result.output.path)

        assert comparison.equal, (
            f"QCPipeline output differs from golden:\n"
            f"  Sample diff: {comparison.sample_diff}\n"
            f"  Variant diff: {comparison.variant_diff}\n"
            f"  {comparison.message}"
        )

    def test_sample_qc_pipeline_sample_counts_match(
        self, test_geno_path: Path, golden_all_sample: Path, tmp_output_dir: Path
    ):
        """QCPipeline removes the same number of samples as golden."""
        from genotools.core.genotypes import GenotypeData
        from genotools.qc import (
            QCPipeline,
            CallrateConfig,
            SexConfig,
            HetConfig,
            RelatedConfig,
        )

        # Load golden output to get expected counts
        golden_output = golden_all_sample / "output"
        golden_data = GenotypeData.from_path(golden_output)

        # Load input and run pipeline
        data = GenotypeData.from_path(test_geno_path)
        pipeline = QCPipeline.sample_qc(
            callrate_config=CallrateConfig(mind=0.05),
            sex_config=SexConfig(),
            het_config=HetConfig(),
            related_config=RelatedConfig(related_cutoff=0.0884),
        )

        output_dir = tmp_output_dir / "sample_qc_counts"
        result = pipeline.run(data, output_dir)

        assert result.output.sample_count == golden_data.sample_count, (
            f"Sample count mismatch: "
            f"expected {golden_data.sample_count}, got {result.output.sample_count}"
        )

    def test_sample_qc_pipeline_is_deterministic(
        self, test_geno_path: Path, tmp_output_dir: Path
    ):
        """QCPipeline.sample_qc() produces identical output on repeated runs."""
        from genotools.core.genotypes import GenotypeData
        from genotools.qc import (
            QCPipeline,
            CallrateConfig,
            SexConfig,
            HetConfig,
            RelatedConfig,
        )

        data = GenotypeData.from_path(test_geno_path)
        pipeline = QCPipeline.sample_qc(
            callrate_config=CallrateConfig(mind=0.05),
            sex_config=SexConfig(),
            het_config=HetConfig(),
            related_config=RelatedConfig(related_cutoff=0.0884),
        )

        # Run twice
        result1 = pipeline.run(data, tmp_output_dir / "run1")
        result2 = pipeline.run(data, tmp_output_dir / "run2")

        comparison = compare_pfiles(result1.output.path, result2.output.path)

        assert comparison.equal, (
            f"Pipeline is not deterministic:\n"
            f"  Sample diff: {comparison.sample_diff}\n"
            f"  {comparison.message}"
        )


class TestQCPipelineVariantQC:
    """Regression tests for QCPipeline.variant_qc() against golden output."""

    def test_variant_qc_pipeline_matches_golden(
        self, golden_related: Path, golden_all_variant: Path, tmp_output_dir: Path
    ):
        """QCPipeline.variant_qc() produces same final output as golden.

        Note: Variant QC runs on sample-QC-filtered data (related output).
        """
        from genotools.core.genotypes import GenotypeData
        from genotools.qc import (
            QCPipeline,
            GenoConfig,
            CaseControlConfig,
            HaplotypeConfig,
            HWEConfig,
        )

        # Input is the sample-QC-filtered output
        related_output = golden_related / "output"
        data = GenotypeData.from_path(related_output)

        # Create pipeline with same settings as --all_variant
        pipeline = QCPipeline.variant_qc(
            geno_config=GenoConfig(geno=0.05),
            hwe_config=HWEConfig(hwe_threshold=1e-4),
            case_control_config=CaseControlConfig(p_threshold=1e-4),
            haplotype_config=HaplotypeConfig(),
        )

        # Run the pipeline
        output_dir = tmp_output_dir / "variant_qc_pipeline"
        result = pipeline.run(data, output_dir)

        # Compare final output to golden
        golden_output = golden_all_variant / "output"
        comparison = compare_pfiles(golden_output, result.output.path)

        assert comparison.equal, (
            f"QCPipeline output differs from golden:\n"
            f"  Sample diff: {comparison.sample_diff}\n"
            f"  Variant diff: {comparison.variant_diff}\n"
            f"  {comparison.message}"
        )


# =============================================================================
# CLI Pipeline Integration Tests (Old CLI - will pass after Phase 4)
# =============================================================================


class TestAllSamplePipeline:
    """Integration tests for the --all_sample pipeline."""

    def test_all_sample_final_output_matches_golden(
        self, test_geno_path: Path, golden_all_sample: Path, tmp_output_dir: Path
    ):
        """Running --all_sample produces same final output as golden.

        This tests the pipeline orchestration in execute_pipeline().
        """
        import subprocess
        import sys

        output_prefix = tmp_output_dir / "all_sample_output"

        # Run the CLI
        cmd = [
            sys.executable, "-m", "genotools",
            "--pfile", str(test_geno_path),
            "--out", str(output_prefix),
            "--all-sample"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        assert result.returncode == 0, (
            f"Pipeline failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )

        # Compare final output to golden
        golden_output = golden_all_sample / "output"
        comparison = compare_pfiles(golden_output, output_prefix)

        assert comparison.equal, (
            f"Pipeline output differs from golden:\n"
            f"  Sample diff: {comparison.sample_diff}\n"
            f"  Variant diff: {comparison.variant_diff}\n"
            f"  {comparison.message}"
        )

    def test_all_sample_is_deterministic(
        self, test_geno_path: Path, tmp_output_dir: Path
    ):
        """Running --all_sample twice produces identical output.

        This ensures the pipeline is deterministic.
        """
        import subprocess
        import sys

        output_prefix_1 = tmp_output_dir / "all_sample_run1"
        output_prefix_2 = tmp_output_dir / "all_sample_run2"

        # Run the CLI twice
        for output_prefix in [output_prefix_1, output_prefix_2]:
            cmd = [
                sys.executable, "-m", "genotools",
                "--pfile", str(test_geno_path),
                "--out", str(output_prefix),
                "--all-sample"
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            assert result.returncode == 0, (
                f"Pipeline failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
            )

        # Compare the two runs
        comparison = compare_pfiles(output_prefix_1, output_prefix_2)

        assert comparison.equal, (
            f"Pipeline is not deterministic:\n"
            f"  Sample diff: {comparison.sample_diff}\n"
            f"  {comparison.message}"
        )


class TestAllVariantPipeline:
    """Integration tests for the --all_variant pipeline."""

    def test_all_variant_final_output_matches_golden(
        self, golden_related: Path, golden_all_variant: Path, tmp_output_dir: Path
    ):
        """Running --all_variant produces same final output as golden.

        Note: Variant QC runs on sample-QC-filtered data.
        """
        import subprocess
        import sys

        # Use related output as input (after sample QC)
        input_path = golden_related / "output"
        output_prefix = tmp_output_dir / "all_variant_output"

        # Run the CLI
        cmd = [
            sys.executable, "-m", "genotools",
            "--pfile", str(input_path),
            "--out", str(output_prefix),
            "--all-variant"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)

        assert result.returncode == 0, (
            f"Pipeline failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )

        # Compare final output to golden
        golden_output = golden_all_variant / "output"
        comparison = compare_pfiles(golden_output, output_prefix)

        assert comparison.equal, (
            f"Pipeline output differs from golden:\n"
            f"  Variant diff: {comparison.variant_diff}\n"
            f"  {comparison.message}"
        )


class TestFullPipeline:
    """Integration tests for running both sample and variant QC."""

    def test_full_qc_pipeline_is_deterministic(
        self, test_geno_path: Path, tmp_output_dir: Path
    ):
        """Running --all_sample --all_variant twice produces identical output.

        This ensures the full pipeline is deterministic.
        """
        import subprocess
        import sys

        output_prefix_1 = tmp_output_dir / "full_qc_run1"
        output_prefix_2 = tmp_output_dir / "full_qc_run2"

        # Run the CLI twice
        for output_prefix in [output_prefix_1, output_prefix_2]:
            cmd = [
                sys.executable, "-m", "genotools",
                "--pfile", str(test_geno_path),
                "--out", str(output_prefix),
                "--all-sample",
                "--all-variant"
            ]
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
            assert result.returncode == 0, (
                f"Pipeline failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
            )

        # Compare the two runs
        comparison = compare_pfiles(output_prefix_1, output_prefix_2)

        assert comparison.equal, (
            f"Full pipeline is not deterministic:\n"
            f"  Sample diff: {comparison.sample_diff}\n"
            f"  Variant diff: {comparison.variant_diff}\n"
            f"  {comparison.message}"
        )

    def test_full_qc_pipeline_completes(
        self, test_geno_path: Path, tmp_output_dir: Path
    ):
        """Running --all_sample --all_variant completes successfully.

        This is a smoke test for the full pipeline.
        """
        import subprocess
        import sys

        output_prefix = tmp_output_dir / "full_qc_output"

        cmd = [
            sys.executable, "-m", "genotools",
            "--pfile", str(test_geno_path),
            "--out", str(output_prefix),
            "--all-sample",
            "--all-variant"
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        assert result.returncode == 0, (
            f"Pipeline failed:\nstdout: {result.stdout}\nstderr: {result.stderr}"
        )

        # Verify output files exist
        assert (output_prefix.with_suffix(".pgen")).exists(), "Output pgen not created"
        assert (output_prefix.with_suffix(".psam")).exists(), "Output psam not created"
        assert (output_prefix.with_suffix(".pvar")).exists(), "Output pvar not created"
        assert (output_prefix.with_suffix(".json")).exists(), "Output json not created"
