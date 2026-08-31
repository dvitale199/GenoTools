#!/usr/bin/env python3
"""
Generate golden reference files using the current GenoTools implementation.

Run this ONCE before starting refactoring to capture expected outputs.
The golden files serve as the reference for regression testing.

Usage:
    python tests/scripts/generate_golden.py \
        --geno tests/data/synthetic/genotools_test \
        --out tests/regression/golden

Requirements:
    - Synthetic test data must exist (run generate_test_data.py first)
    - GenoTools must be installed (pip install -e .)
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

# Add project root to path for imports
PROJECT_ROOT = Path(__file__).parent.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

from genotools.core.genotypes import GenotypeData
from genotools.qc import (
    filter_callrate, CallrateConfig,
    filter_sex, SexConfig,
    filter_heterozygosity, HetConfig,
    filter_relatedness, RelatedConfig,
    filter_variant_missingness, GenoConfig,
    filter_hwe, HWEConfig,
    filter_case_control, CaseControlConfig,
    filter_haplotype, HaplotypeConfig,
    prune_ld, LDConfig,
)


def generate_callrate_golden(geno_path: Path, out_dir: Path, mind: float = 0.05):
    """Generate golden files for callrate filtering."""
    step_out = out_dir / "callrate"
    step_out.mkdir(parents=True, exist_ok=True)
    output_prefix = step_out / "output"

    data = GenotypeData.from_path(geno_path)
    config = CallrateConfig(mind=mind)
    result = filter_callrate(data, config, output_prefix)

    result_dict = result.to_dict()
    with open(step_out / "metrics.json", "w") as f:
        json.dump(result_dict, f, indent=2, default=str)

    return result_dict


def generate_sex_golden(geno_path: Path, out_dir: Path):
    """Generate golden files for sex check."""
    step_out = out_dir / "sex"
    step_out.mkdir(parents=True, exist_ok=True)
    output_prefix = step_out / "output"

    data = GenotypeData.from_path(geno_path)
    config = SexConfig()
    result = filter_sex(data, config, output_prefix)

    result_dict = result.to_dict()
    with open(step_out / "metrics.json", "w") as f:
        json.dump(result_dict, f, indent=2, default=str)

    return result_dict


def generate_het_golden(geno_path: Path, out_dir: Path):
    """Generate golden files for heterozygosity filtering."""
    step_out = out_dir / "het"
    step_out.mkdir(parents=True, exist_ok=True)
    output_prefix = step_out / "output"

    data = GenotypeData.from_path(geno_path)
    config = HetConfig()
    result = filter_heterozygosity(data, config, output_prefix)

    result_dict = result.to_dict()
    with open(step_out / "metrics.json", "w") as f:
        json.dump(result_dict, f, indent=2, default=str)

    return result_dict


def generate_related_golden(geno_path: Path, out_dir: Path, kinship: float = 0.0884):
    """Generate golden files for relatedness filtering."""
    step_out = out_dir / "related"
    step_out.mkdir(parents=True, exist_ok=True)
    output_prefix = step_out / "output"

    data = GenotypeData.from_path(geno_path)
    config = RelatedConfig(related_cutoff=kinship)
    result = filter_relatedness(data, config, output_prefix)

    result_dict = result.to_dict()
    with open(step_out / "metrics.json", "w") as f:
        json.dump(result_dict, f, indent=2, default=str)

    return result_dict


def generate_geno_golden(geno_path: Path, out_dir: Path, geno_threshold: float = 0.05):
    """Generate golden files for variant missingness filtering."""
    step_out = out_dir / "geno"
    step_out.mkdir(parents=True, exist_ok=True)
    output_prefix = step_out / "output"

    data = GenotypeData.from_path(geno_path)
    config = GenoConfig(geno=geno_threshold)
    result = filter_variant_missingness(data, config, output_prefix)

    result_dict = result.to_dict()
    with open(step_out / "metrics.json", "w") as f:
        json.dump(result_dict, f, indent=2, default=str)

    return result_dict


def generate_hwe_golden(geno_path: Path, out_dir: Path, hwe_threshold: float = 1e-4):
    """Generate golden files for HWE filtering."""
    step_out = out_dir / "hwe"
    step_out.mkdir(parents=True, exist_ok=True)
    output_prefix = step_out / "output"

    data = GenotypeData.from_path(geno_path)
    config = HWEConfig(hwe_threshold=hwe_threshold)
    result = filter_hwe(data, config, output_prefix)

    result_dict = result.to_dict()
    with open(step_out / "metrics.json", "w") as f:
        json.dump(result_dict, f, indent=2, default=str)

    return result_dict


def generate_case_control_golden(geno_path: Path, out_dir: Path, p_threshold: float = 1e-4):
    """Generate golden files for case-control differential missingness filtering."""
    step_out = out_dir / "case_control"
    step_out.mkdir(parents=True, exist_ok=True)
    output_prefix = step_out / "output"

    data = GenotypeData.from_path(geno_path)
    config = CaseControlConfig(p_threshold=p_threshold)
    result = filter_case_control(data, config, output_prefix)

    result_dict = result.to_dict()
    with open(step_out / "metrics.json", "w") as f:
        json.dump(result_dict, f, indent=2, default=str)

    return result_dict


def generate_haplotype_golden(geno_path: Path, out_dir: Path, p_threshold: float = 1e-4):
    """Generate golden files for haplotype filtering."""
    step_out = out_dir / "haplotype"
    step_out.mkdir(parents=True, exist_ok=True)
    output_prefix = step_out / "output"

    data = GenotypeData.from_path(geno_path)
    config = HaplotypeConfig(p_threshold=p_threshold)
    result = filter_haplotype(data, config, output_prefix)

    result_dict = result.to_dict()
    with open(step_out / "metrics.json", "w") as f:
        json.dump(result_dict, f, indent=2, default=str)

    return result_dict


def generate_ld_golden(geno_path: Path, out_dir: Path, window_size: int = 50, step_size: int = 5, r2_threshold: float = 0.5):
    """Generate golden files for LD pruning."""
    step_out = out_dir / "ld"
    step_out.mkdir(parents=True, exist_ok=True)
    output_prefix = step_out / "output"

    data = GenotypeData.from_path(geno_path)
    config = LDConfig(window_size=window_size, step_size=step_size, r2_threshold=r2_threshold)
    result = prune_ld(data, config, output_prefix)

    result_dict = result.to_dict()
    with open(step_out / "metrics.json", "w") as f:
        json.dump(result_dict, f, indent=2, default=str)

    return result_dict


def generate_all_sample_golden(geno_path: Path, out_dir: Path):
    """Generate golden files for --all-sample pipeline run."""
    step_out = out_dir / "all_sample"
    step_out.mkdir(parents=True, exist_ok=True)
    output_prefix = step_out / "output"

    # Clean up any existing output files (pipeline refuses to overwrite)
    for pattern in ["output*", "*.log", "*.json"]:
        for f in step_out.glob(pattern):
            # A pipeline run can leave a directory behind (--all-variant writes
            # an output_bfile/), and unlink() raises IsADirectoryError on one -
            # after having already deleted the goldens ahead of it in the glob.
            if f.is_dir():
                shutil.rmtree(f)
            else:
                f.unlink()

    cmd = [
        sys.executable, "-m", "genotools",
        "--pfile", str(geno_path),
        "--out", str(output_prefix),
        "--all-sample"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        # Recording the failure as a golden is worse than not regenerating:
        # the pipeline goldens went stale that way, when the command here was
        # "-m genotools.cli" (the package, not the entry point) and every run
        # since exited with "cannot be directly executed".
        raise SystemExit(
            f"Pipeline run failed (exit {result.returncode}); refusing to "
            f"write a golden from it:\n{result.stderr}"
        )

    metrics = {
        "pass": result.returncode == 0,
        "step": "all_sample_pipeline",
        "returncode": result.returncode,
        "stdout": result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout,
        "stderr": result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr,
    }

    with open(step_out / "metrics.json", "w") as f:
        json.dump(metrics, f, indent=2, default=str)

    return {"pass": result.returncode == 0, "output": {"plink_out": str(output_prefix)}}


def generate_all_variant_golden(geno_path: Path, out_dir: Path):
    """Generate golden files for --all-variant pipeline run."""
    step_out = out_dir / "all_variant"
    step_out.mkdir(parents=True, exist_ok=True)
    output_prefix = step_out / "output"

    # Clean up any existing output files (pipeline refuses to overwrite)
    for pattern in ["output*", "*.log", "*.json"]:
        for f in step_out.glob(pattern):
            # A pipeline run can leave a directory behind (--all-variant writes
            # an output_bfile/), and unlink() raises IsADirectoryError on one -
            # after having already deleted the goldens ahead of it in the glob.
            if f.is_dir():
                shutil.rmtree(f)
            else:
                f.unlink()

    cmd = [
        sys.executable, "-m", "genotools",
        "--pfile", str(geno_path),
        "--out", str(output_prefix),
        "--all-variant"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if result.returncode != 0:
        # Recording the failure as a golden is worse than not regenerating:
        # the pipeline goldens went stale that way, when the command here was
        # "-m genotools.cli" (the package, not the entry point) and every run
        # since exited with "cannot be directly executed".
        raise SystemExit(
            f"Pipeline run failed (exit {result.returncode}); refusing to "
            f"write a golden from it:\n{result.stderr}"
        )

    metrics = {
        "pass": result.returncode == 0,
        "step": "all_variant_pipeline",
        "returncode": result.returncode,
        "stdout": result.stdout[-2000:] if len(result.stdout) > 2000 else result.stdout,
        "stderr": result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr,
    }

    with open(step_out / "metrics.json", "w") as f:
        json.dump(metrics, f, indent=2, default=str)

    return {"pass": result.returncode == 0, "output": {"plink_out": str(output_prefix)}}


def main():
    parser = argparse.ArgumentParser(
        description="Generate golden reference files for regression testing"
    )
    parser.add_argument(
        "--geno",
        type=Path,
        required=True,
        help="Test genotype file prefix (e.g., tests/data/synthetic/genotools_test)"
    )
    parser.add_argument(
        "--out",
        type=Path,
        required=True,
        help="Output directory for golden files (e.g., tests/regression/golden)"
    )
    parser.add_argument(
        "--steps",
        nargs="+",
        default=["callrate", "sex", "het", "related", "geno", "hwe", "case_control", "haplotype", "ld", "all_sample", "all_variant"],
        help="QC steps to generate golden files for"
    )
    args = parser.parse_args()

    # Validate input
    if not args.geno.with_suffix(".pgen").exists():
        print(f"ERROR: Input pfile not found: {args.geno}")
        print("Run generate_test_data.py first to create synthetic test data.")
        sys.exit(1)

    # Create output directory
    args.out.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Golden File Generator")
    print("=" * 60)
    print(f"\nInput: {args.geno}")
    print(f"Output: {args.out}")
    print(f"Steps: {', '.join(args.steps)}")

    # Track results
    results = {}

    # Run each step
    print("\nGenerating golden files...")

    current_input = args.geno

    # Track sample-QC output for pipeline tests
    sample_qc_output = None

    if "callrate" in args.steps:
        print("\n[1/11] Callrate filtering...")
        results["callrate"] = generate_callrate_golden(current_input, args.out)
        if results["callrate"]["pass"]:
            current_input = Path(results["callrate"]["output"]["plink_out"])

    if "sex" in args.steps:
        print("\n[2/11] Sex check...")
        results["sex"] = generate_sex_golden(current_input, args.out)
        if results["sex"]["pass"]:
            current_input = Path(results["sex"]["output"]["plink_out"])

    if "het" in args.steps:
        print("\n[3/11] Heterozygosity filtering...")
        results["het"] = generate_het_golden(current_input, args.out)
        if results["het"]["pass"]:
            current_input = Path(results["het"]["output"]["plink_out"])

    if "related" in args.steps:
        print("\n[4/11] Relatedness filtering...")
        results["related"] = generate_related_golden(current_input, args.out)
        if results["related"]["pass"]:
            current_input = Path(results["related"]["output"]["plink_out"])
            sample_qc_output = current_input  # Save for pipeline tests

    if "geno" in args.steps:
        print("\n[5/11] Variant missingness filtering...")
        results["geno"] = generate_geno_golden(current_input, args.out)
        if results["geno"]["pass"]:
            current_input = Path(results["geno"]["output"]["plink_out"])

    if "hwe" in args.steps:
        print("\n[6/11] HWE filtering...")
        results["hwe"] = generate_hwe_golden(current_input, args.out)
        if results["hwe"]["pass"]:
            current_input = Path(results["hwe"]["output"]["plink_out"])

    if "case_control" in args.steps:
        print("\n[7/11] Case-control filtering...")
        results["case_control"] = generate_case_control_golden(current_input, args.out)
        if results["case_control"]["pass"]:
            current_input = Path(results["case_control"]["output"]["plink_out"])

    if "haplotype" in args.steps:
        print("\n[8/11] Haplotype filtering...")
        results["haplotype"] = generate_haplotype_golden(current_input, args.out)
        if results["haplotype"]["pass"]:
            current_input = Path(results["haplotype"]["output"]["plink_out"])

    if "ld" in args.steps:
        print("\n[9/11] LD pruning...")
        results["ld"] = generate_ld_golden(current_input, args.out)
        if results["ld"]["pass"]:
            current_input = Path(results["ld"]["output"]["plink_out"])

    # Pipeline golden files (run from raw input)
    if "all_sample" in args.steps:
        print("\n[10/11] --all-sample pipeline...")
        results["all_sample"] = generate_all_sample_golden(args.geno, args.out)

    if "all_variant" in args.steps:
        print("\n[11/11] --all-variant pipeline...")
        # Use sample-QC-filtered output if available, otherwise use related output
        all_variant_input = sample_qc_output if sample_qc_output else (args.out / "related" / "output")
        if all_variant_input and Path(str(all_variant_input) + ".pgen").exists():
            results["all_variant"] = generate_all_variant_golden(all_variant_input, args.out)
        else:
            print("    SKIP: No sample-QC output available for all_variant")

    # Save summary
    summary_path = args.out / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2, default=str)

    print("\n" + "=" * 60)
    print("DONE!")
    print("=" * 60)

    print("\nGenerated golden files for:")
    for step, result in results.items():
        status = "PASS" if result.get("pass", False) else "FAIL"
        print(f"  {step}: {status}")

    print(f"\nSummary saved to: {summary_path}")


if __name__ == "__main__":
    main()
