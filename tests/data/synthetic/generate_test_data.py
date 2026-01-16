#!/usr/bin/env python3
"""
Generate synthetic test data for GenoTools regression testing.

This script creates a synthetic genotype dataset with intentional outliers
to test each QC step in the GenoTools pipeline.

Usage:
    python generate_test_data.py [--output-dir OUTPUT_DIR]

Output:
    - genotools_test.pgen/pvar/psam: Synthetic genotype data
    - outlier_manifest.csv: Records of all injected outliers
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pgenlib

# =============================================================================
# Configuration
# =============================================================================

RANDOM_SEED = 42
N_SAMPLES = 500
N_AUTOSOMAL_VARIANTS = 40000
N_X_VARIANTS = 500
N_TOTAL_VARIANTS = N_AUTOSOMAL_VARIANTS + N_X_VARIANTS

# Reference panel location
REF_BIM = Path.home() / ".genotools/ref/ref_panel/ref_panel_gp2_prune_rm_underperform_pos_update.bim"

# X chromosome valid range for sex check
X_START = 2781479
X_END = 155701383

# Outlier counts
OUTLIER_CONFIG = {
    'callrate': 12,           # Samples with >10% missingness
    'sex_discordance': 12,    # Samples with wrong genetic sex
    'het_high': 4,            # Samples with excess heterozygosity
    'het_low': 4,             # Samples with deficient heterozygosity
    'duplicate_pairs': 4,     # Duplicate sample pairs
    'related_pairs': 8,       # Related sample pairs (1st/2nd degree)
    'geno_miss_variants': 150,      # Variants with >5% missingness
    'case_control_variants': 75,    # Variants with differential missingness
    'hwe_violation_variants': 75,   # Variants violating HWE
}


# =============================================================================
# Utility Functions
# =============================================================================

def set_seed(offset=0):
    """Set numpy random seed with optional offset."""
    np.random.seed(RANDOM_SEED + offset)


def run_command(cmd, description=""):
    """Run a shell command and check for errors."""
    print(f"  Running: {description or cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  ERROR: {result.stderr}")
        raise RuntimeError(f"Command failed: {cmd}")
    return result


# =============================================================================
# Step 1: Extract Variants from Reference Panel
# =============================================================================

def extract_autosomal_variants(output_dir):
    """
    Extract N_AUTOSOMAL_VARIANTS from reference panel BIM file,
    stratified by chromosome to maintain distribution.
    """
    print("\nStep 1: Extracting variants from reference panel...")

    if not REF_BIM.exists():
        raise FileNotFoundError(f"Reference BIM not found: {REF_BIM}")

    set_seed(0)

    # Read reference BIM
    ref_bim = pd.read_csv(
        REF_BIM,
        sep='\t',
        header=None,
        names=['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2'],
        dtype={'CHR': str}
    )

    print(f"  Reference panel has {len(ref_bim)} variants")

    # Sample proportionally from each chromosome
    sampling_fraction = N_AUTOSOMAL_VARIANTS / len(ref_bim)
    sampled = ref_bim.groupby('CHR', group_keys=False).apply(
        lambda x: x.sample(frac=sampling_fraction, random_state=RANDOM_SEED)
    )

    # Ensure exactly N_AUTOSOMAL_VARIANTS
    if len(sampled) > N_AUTOSOMAL_VARIANTS:
        sampled = sampled.head(N_AUTOSOMAL_VARIANTS)
    elif len(sampled) < N_AUTOSOMAL_VARIANTS:
        # Sample additional variants if needed
        remaining = ref_bim[~ref_bim['SNP'].isin(sampled['SNP'])]
        additional = remaining.sample(n=N_AUTOSOMAL_VARIANTS - len(sampled), random_state=RANDOM_SEED)
        sampled = pd.concat([sampled, additional])

    # Sort by chromosome and position
    sampled['CHR_NUM'] = sampled['CHR'].astype(int)
    sampled = sampled.sort_values(['CHR_NUM', 'POS']).drop(columns=['CHR_NUM'])
    sampled = sampled.reset_index(drop=True)

    # Save variant list
    variant_list_path = output_dir / "autosomal_variants.txt"
    sampled['SNP'].to_csv(variant_list_path, index=False, header=False)

    print(f"  Selected {len(sampled)} autosomal variants")
    return sampled


# =============================================================================
# Step 2: Generate Synthetic X Chromosome Variants
# =============================================================================

def generate_x_variants(output_dir):
    """
    Generate synthetic X chromosome variant definitions for sex checks.
    """
    print("\nStep 2: Generating X chromosome variants...")

    set_seed(1)

    # Generate positions across X chromosome
    x_positions = np.sort(np.random.randint(X_START, X_END, N_X_VARIANTS))

    # Generate alleles
    alleles = ['A', 'C', 'G', 'T']
    a1 = np.random.choice(alleles, N_X_VARIANTS)
    a2 = np.random.choice(alleles, N_X_VARIANTS)

    # Ensure A1 != A2
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    for i in range(N_X_VARIANTS):
        if a1[i] == a2[i]:
            a2[i] = complement[a1[i]]

    x_variants = pd.DataFrame({
        'CHR': 'X',  # Use 'X' instead of '23' to match standard PLINK/VCF conventions
        'SNP': [f'X_SIM_{i:05d}' for i in range(N_X_VARIANTS)],
        'CM': 0,
        'POS': x_positions,
        'A1': a1,
        'A2': a2
    })

    print(f"  Generated {len(x_variants)} X chromosome variants")
    return x_variants


# =============================================================================
# Step 3: Generate Base Genotypes with PLINK2
# =============================================================================

def generate_base_genotypes(output_dir):
    """
    Use PLINK2 --dummy to generate initial random genotypes.
    """
    print("\nStep 3: Generating base genotypes with PLINK2...")

    raw_prefix = output_dir / "genotools_test_raw"

    cmd = (
        f"plink2 --dummy {N_SAMPLES} {N_TOTAL_VARIANTS} 0.01 "
        f"--seed {RANDOM_SEED} "
        f"--out {raw_prefix}"
    )

    run_command(cmd, "PLINK2 --dummy")

    print(f"  Generated {N_SAMPLES} samples x {N_TOTAL_VARIANTS} variants")
    return raw_prefix


# =============================================================================
# Step 4: Create Combined Variant File (PVAR)
# =============================================================================

def create_variant_file(autosomal_variants, x_variants, output_dir):
    """
    Create PVAR file with reference panel IDs for autosomal
    and synthetic IDs for X chromosome.
    """
    print("\nStep 4: Creating variant file with reference panel IDs...")

    # Combine autosomal and X variants
    all_variants = pd.concat([autosomal_variants, x_variants], ignore_index=True)

    # Create PVAR format
    pvar = pd.DataFrame({
        '#CHROM': all_variants['CHR'],
        'POS': all_variants['POS'],
        'ID': all_variants['SNP'],
        'REF': all_variants['A2'],  # A2 is typically REF in PLINK
        'ALT': all_variants['A1']   # A1 is typically ALT
    })

    # Save PVAR
    pvar_path = output_dir / "genotools_test.pvar"
    with open(pvar_path, 'w') as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\n")
    pvar.to_csv(pvar_path, sep='\t', index=False, header=False, mode='a')

    print(f"  Created PVAR with {len(pvar)} variants")
    return pvar_path, all_variants


# =============================================================================
# Step 5: Create Sample Metadata (PSAM)
# =============================================================================

def create_sample_metadata(output_dir):
    """
    Generate PSAM file with proper sample structure.
    """
    print("\nStep 5: Creating sample metadata...")

    set_seed(2)

    # Create sample IDs
    sample_ids = [f'SAMP_{i:04d}' for i in range(N_SAMPLES)]
    family_ids = [f'FAM_{i:04d}' for i in range(N_SAMPLES)]

    # Assign sex (roughly 50/50)
    sex = np.zeros(N_SAMPLES, dtype=int)
    sex[:N_SAMPLES // 2] = 1  # Males
    sex[N_SAMPLES // 2:] = 2  # Females
    np.random.shuffle(sex)

    # Assign phenotype (250 cases, 250 controls)
    pheno = np.zeros(N_SAMPLES, dtype=int)
    pheno[:250] = 2  # Cases
    pheno[250:] = 1  # Controls
    np.random.shuffle(pheno)

    psam = pd.DataFrame({
        '#FID': family_ids,
        'IID': sample_ids,
        'PAT': 0,
        'MAT': 0,
        'SEX': sex,
        'PHENO1': pheno
    })

    # Save PSAM
    psam_path = output_dir / "genotools_test.psam"
    psam.to_csv(psam_path, sep='\t', index=False)

    print(f"  Created PSAM with {len(psam)} samples")
    print(f"    Males: {sum(sex == 1)}, Females: {sum(sex == 2)}")
    print(f"    Cases: {sum(pheno == 2)}, Controls: {sum(pheno == 1)}")

    return psam_path, psam


# =============================================================================
# Step 6: Inject Outliers
# =============================================================================

def load_genotypes(pgen_path, n_samples, n_variants):
    """Load genotypes from PGEN file into numpy array."""
    genotypes = np.zeros((n_variants, n_samples), dtype=np.int8)

    with pgenlib.PgenReader(str(pgen_path).encode()) as reader:
        for var_idx in range(n_variants):
            buf = np.zeros(n_samples, dtype=np.int32)
            reader.read(var_idx, buf)
            genotypes[var_idx] = buf.astype(np.int8)

    return genotypes


def save_genotypes_pgenlib(genotypes, output_prefix, pvar_path, psam_path):
    """
    Save genotypes directly using pgenlib's PgenWriter.

    This avoids the VCF intermediate which had conversion issues.
    Uses append_biallelic which properly handles missing values as -9 in int8.
    """
    n_variants, n_samples = genotypes.shape

    print(f"    Writing pgen file with {n_variants} variants and {n_samples} samples...")

    # Ensure int8 for append_biallelic (handles -9 as missing)
    geno_for_write = genotypes.astype(np.int8)

    pgen_path = str(output_prefix) + '.pgen'

    # Write using PgenWriter with append_biallelic for proper missing handling
    with pgenlib.PgenWriter(
        filename=pgen_path.encode(),
        sample_ct=n_samples,
        variant_ct=n_variants,
        nonref_flags=False
    ) as writer:
        for var_idx in range(n_variants):
            writer.append_biallelic(geno_for_write[var_idx])

            # Progress indicator
            if (var_idx + 1) % 10000 == 0:
                print(f"      Written {var_idx + 1}/{n_variants} variants...")

    # Copy pvar and psam to output location if different from source
    import shutil
    output_pvar = Path(str(output_prefix) + '.pvar')
    output_psam = Path(str(output_prefix) + '.psam')

    if Path(pvar_path).resolve() != output_pvar.resolve():
        shutil.copy(pvar_path, output_pvar)
    if Path(psam_path).resolve() != output_psam.resolve():
        shutil.copy(psam_path, output_psam)

    print(f"    Pgen files created at: {output_prefix}")
    return output_prefix


def inject_callrate_outliers(genotypes, outlier_manifest, used_samples):
    """
    Inject high missingness into selected samples (>10% missing).
    """
    print("\n  6.1: Injecting callrate outliers...")

    set_seed(10)
    n_outliers = OUTLIER_CONFIG['callrate']
    n_variants = genotypes.shape[0]

    # Select samples not already used
    available = [i for i in range(N_SAMPLES) if i not in used_samples]
    outlier_indices = np.random.choice(available, n_outliers, replace=False)

    for idx in outlier_indices:
        # Set 15-25% of genotypes to missing
        miss_rate = np.random.uniform(0.15, 0.25)
        n_to_miss = int(n_variants * miss_rate)
        miss_var_indices = np.random.choice(n_variants, n_to_miss, replace=False)
        genotypes[miss_var_indices, idx] = -9  # Missing code

    outlier_ids = [f'SAMP_{i:04d}' for i in outlier_indices]
    outlier_manifest['callrate_samples'] = outlier_ids
    used_samples.update(outlier_indices)

    print(f"    Added {n_outliers} samples with high missingness")
    return genotypes


def inject_sex_discordance(genotypes, psam, outlier_manifest, used_samples):
    """
    Inject sex discordance by modifying X chromosome genotypes.
    Males appear female (heterozygous X) and vice versa.
    """
    print("\n  6.2: Injecting sex discordance outliers...")

    set_seed(11)
    n_outliers = OUTLIER_CONFIG['sex_discordance']

    # X variants are at the end
    x_start_idx = N_AUTOSOMAL_VARIANTS

    # Select samples not already used
    available = [i for i in range(N_SAMPLES) if i not in used_samples]
    outlier_indices = np.random.choice(available, n_outliers, replace=False)

    for idx in outlier_indices:
        reported_sex = psam.loc[idx, 'SEX']
        x_genos = genotypes[x_start_idx:, idx].copy()

        # Get non-missing indices
        non_missing = x_genos != -9

        if reported_sex == 1:  # Reported male, make appear female (heterozygous)
            x_genos[non_missing] = 1  # Heterozygous
        else:  # Reported female, make appear male (homozygous)
            x_genos[non_missing] = np.random.choice([0, 2], size=np.sum(non_missing))

        genotypes[x_start_idx:, idx] = x_genos

    outlier_ids = [f'SAMP_{i:04d}' for i in outlier_indices]
    outlier_manifest['sex_discordance_samples'] = outlier_ids
    used_samples.update(outlier_indices)

    print(f"    Added {n_outliers} samples with sex discordance")
    return genotypes


def inject_het_outliers(genotypes, outlier_manifest, used_samples):
    """
    Inject heterozygosity outliers by modifying genotype patterns.
    """
    print("\n  6.3: Injecting heterozygosity outliers...")

    set_seed(12)
    n_high = OUTLIER_CONFIG['het_high']
    n_low = OUTLIER_CONFIG['het_low']

    # Only modify autosomal variants
    auto_end = N_AUTOSOMAL_VARIANTS

    # Select samples not already used
    available = [i for i in range(N_SAMPLES) if i not in used_samples]

    # High heterozygosity (excess het, low F)
    high_het_indices = np.random.choice(available, n_high, replace=False)
    available = [i for i in available if i not in high_het_indices]

    # Low heterozygosity (deficient het, high F)
    low_het_indices = np.random.choice(available, n_low, replace=False)

    # Inject high het (convert hom to het)
    for idx in high_het_indices:
        auto_genos = genotypes[:auto_end, idx].copy()
        hom_mask = (auto_genos == 0) | (auto_genos == 2)
        hom_mask &= (auto_genos != -9)  # Exclude missing
        n_hom = np.sum(hom_mask)
        n_to_convert = int(n_hom * 0.8)
        convert_indices = np.random.choice(np.where(hom_mask)[0], n_to_convert, replace=False)
        auto_genos[convert_indices] = 1
        genotypes[:auto_end, idx] = auto_genos

    # Inject low het (convert het to hom)
    for idx in low_het_indices:
        auto_genos = genotypes[:auto_end, idx].copy()
        het_mask = (auto_genos == 1)
        n_het = np.sum(het_mask)
        n_to_convert = int(n_het * 0.8)
        convert_indices = np.random.choice(np.where(het_mask)[0], n_to_convert, replace=False)
        auto_genos[convert_indices] = np.random.choice([0, 2], n_to_convert)
        genotypes[:auto_end, idx] = auto_genos

    outlier_ids_high = [f'SAMP_{i:04d}' for i in high_het_indices]
    outlier_ids_low = [f'SAMP_{i:04d}' for i in low_het_indices]
    outlier_manifest['het_high_samples'] = outlier_ids_high
    outlier_manifest['het_low_samples'] = outlier_ids_low
    used_samples.update(high_het_indices)
    used_samples.update(low_het_indices)

    print(f"    Added {n_high} high-het and {n_low} low-het samples")
    return genotypes


def inject_related_pairs(genotypes, outlier_manifest, used_samples):
    """
    Create related sample pairs by copying genotypes with controlled differences.
    """
    print("\n  6.4: Injecting related sample pairs...")

    set_seed(13)
    n_duplicates = OUTLIER_CONFIG['duplicate_pairs']
    n_related = OUTLIER_CONFIG['related_pairs']
    n_variants = genotypes.shape[0]

    # Select samples for pairs
    available = [i for i in range(N_SAMPLES) if i not in used_samples]

    # Need 2 samples per pair
    needed = (n_duplicates + n_related) * 2
    if len(available) < needed:
        raise ValueError(f"Not enough available samples for related pairs. Need {needed}, have {len(available)}")

    selected = np.random.choice(available, needed, replace=False)

    duplicate_pairs = []
    related_pairs = []

    # Create duplicate pairs (>99% shared)
    for i in range(n_duplicates):
        src_idx = selected[i * 2]
        tgt_idx = selected[i * 2 + 1]

        # Copy genotypes with ~1% difference
        genotypes[:, tgt_idx] = genotypes[:, src_idx].copy()
        n_diff = int(n_variants * 0.01)
        diff_indices = np.random.choice(n_variants, n_diff, replace=False)
        genotypes[diff_indices, tgt_idx] = np.random.choice([0, 1, 2], n_diff)

        duplicate_pairs.append((f'SAMP_{src_idx:04d}', f'SAMP_{tgt_idx:04d}'))
        used_samples.update([src_idx, tgt_idx])

    # Create related pairs (~50% shared for first-degree)
    offset = n_duplicates * 2
    for i in range(n_related):
        src_idx = selected[offset + i * 2]
        tgt_idx = selected[offset + i * 2 + 1]

        # Share ~50% of genotypes (first-degree relative simulation)
        genotypes[:, tgt_idx] = genotypes[:, src_idx].copy()
        n_diff = int(n_variants * 0.50)  # 50% different
        diff_indices = np.random.choice(n_variants, n_diff, replace=False)
        genotypes[diff_indices, tgt_idx] = np.random.choice([0, 1, 2], n_diff)

        related_pairs.append((f'SAMP_{src_idx:04d}', f'SAMP_{tgt_idx:04d}'))
        used_samples.update([src_idx, tgt_idx])

    outlier_manifest['duplicate_pairs'] = duplicate_pairs
    outlier_manifest['related_pairs'] = related_pairs

    print(f"    Created {n_duplicates} duplicate pairs and {n_related} related pairs")
    return genotypes


def inject_variant_missingness(genotypes, outlier_manifest, used_variants):
    """
    Add high missingness to selected variants (>5% missing).
    """
    print("\n  6.5: Injecting variant missingness outliers...")

    set_seed(14)
    n_outliers = OUTLIER_CONFIG['geno_miss_variants']

    # Only modify autosomal variants
    available = [i for i in range(N_AUTOSOMAL_VARIANTS) if i not in used_variants]
    outlier_var_indices = np.random.choice(available, n_outliers, replace=False)

    for var_idx in outlier_var_indices:
        # Set 10-15% of samples to missing
        miss_rate = np.random.uniform(0.10, 0.15)
        n_to_miss = int(N_SAMPLES * miss_rate)
        miss_sample_indices = np.random.choice(N_SAMPLES, n_to_miss, replace=False)
        genotypes[var_idx, miss_sample_indices] = -9

    used_variants.update(outlier_var_indices)
    outlier_manifest['geno_miss_variants'] = list(outlier_var_indices)

    print(f"    Added {n_outliers} variants with high missingness")
    return genotypes


def inject_case_control_missingness(genotypes, psam, outlier_manifest, used_variants):
    """
    Add differential missingness between cases and controls.
    """
    print("\n  6.6: Injecting case-control differential missingness...")

    set_seed(15)
    n_outliers = OUTLIER_CONFIG['case_control_variants']

    case_indices = psam[psam['PHENO1'] == 2].index.tolist()
    control_indices = psam[psam['PHENO1'] == 1].index.tolist()

    # Select variants not already used
    available = [i for i in range(N_AUTOSOMAL_VARIANTS) if i not in used_variants]
    outlier_var_indices = np.random.choice(available, n_outliers, replace=False)

    for var_idx in outlier_var_indices:
        # High missingness in cases (20%), low in controls (2%)
        n_case_miss = int(len(case_indices) * 0.20)
        n_ctrl_miss = int(len(control_indices) * 0.02)

        case_miss_samples = np.random.choice(case_indices, n_case_miss, replace=False)
        ctrl_miss_samples = np.random.choice(control_indices, n_ctrl_miss, replace=False)

        genotypes[var_idx, case_miss_samples] = -9
        genotypes[var_idx, ctrl_miss_samples] = -9

    used_variants.update(outlier_var_indices)
    outlier_manifest['case_control_variants'] = list(outlier_var_indices)

    print(f"    Added {n_outliers} variants with differential missingness")
    return genotypes


def inject_hwe_violations(genotypes, outlier_manifest, used_variants):
    """
    Create HWE violations by converting heterozygotes to homozygotes.
    """
    print("\n  6.7: Injecting HWE violation variants...")

    set_seed(16)
    n_outliers = OUTLIER_CONFIG['hwe_violation_variants']

    # Select variants not already used
    available = [i for i in range(N_AUTOSOMAL_VARIANTS) if i not in used_variants]
    outlier_var_indices = np.random.choice(available, n_outliers, replace=False)

    for var_idx in outlier_var_indices:
        var_genos = genotypes[var_idx].copy()
        # Convert all heterozygotes to homozygotes (creates excess homozygosity)
        het_mask = var_genos == 1
        var_genos[het_mask] = np.random.choice([0, 2], np.sum(het_mask))
        genotypes[var_idx] = var_genos

    used_variants.update(outlier_var_indices)
    outlier_manifest['hwe_violation_variants'] = list(outlier_var_indices)

    print(f"    Added {n_outliers} variants with HWE violations")
    return genotypes


def inject_all_outliers(genotypes, psam, outlier_manifest):
    """
    Orchestrate injection of all outlier types.
    """
    print("\nStep 6: Injecting outliers...")

    used_samples = set()
    used_variants = set()

    # Sample-level outliers
    genotypes = inject_callrate_outliers(genotypes, outlier_manifest, used_samples)
    genotypes = inject_sex_discordance(genotypes, psam, outlier_manifest, used_samples)
    genotypes = inject_het_outliers(genotypes, outlier_manifest, used_samples)
    genotypes = inject_related_pairs(genotypes, outlier_manifest, used_samples)

    # Variant-level outliers
    genotypes = inject_variant_missingness(genotypes, outlier_manifest, used_variants)
    genotypes = inject_case_control_missingness(genotypes, psam, outlier_manifest, used_variants)
    genotypes = inject_hwe_violations(genotypes, outlier_manifest, used_variants)

    return genotypes


# =============================================================================
# Step 7: Finalize Dataset
# =============================================================================

def save_outlier_manifest(outlier_manifest, output_dir):
    """Save the outlier manifest to CSV."""
    print("\nStep 7: Saving outlier manifest...")

    manifest_path = output_dir / "outlier_manifest.csv"

    rows = []
    for outlier_type, values in outlier_manifest.items():
        if isinstance(values, list) and len(values) > 0:
            if isinstance(values[0], tuple):
                # Handle pairs
                for pair in values:
                    rows.append({
                        'outlier_type': outlier_type,
                        'sample_1': pair[0],
                        'sample_2': pair[1],
                        'variant_idx': None
                    })
            elif isinstance(values[0], str):
                # Handle sample IDs
                for sample in values:
                    rows.append({
                        'outlier_type': outlier_type,
                        'sample_1': sample,
                        'sample_2': None,
                        'variant_idx': None
                    })
            elif isinstance(values[0], (int, np.integer)):
                # Handle variant indices
                for var_idx in values:
                    rows.append({
                        'outlier_type': outlier_type,
                        'sample_1': None,
                        'sample_2': None,
                        'variant_idx': int(var_idx)
                    })

    manifest_df = pd.DataFrame(rows)
    manifest_df.to_csv(manifest_path, index=False)

    print(f"  Saved manifest with {len(rows)} outlier records to {manifest_path}")
    return manifest_path


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Generate synthetic test data for GenoTools")
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path(__file__).parent,
        help='Output directory for generated files'
    )
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("GenoTools Synthetic Test Data Generator")
    print("=" * 60)
    print(f"\nConfiguration:")
    print(f"  Samples: {N_SAMPLES}")
    print(f"  Autosomal variants: {N_AUTOSOMAL_VARIANTS}")
    print(f"  X chromosome variants: {N_X_VARIANTS}")
    print(f"  Random seed: {RANDOM_SEED}")
    print(f"  Output directory: {output_dir}")

    # Initialize outlier manifest
    outlier_manifest = {}

    # Step 1: Extract autosomal variants from reference
    autosomal_variants = extract_autosomal_variants(output_dir)

    # Step 2: Generate X chromosome variants
    x_variants = generate_x_variants(output_dir)

    # Step 3: Generate base genotypes
    raw_prefix = generate_base_genotypes(output_dir)

    # Step 4: Create combined variant file
    pvar_path, all_variants = create_variant_file(autosomal_variants, x_variants, output_dir)

    # Step 5: Create sample metadata
    psam_path, psam = create_sample_metadata(output_dir)

    # Load genotypes from raw pgen
    print("\n  Loading genotypes for modification...")
    raw_pgen = Path(str(raw_prefix) + ".pgen")
    genotypes = load_genotypes(raw_pgen, N_SAMPLES, N_TOTAL_VARIANTS)
    print(f"    Loaded genotype matrix: {genotypes.shape}")

    # Step 6: Inject outliers
    genotypes = inject_all_outliers(genotypes, psam, outlier_manifest)

    # Step 7: Save modified genotypes using pgenlib directly
    print("\n  Saving modified genotypes...")
    final_prefix = output_dir / "genotools_test"
    save_genotypes_pgenlib(genotypes, final_prefix, pvar_path, psam_path)

    # Step 8: Save outlier manifest
    save_outlier_manifest(outlier_manifest, output_dir)

    # Cleanup temporary files
    print("\n  Cleaning up temporary files...")
    for pattern in ['*_raw.*', 'autosomal_variants.txt', '*.log']:
        for f in output_dir.glob(pattern):
            if f.is_file():
                f.unlink()

    print("\n" + "=" * 60)
    print("DONE!")
    print("=" * 60)
    print(f"\nGenerated files:")
    print(f"  {final_prefix}.pgen")
    print(f"  {final_prefix}.pvar")
    print(f"  {final_prefix}.psam")
    print(f"  {output_dir / 'outlier_manifest.csv'}")

    print(f"\nOutlier summary:")
    for key, value in outlier_manifest.items():
        if isinstance(value, list):
            print(f"  {key}: {len(value)}")

    return final_prefix


if __name__ == '__main__':
    main()
