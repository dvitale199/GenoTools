# Test Data

This directory contains test data for GenoTools regression testing.

## Directory Structure

```
data/
├── synthetic/                  # Synthetic test data with known outliers
│   ├── generate_test_data.py  # Script to generate synthetic data
│   ├── genotools_test.pgen    # Genotype data (binary)
│   ├── genotools_test.pvar    # Variant info
│   ├── genotools_test.psam    # Sample info
│   └── outlier_manifest.csv   # Records of injected outliers
└── README.md                   # This file
```

## Synthetic Test Data

The synthetic dataset contains **500 samples** and **40,500 variants** (40,000 autosomal + 500 X chromosome).

### Injected Outliers

The following outliers are intentionally injected for testing QC step detection:

| Outlier Type | Count | Description |
|-------------|-------|-------------|
| callrate_samples | 12 | Samples with >15% missingness |
| sex_discordance_samples | 12 | Samples with genetic sex ≠ reported sex |
| het_high_samples | 4 | Samples with excess heterozygosity |
| het_low_samples | 4 | Samples with deficient heterozygosity |
| duplicate_pairs | 4 | Sample pairs that are nearly identical |
| related_pairs | 8 | Sample pairs with ~50% genotype sharing |
| geno_miss_variants | 150 | Variants with >10% missingness |
| case_control_variants | 75 | Variants with differential missingness |
| hwe_violation_variants | 75 | Variants violating Hardy-Weinberg equilibrium |

### Regenerating Test Data

To regenerate the synthetic test data:

```bash
# Activate the development environment
source .venv/bin/activate

# Run the generator
python tests/data/synthetic/generate_test_data.py
```

**Requirements:**
- Reference panel at `~/.genotools/ref/ref_panel/ref_panel_gp2_prune_rm_underperform_pos_update.bim`
- PLINK2 installed and accessible
- pgenlib Python package

### Outlier Manifest Format

The `outlier_manifest.csv` file contains:

| Column | Description |
|--------|-------------|
| outlier_type | Type of outlier (matches categories above) |
| sample_1 | Sample ID (for sample-level outliers) |
| sample_2 | Second sample ID (for pairs) |
| variant_idx | Variant index (for variant-level outliers) |

## Using Test Data

### In Pytest

```python
@pytest.fixture
def test_geno_path():
    return Path("tests/data/synthetic/genotools_test")
```

### Direct Usage

```bash
# Run GenoTools on test data
genotools --pfile tests/data/synthetic/genotools_test --out /tmp/test_output --callrate 0.02
```

## Notes

- The synthetic data uses variant IDs from the GP2 reference panel for realistic overlap testing
- X chromosome variants are synthetic with positions in the standard X range
- Random seed is fixed (42) for reproducibility
