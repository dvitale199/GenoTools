# GenoTools AI Development Guide

**AI-Assisted Development Guide for GenoTools**

This guide provides patterns, conventions, and requirements for AI coding assistants working on the GenoTools codebase.

---

## Quick Reference

### Before You Start
1. Read this entire document
2. Understand the pipeline architecture (CLI → QC classes → PLINK execution)
3. Know the standard return dictionary format
4. Verify external tool availability (PLINK, PLINK2, KING)

### Critical Rules
- **Always validate inputs** at the start of every method
- **Always return the standard output dictionary** format
- **Never break the pipeline chain** - methods must produce valid PLINK2 pfiles
- **Log everything** via `concat_logs()`
- **Test with real PLINK files** - the pipeline depends on external executables

---

## Project Overview

GenoTools is a command-line tool for genotype quality control (QC) and ancestry prediction in genetic studies. It wraps PLINK/PLINK2 commands in a Python pipeline with ML-based ancestry inference.

### Key Components

| Component | Purpose |
|-----------|---------|
| `genotools/pipeline.py` | CLI argument parsing, pipeline orchestration |
| `genotools/qc.py` | `SampleQC` and `VariantQC` classes |
| `genotools/ancestry.py` | `Ancestry` class for ML predictions |
| `genotools/gwas.py` | `Assoc` class for GWAS/PCA |
| `genotools/utils.py` | Shell execution, file conversion helpers |
| `genotools/dependencies.py` | External tool management |

### Entry Points
```python
# Main CLI
genotools → genotools.__main__:handle_main

# Reference download
genotools-download → genotools.download_refs:handle_download
```

---

## Architecture

### Pipeline Flow
```
Input (bfile/pfile/vcf)
    ↓
Format Conversion → PLINK2 pfiles
    ↓
Ancestry Prediction (optional)
    ↓
Split by Ancestry
    ↓
┌─────────────────────────────────────┐
│ QC Pipeline (per ancestry group)    │
│   callrate → sex → het → related →  │
│   case_control → haplotype → hwe →  │
│   geno → ld → assoc                 │
└─────────────────────────────────────┘
    ↓
JSON Output + Cleaned Files
```

### Class Responsibilities

**SampleQC** - Sample-level quality control
- `run_callrate_prune()` - Remove samples with low call rates
- `run_sex_prune()` - Remove samples with sex discrepancies
- `run_het_prune()` - Remove samples with extreme heterozygosity
- `run_related_prune()` - Handle related/duplicate samples
- `run_confirming_kinship()` - Verify family relationships

**VariantQC** - Variant-level quality control
- `run_geno_prune()` - Remove variants with high missingness
- `run_case_control_prune()` - Remove variants with case/control differences
- `run_haplotype_prune()` - Remove haplotype-inconsistent variants
- `run_hwe_prune()` - Remove Hardy-Weinberg violating variants
- `run_ld_prune()` - Prune variants in linkage disequilibrium

**Ancestry** - Ancestry prediction
- PCA calculation and projection
- UMAP + XGBoost classification
- Container/cloud inference support
- Admixture handling

**Assoc** - Association analysis
- PCA preparation and execution
- GWAS execution
- Lambda/inflation calculation

---

## Code Conventions

### Standard Return Dictionary

**Every QC method must return this structure:**

```python
{
    'pass': bool,           # True if step completed successfully
    'step': str,            # Step identifier (e.g., 'callrate_prune')
    'metrics': {
        'outlier_count': int,  # Number of samples/variants pruned
        # ... other step-specific metrics
    },
    'output': {
        'pruned_samples': str,  # Path to pruned sample IDs (or None)
        'plink_out': str,       # Path to output pfiles (without extension)
        # ... other output files
    }
}
```

**Example:**
```python
out_dict = {
    'pass': process_complete,
    'step': step,
    'metrics': metrics_dict,
    'output': outfiles_dict
}
return out_dict
```

### Input Validation Pattern

**Every method must validate inputs before processing:**

```python
def run_some_prune(self, threshold=0.05):
    geno_path = self.geno_path
    out_path = self.out_path

    # 1. Check paths are set
    if geno_path is None or out_path is None:
        raise ValueError("Both geno_path and out_path must be set before calling this method.")

    # 2. Check input files exist
    if not os.path.exists(f'{geno_path}.pgen'):
        raise FileNotFoundError(f"{geno_path} does not exist.")

    # 3. Check parameter types
    if not isinstance(threshold, (int, float)):
        raise TypeError("threshold should be of type int or float.")

    # 4. Check parameter bounds
    if threshold < 0 or threshold > 1:
        raise ValueError("threshold should be between 0 and 1.")

    # ... proceed with implementation
```

### Shell Command Execution

**Use `shell_do()` for all external commands:**

```python
from genotools.utils import shell_do, concat_logs

# Execute PLINK command
plink_cmd = f"{plink2_exec} --pfile {geno_path} --mind {mind} --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}"
shell_do(plink_cmd)

# Always log after execution
listOfFiles = [f'{out_path}.log']
concat_logs(step, out_path, listOfFiles)
```

### File Path Conventions

```python
# Input/output paths never include extensions
geno_path = '/path/to/data'      # Files: data.pgen, data.pvar, data.psam
out_path = '/path/to/output'     # Will create: output.pgen, output.pvar, output.psam

# Intermediate files use step suffix
step_output = f'{out_path}_{step}'  # e.g., output_callrate_prune

# Outlier files
outliers_out = f'{out_path}.outliers'
```

### PLINK2 psam Column Standard

**Always preserve sample metadata columns:**

```python
# Use this flag for all --make-pgen commands
--make-pgen psam-cols=fid,parents,sex,pheno1,phenos
```

---

## Testing

### Current State
GenoTools does not have a formal test suite. Test data exists in `/data/` but no pytest infrastructure.

### When Adding Tests

**Recommended structure:**
```
tests/
├── conftest.py           # Fixtures for test data paths
├── test_sample_qc.py     # SampleQC method tests
├── test_variant_qc.py    # VariantQC method tests
├── test_ancestry.py      # Ancestry prediction tests
├── test_pipeline.py      # Integration tests
└── test_utils.py         # Utility function tests
```

**Test pattern:**
```python
import pytest
from genotools.qc import SampleQC

class TestSampleQC:
    """Tests for SampleQC class."""

    def test_callrate_prune_valid_input(self, sample_pfiles):
        """Test callrate pruning with valid input."""
        # Arrange
        qc = SampleQC(geno_path=sample_pfiles, out_path='/tmp/test_out')

        # Act
        result = qc.run_callrate_prune(mind=0.02)

        # Assert
        assert result['pass'] is True
        assert result['step'] == 'callrate_prune'
        assert 'outlier_count' in result['metrics']

    def test_callrate_prune_invalid_mind(self, sample_pfiles):
        """Test that invalid mind value raises ValueError."""
        qc = SampleQC(geno_path=sample_pfiles, out_path='/tmp/test_out')

        with pytest.raises(ValueError, match="mind should be between 0 and 1"):
            qc.run_callrate_prune(mind=1.5)
```

---

## Common Pitfalls

### 1. Missing await/async
GenoTools is **synchronous**. Do not add async/await unless refactoring the entire pipeline.

### 2. Breaking the pfile chain
Every QC step must:
- Read from `{geno_path}.pgen/.pvar/.psam`
- Write to `{out_path}.pgen/.pvar/.psam`

```python
# ❌ BAD - Breaks the chain
result = run_analysis(geno_path)
return result  # No pfiles created

# ✅ GOOD - Maintains the chain
plink_cmd = f"{plink2_exec} --pfile {geno_path} ... --make-pgen ... --out {out_path}"
shell_do(plink_cmd)
# Now out_path.pgen exists for next step
```

### 3. Forgetting to log
```python
# ❌ BAD - No logging
shell_do(plink_cmd)
return out_dict

# ✅ GOOD - Always log
shell_do(plink_cmd)
listOfFiles = [f'{out_path}.log']
concat_logs(step, out_path, listOfFiles)
return out_dict
```

### 4. Incorrect outlier file format
```python
# Outlier files must be tab-separated with #FID header
# ❌ BAD
df.to_csv(outliers_out, sep=',')

# ✅ GOOD
df = df.rename({'FID': '#FID'}, axis=1)
df.to_csv(outliers_out, sep='\t', header=True, index=False)
```

### 5. Platform-specific code without checks
```python
# ❌ BAD - KING only works on Linux
result = run_king_analysis()

# ✅ GOOD - Check platform first
import platform
if platform.system() != 'Linux':
    print('This analysis can only run on Linux!')
    return None
result = run_king_analysis()
```

---

## Adding New QC Steps

### Template for Sample-Level QC

```python
def run_new_sample_prune(self, param1=default1, param2=default2):
    """
    Execute new sample pruning on genotype data.

    Parameters:
    - param1 (type): Description. Default is default1.
    - param2 (type): Description. Default is default2.

    Returns:
    - dict: Standard output dictionary with 'pass', 'step', 'metrics', 'output'.
    """
    geno_path = self.geno_path
    out_path = self.out_path

    # Input validation
    if geno_path is None or out_path is None:
        raise ValueError("Both geno_path and out_path must be set.")
    if not os.path.exists(f'{geno_path}.pgen'):
        raise FileNotFoundError(f"{geno_path} does not exist.")
    # ... validate param1, param2

    step = "new_sample_prune"
    outliers_out = f'{out_path}.outliers'

    # Execute PLINK commands
    plink_cmd = f"{plink2_exec} --pfile {geno_path} ... --make-pgen psam-cols=fid,parents,sex,pheno1,phenos --out {out_path}"
    shell_do(plink_cmd)

    # Log
    listOfFiles = [f'{out_path}.log']
    concat_logs(step, out_path, listOfFiles)

    # Process results and count outliers
    if os.path.isfile(f'{out_path}.some_output'):
        outliers = pd.read_csv(f'{out_path}.some_output', sep='\s+')
        outliers = outliers.rename({'FID': '#FID'}, axis=1)
        outliers.to_csv(outliers_out, sep='\t', header=True, index=False)
        outlier_count = outliers.shape[0]
    else:
        outlier_count = 0

    process_complete = True

    # Build return dictionary
    outfiles_dict = {
        'pruned_samples': outliers_out,
        'plink_out': out_path,
    }
    metrics_dict = {
        'outlier_count': outlier_count
    }
    out_dict = {
        'pass': process_complete,
        'step': step,
        'metrics': metrics_dict,
        'output': outfiles_dict
    }

    return out_dict
```

### Registering New Steps in Pipeline

After adding a new method, update `pipeline.py`:

```python
# In handle_main() or execute_pipeline()

# 1. Add to appropriate step list
samp_steps = ['callrate', 'sex', 'het', 'related', 'kinship_check', 'new_step']

# 2. Add to steps_dict mapping
steps_dict = {
    # ... existing steps
    'new_step': samp_qc.run_new_sample_prune,
}

# 3. Add argument in gt_argparse()
parser.add_argument('--new_step', type=float, nargs='?', default=None, const=0.05,
                    help='Description of new step')
```

---

## External Dependencies

### Required Executables
- **PLINK 1.9** - Legacy format support, sex checks
- **PLINK2** - Primary tool for QC operations
- **KING** - Relatedness analysis (Linux only)

### Checking Dependencies
```python
from genotools.dependencies import check_plink, check_plink2, check_king

plink_exec = check_plink()    # Returns path or downloads
plink2_exec = check_plink2()  # Returns path or downloads
king_exec = check_king()      # Returns path (Linux) or None
```

### Python Dependencies
```
pandas, numpy              # Data manipulation
scikit-learn, xgboost     # ML models
umap-learn==0.5.3         # Dimensionality reduction (pinned version)
scipy, statsmodels        # Statistics
matplotlib, seaborn       # Visualization
google-cloud-aiplatform   # Cloud predictions
```

---

## File Formats

### Input Formats
- **bfile**: PLINK 1.9 (`.bed`, `.bim`, `.fam`)
- **pfile**: PLINK 2 (`.pgen`, `.pvar`, `.psam`)
- **VCF**: Variant Call Format (`.vcf`)

### Internal Format
All processing uses PLINK2 pfiles. Conversion happens automatically at pipeline start.

### Required psam Columns
```
#FID    IID    SEX    PHENO1    [additional columns]
```
- `SEX`: 0=unknown, 1=male, 2=female
- `PHENO1`: -9=missing, 1=control, 2=case

### Output JSON Structure
```json
{
    "ancestry_counts": {"EUR": 100, "AFR": 50, ...},
    "ancestry_labels": [{"#FID": "...", "IID": "...", "label": "EUR"}, ...],
    "QC": [
        {"step": "callrate_prune", "pruned_count": 5, "metric": "outlier_count", ...}
    ],
    "GWAS": [
        {"value": 1.02, "metric": "lambda", "ancestry": "EUR"}
    ],
    "pruned_samples": [...],
    "related_samples": [...]
}
```

---

## Ancestry Prediction

### Model Architecture
- **Step 1**: PCA on merged reference + input samples
- **Step 2**: UMAP transformation
- **Step 3**: XGBoost classification
- **Supported labels**: AFR, SAS, EAS, EUR, AMR, AJ, CAS, MDE, FIN, AAC

### Container Support
```bash
# Docker
genotools --pfile input --out output --ancestry --container

# Singularity
genotools --pfile input --out output --ancestry --singularity

# Google Cloud
genotools --pfile input --out output --ancestry --cloud
```

---

## Development Workflow

### Virtual Environment Setup

Use separate virtual environments to keep stable and development versions isolated:

| Environment | Install Command | Purpose |
|-------------|-----------------|---------|
| `.venv` | `pip install -e .` | Active development, code changes reflected immediately |
| `.venv-stable` | `pip install .` | Frozen baseline snapshot for regression comparison |

```bash
# Development environment (editable - use this for active work)
python -m venv .venv
source .venv/bin/activate
pip install -e .

# Stable baseline (frozen snapshot - use for comparison)
python -m venv .venv-stable
source .venv-stable/bin/activate
pip install .
```

**Note:** Both are installed from the local repo, not PyPI, since PyPI may be outdated.

### 1. Making Changes
```bash
# Install in development mode
pip install -e .

# Run linting (recommended)
flake8 genotools/
black genotools/

# Test manually
genotools --pfile test_data --out test_output --callrate
```

### 2. Testing Changes
```bash
# Basic pipeline test
genotools --pfile data/test --out output/test --all_sample --all_variant

# With ancestry
genotools --pfile data/test --out output/test --ancestry --ref_panel /path/to/ref --ref_labels /path/to/labels
```

### 3. Committing
```bash
# Version bump in setup.py if needed
# Update CITATION.cff if adding authors
git add .
git commit -m "Description of changes"
```

---

## Error Handling Checklist

When AI generates code, verify:

- [ ] All inputs validated at method start
- [ ] FileNotFoundError raised for missing files
- [ ] ValueError raised for invalid parameters
- [ ] TypeError raised for wrong parameter types
- [ ] Platform checks for OS-specific features
- [ ] Proper error messages guide users to logs
- [ ] `--warn` flag behavior respected (continue on failure)

---

## Quick Debugging

### Check PLINK availability
```python
from genotools.dependencies import check_plink, check_plink2
print(check_plink())   # Should print path
print(check_plink2())  # Should print path
```

### Inspect pipeline state
```python
# After pipeline runs, check pass_fail dict
for step, status in out_dict['pass_fail'].items():
    print(f"{step}: {'PASS' if status['status'] else 'FAIL'}")
```

### Read logs
```bash
cat output_all_logs.log      # Full PLINK output
cat output_cleaned_logs.log  # Formatted summary
```

---

## Summary

| Rule | Reason |
|------|--------|
| Validate all inputs | Catch errors early with clear messages |
| Use standard return dict | Pipeline depends on consistent structure |
| Always log with concat_logs() | Debugging and audit trail |
| Maintain pfile chain | Each step feeds the next |
| Check platform for KING | Only works on Linux |
| Pin umap-learn==0.5.3 | Model compatibility |
