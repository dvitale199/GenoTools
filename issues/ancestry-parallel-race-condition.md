# Race condition: concurrent ancestry runs corrupt shared reference panel intermediates

## Status
Open — tracked for future fix

## Summary

Running two ancestry predictions concurrently against the same `--ref-panel` causes race conditions. The legacy `get_raw_files()` method writes intermediate files derived from the reference panel path, so parallel invocations overwrite each other's intermediates.

This is why the A/B test script (`tests/scripts/test_ancestry_ab.py`) and the previous `test_old_vs_new.py` run ancestry jobs **sequentially** rather than in parallel.

## Affected code

`genotools/ancestry.py` — `Ancestry.get_raw_files()` creates intermediate files using paths like:

```
{out_path}_variant_pruned.*
{out_path}_common_snps.*
{out_path}_umap_linearsvc_ancestry_model.*
```

While the `out_path`-prefixed files are unique per run, the method also operates on the reference panel bfiles directly:

- `plink2 --bfile {ref_panel} --extract ... --make-bed --out {ref_common_snps}` — writes new bfiles from ref panel
- Allele harmonization reads/writes files adjacent to the ref panel path
- `get_common_snps()` in `utils.py` creates merge IDs and writes `.common_snps` files

When two processes run simultaneously with the same `--ref-panel`, these PLINK operations can collide on the same output paths.

## Reproduction

```bash
# These two commands running simultaneously will corrupt each other:
genotools --pfile data1 --out results/run1 --ref-panel shared_ref --ref-labels labels --ancestry &
genotools --pfile data2 --out results/run2 --ref-panel shared_ref --ref-labels labels --ancestry &
wait
```

## Current workaround

Run ancestry jobs sequentially. This is enforced in test scripts but not in user-facing documentation.

## Proposed fix

Isolate all reference panel intermediate files into per-run temporary directories:

1. At the start of `get_raw_files()`, copy (or symlink) the reference panel bfiles into the run's tmp directory
2. Derive all intermediate paths from `{tmp_dir}/{ref_panel_name}` instead of the original ref panel path
3. This makes each run fully independent — no shared writable state

This would also enable parallel ancestry runs in production (e.g., processing multiple cohorts against the same reference panel).

### Implementation notes

- The copy/symlink should happen in `_run_ancestry_prediction()` (runner.py) before calling `get_raw_files()`
- Symlinks are preferred over copies for the ~2GB reference panel bfiles
- The tmp directory is already created per run (`self.state.tmp_dir`)
- After the fix, remove the sequential constraint from test scripts
