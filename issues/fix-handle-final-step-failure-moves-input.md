# Bug: `_handle_final_step_failure` moves original input files

## Summary

When running without `--full-output` (default) and with `--warn` (default), if all QC steps fail, `_handle_final_step_failure` calls `os.rename()` on the user's original input pfiles, moving them to the output location. The original files are destroyed.

## Reproduction

```bash
genotools --pfile /my/data --out /tmp/result --all-sample
```

After this command, `/my/data.pgen`, `/my/data.psam`, `/my/data.pvar` no longer exist — they've been moved to `/tmp/result.pgen`, etc.

## Root cause

`runner.py`, `_handle_final_step_failure` (lines ~913-922):

```python
move_path = (
    geno_path  # ← the user's ORIGINAL input path
    if self.args.full_output or not self.args.ancestry.run_ancestry
    else working_geno
)
if os.path.isfile(f"{move_path}.pgen"):
    os.rename(f"{move_path}.pgen", f"{out_path}.pgen")
    os.rename(f"{move_path}.psam", f"{out_path}.psam")
    os.rename(f"{move_path}.pvar", f"{out_path}.pvar")
```

The failure chain:

1. `full_output=False` (default) — pipeline sets `working_geno` to a tmp dir path
2. Input files are NOT copied to the tmp dir
3. `_compute_step_paths` uses `working_geno` as the first step's input
4. `warn_only=True` (default) — steps are skipped (not errored) when `.pgen` doesn't exist
5. All steps fail → `_handle_final_step_failure` fires
6. `move_path` resolves to the original `geno_path` (not `working_geno`) because `full_output=False` and `ancestry=False`
7. `os.rename` moves the original input files

## Conditions

All three must be true (all are defaults):
- `--full-output` is **not** passed
- `--warn` mode is active (default `warn_only=True`)
- Ancestry is **not** run

## Related: `_cleanup_intermediate_files`

This method also calls `os.remove()` on intermediate pfiles, but it has a guard:
```python
if remove_path != geno_path:
```
So it does NOT affect original input. Only `_handle_final_step_failure` is dangerous.

## Fix options

1. **Use `shutil.copy2` instead of `os.rename`** when the source is outside the tmp dir
2. **Never rename the original input** — add a guard like `_cleanup_intermediate_files` has
3. **Copy input to tmp dir** at the start of `_run_qc_pipeline` when `full_output=False`, so all operations work on copies

Option 2 is the simplest. Option 3 fixes the deeper issue (steps failing because input isn't in the tmp dir).
