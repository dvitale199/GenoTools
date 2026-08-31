# GenoTools JSON Output Documentation

## Overview
This documentation provides a detailed description of the keys in JSON output for the `GenoTools` package, a comprehensive tool for genetic data analysis. This output file provides users with important output and summary statistics for processes like quality control (QC), ancestry analysis, and genome-wide association studies (GWAS). Please note that output associated with steps that are not completed in a `GenoTools` run will not be output to the JSON file.

---

### Loading JSON Components into a Pandas DataFrame in Python

```
import json
import pandas as pd

json_path = "/path/to/output.json"

f = open(json_path)
json_file = json.load(f)

input_samples = pd.DataFrame(json_file['input_samples']) # replace 'input_samples' with any JSON key!
```

---

### JSON Dictionary Keys

- **`input_samples`**   
  - *Process*: `All`  
  - *Description*: Input .psam file.

- **`ancestry_counts`**  
  - *Process*: `Ancestry`  
  - *Description*: Counts of the number of samples in each predicted ancestry group.

- **`ancestry_labels`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing IDs and predicted label for each sample.

- **`confusion_matrix`**  
  - *Process*: `Ancestry`  
  - *Description*: Reference panel test set confusion matrix.
 
- **`test_accuracy`**  
  - *Process*: `Ancestry`  
  - *Description*: Reference panel test set accuracy.
 
- **`ref_pcs`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing the reference panel prinicpal components.

- **`projected_pca`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing the input samples prinicpal components.

- **`total_umap`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing the UMAP representation for the reference panel and input samples.

- **`ref_umap`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing the UMAP representation for the reference panel.

- **`new_samples_umap`**  
  - *Process*: `Ancestry`  
  - *Description*: DataFrame containing the UMAP representation for the input samples.

- **`run_info`**  
  - *Process*: `All`  
  - *Description*: *New in 2.0.1.* Which build produced the report, and from
    what command.

    | Key | Meaning |
    |---|---|
    | `version` | GenoTools version |
    | `invocation` | The command line, quoted so it can be pasted back into a shell. Absent when the pipeline was driven from Python rather than the CLI, where there is no command line to record |

- **`parameters`**  
  - *Process*: `Quality Control (QC)`  
  - *Description*: *New in 2.0.1.* The settings each step ran with. One row per
    (step, ancestry, parameter), with these columns:

    | Column | Meaning |
    |---|---|
    | `step` | Reported step name, matching the `step` column of `QC` |
    | `parameter` | Setting name, e.g. `mind`, `auto_sd`, `het_lower`. Nested settings use dotted names (`pca.n_pcs`) |
    | `value` | Its value |
    | `ancestry` | Ancestry group, or `all` for a run without `--ancestry` |
    | `source` | `requested` or `resolved` — see below |

    Before 2.0.1 the report recorded outcomes only: every top-level key was a
    count, a label or a pass/fail, and not one was a threshold. That is
    survivable while every group is cut by the same rule. `--het sd` is not —
    it derives its bounds from each dataset, so `--het -0.2 0.2 --het-ancestry
    AMR sd` produces `outlier_count` rows cut by two different rules, and
    nothing in the file said which was which.

    **`source` separates the request from the result.** `requested` is the
    configuration the step was handed; `resolved` is what the step worked out
    at runtime. `--het sd 2` is a request that resolves to different bounds in
    every ancestry group, and neither half is recoverable from the other: the
    multiplier does not predict the bounds, and the bounds do not reveal the
    multiplier. `requested` rows are the config as constructed, including
    defaults a given mode never consults — in `sd` mode `f_lower`/`f_upper`
    are present but inert, and the `resolved` row `het_mode` says which pair
    was live.

    Only the heterozygosity step derives anything today, so it is the only one
    with `resolved` rows. Every step contributes `requested` rows.

- **`QC`**  
  - *Process*: `Quality Control (QC)`  
  - *Description*: DataFrame containing summary statistics from QC steps. One row
    per (step, ancestry, metric), with these columns:

    | Column | Meaning |
    |---|---|
    | `step` | Reported step name, e.g. `callrate_prune` |
    | `pruned_count` | Samples or variants removed |
    | `metric` | Which count `pruned_count` holds, e.g. `outlier_count` |
    | `ancestry` | Ancestry group, or `all` for a run without `--ancestry` |
    | `level` | `sample` or `variant` |
    | `pass` | `true` only when the step ran and succeeded |
    | `outcome` | *New in 2.0.* `pass`, `fail`, or `skipped` |
    | `reason` | *New in 2.0.* Why, for `fail` and `skipped`; `null` otherwise |

    `pruned_count` holds a count and nothing else, which is why the settings
    live in `parameters` rather than here: a mode name or a derived bound in
    this table would sit under a column meaning "how many were pruned".

    A step the run did not request has no row. A step that was requested but
    could not run does: `outcome` distinguishes a step that **failed** from one
    the data **ruled out** (for example het on an ancestry group below PLINK's
    50-sample floor for LD estimation). Counts are zeroed for both. Before 2.0
    these rows were either absent (skips) or present without any indication of
    why (failures), so a missing row could not be told from a step that was
    never asked for.

- **`GWAS`**  
  - *Process*: `Genome-wide Assoication Study (GWAS)`  
  - *Description*: DataFrame containing lambda and lambda1000 values from completed association analysis/analyses.

- **`pruned_samples`**  
  - *Process*: `Ancestry/Quality Control (QC)`  
  - *Description*: DataFrame containing samples pruned due to insufficient ancestry count or sample-level QC.

- **`related_samples`**  
  - *Process*: `Quality Control (QC)`  
  - *Description*: DataFrame containing relatedness information output by PLINK/KING from completed relatedness pruning.

- **`pass_fail`** / **`{ANCESTRY}_pass_fail`**  
  - *Process*: `Quality Control (QC)`  
  - *Description*: Per-step status. `pass_fail` for a run without `--ancestry`;
    one `{ANCESTRY}_pass_fail` key per group otherwise. Each entry holds:

    | Key | Meaning |
    |---|---|
    | `status` | Boolean; `true` only when the step ran and succeeded |
    | `outcome` | *New in 2.0.* `pass`, `fail`, or `skipped` |
    | `reason` | *New in 2.0.* Why, for `fail` and `skipped`; `null` otherwise |
    | `input` | Pfile prefix the step read |
    | `output` | Pfile prefix the step wrote (equal to `input` for a skip, which passes the chain through untouched) |

---
