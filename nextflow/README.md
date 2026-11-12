# EGA Hail / Nextflow pipeline — AC0 embedded in M4

This package updates the supplied Nextflow Hail pipeline. Most importantly,
**AC0 is strictly an optional VCF FILTER annotation; variants are never removed
by the AC0 step**. The rest of the AF pipeline works as before.

## Run

1. From this directory, edit **`config.yaml`** (the one alongside `main.nf`):

   ```yaml
   ac0_flagging: true
   ```

   Set `false` to disable. This switch is independent of the existing
   `genotype_filtering` flag. The old standalone `AF_hail_pipeline/config.yaml`
   does **not** control the Nextflow workflow.

2. Check `vcf_dir`, `ref_gen`, and the other existing QC parameters. If using
   the example input included with the ZIP, the existing defaults point at it.
   The included example only declares `GT`, so **AC0 is expected to be skipped**.

3. Rebuild the Docker image:

   ```bash
   docker build -t af_hail_pipeline:latest ./AF_hail_pipeline
   ```

4. Launch with the installed Nextflow CLI (same as the original workflow):

   ```bash
   nextflow run main.nf -c nextflow.config
   ```

## AC0 behavior

The AC0 check runs in `M4_AF_ANNOTATION`, on the MatrixTable **after** existing
preprocessing, sample filtering, relatedness and ancestry steps, and **before**
M4 allele-frequency aggregation. It does not change `GT` or recalculated
`AF_total_recalc`, `AC_total_recalc`, or `AN_total_recalc`.

A genotype is considered high-quality when it has `GQ >= 20` and
`DP >= 10` (diploid) or `DP >= 5` (haploid). For heterozygous calls, the
allele depth `AD` divided by `DP` must be at least `0.2` for each called ALT
allele (the gnomAD adjusted-genotype criterion). Phased `1|0` is handled like
`0|1`. `0/0` contributes no ALT copies even when it passes these checks.

For each VCF site, if the **sum of all high-quality alternative allele copies
is zero**, the pipeline adds `AC0` to its existing **FILTER** set:

```
#CHROM  POS   ...  FILTER        INFO   ...
1       10001 ...  AC0           ...
1       10002 ...  LowQual       ...
```

All original FILTER flags are preserved. For a site with multiple ALT alleles
and `split_multiallelic: false`, AC0 is *site-level*: it is added only when
**none** of the ALTs has a high-quality genotype. For allele-specific behavior,
set `split_multiallelic: true` (pre-existing pipeline option).

**Important:** the existing AF/AC calculations remain based on the pipeline's
normal GT values. AC0 uses an *independent* adjusted-genotype check and does
not replace or recalculate those INFO fields. A site may therefore have
`FILTER=AC0` while `INFO/AC_total_recalc` is nonzero; that is expected and
means all ALT support was low-quality by AC0 criteria.

## Missing FORMAT policy: always continue

The required **VCF FORMAT** fields are `GT`, `GQ`, `DP`, and `AD`.
When `ac0_flagging: true`:

* The pipeline inspects the FORMAT declarations in **every input VCF**. If any
  required field is missing from any file, **AC0 is skipped for the entire
  merged dataset**, and the missing fields and files are reported. The M4 AF
  processing continues unchanged.
* When the required declarations exist, missing **genotype values** are
  handled **per variant**, not per dataset. If a called GT has missing `GQ`
  or `DP`, or a heterozygous GT has missing/invalid `AD`, AC0 is **not
  evaluated for that variant**. Other complete variants are still evaluated
  and may receive AC0. If no variant is evaluable, AC0 is skipped entirely.
  Incomplete variants retain their original FILTER values.
* Missing `GT` values (`./.` / no calls) are normal and are counted in the
  report, but do not invalidate other genotype calls.
* If the optional AC0 check encounters another error, it logs the error and
  falls back to the original MatrixTable without interrupting M4.

The missing-FORMAT policy only applies to this **optional** AC0 feature.
Existing required pipeline steps can still fail for independent reasons such as
malformed VCFs, incompatible multi-VCF headers, or missing fields required by
pre-existing QC settings. This change does not disable their error handling.

### Where to find the result

The following files are **always published**, even with `save_logs=false` or
`ac0_flagging=false`:

```
results/qc/ac0/ac0_status.tsv
results/qc/ac0/ac0_status.json
```

The same AC0 status appears in the `M4_AF_ANNOTATION` process log and the
existing `WorkflowStats_*.csv` (when `save_logs=true`). Possible statuses:
`DISABLED`, `APPLIED`, `PARTIALLY_APPLIED_MISSING_FORMAT_VALUES`,
`SKIPPED_MISSING_FORMAT_FIELDS`, `SKIPPED_MATRIX_SCHEMA`,
`SKIPPED_MISSING_FORMAT_VALUES`, `SKIPPED_INPUT_AUDIT`, `SKIPPED_AC0_ERROR`.

`ac0_status.json` includes per-VCF missing FORMAT fields, missing-value
counts, numbers of variants evaluated / skipped / flagged, and the thresholds.
`ac0_status.tsv` includes the same totals for straightforward inspection.

## Test

Fast header/report unit tests (Hail not required):

```bash
python -m unittest discover -s tests -p 'test_ac0_header.py' -v
```

Optional Hail integration smoke test after rebuilding the Docker image:

```bash
docker run --rm \
    -v "$PWD/tests:/tests:ro" \
    -v "$PWD/config.yaml:/app/config.yaml:ro" \
    af_hail_pipeline:latest python /tests/test_ac0_hail.py
```

This tests `0/0` REF behavior, failed ALT calls, existing FILTER preservation,
phased heterozygotes, haploid DP, retained row count, and per-variant missing GQ handling.

Note: the tests in `tests/` are code-level checks maintained with the pipeline;
no standalone AC0 validation pack or its synthetic test-input folder is included.

## AC0 implementation location

The AC0 inspection, filtering and reporting functions are now integrated in
`AF_hail_pipeline/M4_af_annotation.py` (under its *Optional AC0 FILTER flagging*
section). `nextflow_stage.py` calls `m4.maybe_flag_ac0()` before M4 AF
annotation and uses `m4.write_ac0_reports()` and
`m4.append_to_workflow_csv()` for outputs. No separate `ac0_flagging.py`
module is needed, and AC0 is performed only once.

The Hail per-variant eligibility mask and quality-adjusted allele count use
expressions from the current MatrixTable at every annotation step. This avoids
Hail's `Cannot combine expressions from different source objects` error.

## Packaging notes

This ZIP includes the pipeline source, both original configs (the root config
is the one used by Nextflow), Dockerfile, GrafAnc binary/annotations, the
original small example input VCF, tests and documentation. No AC0 validation
pack, synthetic AC0 VCF fixtures or `ac0_test/` folder is included.
Generated `work/`, `results/`, `.nextflow*`, Python caches and the inner Git
history were excluded to avoid distributing large stale intermediate files.

For reference, the thresholds are derived from the official
[gnomAD `get_adj_expr` documentation](https://broadinstitute.github.io/gnomad_methods/api_reference/utils/annotations.html)
and VCF FILTER handling follows
[Hail export VCF documentation](https://hail.is/docs/0.2/methods/impex.html).

## Configuration-file mount fix (M1 FileNotFoundError)

If `utils.py` raised `FileNotFoundError: config.yaml` at the first Nextflow
process, this build contains a two-part fix:

- `nextflow.config` bind-mounts the **selected** `--hail_config` YAML to
  `/app/config.yaml` (read-only) in every Docker task.
- `AF_hail_pipeline/utils.py` tries the normal staged `config.yaml` first,
  then falls back to `/app/config.yaml` if its task-local path cannot be read.

This is independent of AC0 and does not change your quality thresholds or
filtering logic.
