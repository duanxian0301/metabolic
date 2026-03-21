# Nonlipid Adaptation Guide

This guide explains how to adapt the clean lipid pipeline to `nonlipid` or another trait module without carrying over the failed exploratory branches.

## Core principle

Do not start by editing late-stage `final8` scripts.

For a new module, the correct route is:

1. Start from the same 112-trait main-analysis universe.
2. Define a new candidate universe for the module of interest.
3. Create compact review and ultra-pure review tables for that module.
4. Re-run LDSC and EFA/ESEM from the module-specific manifests.
5. Only then prepare factor GWAS and final LDSC validation.

## Files that are module-specific

For a new module, these input tables must be newly created:

- candidate trait table
- compact review table
- compact kept/dropped tables
- ultra-pure review table
- ultra-pure kept/dropped tables
- final review table
- final kept/dropped tables
- final trait manifest
- final factor structure summary

The lipid versions of these files are examples only; they should not be reused directly for a nonlipid analysis.

## Minimal schema to preserve

The downstream scripts assume review tables keep stable columns such as:

- `study_accession`
- `trait_code`
- `biomarker_name`
- `group`
- `sumstats_file`

Depending on stage, additional columns are expected:

- compact stage:
  - `selection_status`
  - `selection_reason`
  - `proposed_module`
  - `marker_role`
- ultra-pure stage:
  - `selection_status`
  - `selection_reason`
  - `ultra_pure_factor`
- final stage:
  - `selection_status_final8` or analogous final keep/drop field
  - final-stage selection reason

If you rename stage-specific columns, also update the supplement/workbook scripts.

## Recommended duplication strategy

1. Copy `lipid_final8` to a new sibling folder such as `nonlipid_module`.
2. Rename the script files only after the module manifests are finalized.
3. Update all `panel_path`, `output_dir`, `root_dir`, and model text blocks at the top of each copied script.
4. Keep the WSL-native factor GWAS route.
5. If any factor has only two indicators, do not split factors into separate `commonfactorGWAS` runs; use the full multi-factor `userGWAS` model instead.

## What to avoid

- Do not upload exploratory failed models as the public pipeline.
- Do not rely on temporary `tmp_*` scripts.
- Do not reuse lipid review tables as if they were generic.
- Do not split two-indicator factors into isolated `commonfactorGWAS` models.
