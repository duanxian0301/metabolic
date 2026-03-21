# GenomicGEM Metabolic Clean Upload

This upload package contains only the successful lipid-module pipeline used to produce the `lipid_final8` results.

It is intentionally separated from the larger working directory because the original worktree also contains many failed model searches, temporary tests, debug scripts, and one-off diagnostics that should not be uploaded as the public reproducible pipeline.

## What is included

- `wsl_native_environment.md`
  - The working WSL-native GenomicSEM environment notes inherited from the original `genomicGEM` workflow.
- `lipid_final8/scripts`
  - The successful scripts needed to reproduce the lipid-only pipeline from the initial 112-trait universe to final factor GWAS and factor-level LDSC validation.
- `lipid_final8/inputs/selection`
  - The manually curated review tables required to reproduce the successful `112 -> 94 -> 20 -> 10 -> 8` route.
- `lipid_final8/inputs/final_model`
  - The final 8-trait manifest and factor structure used for downstream factor GWAS.
- `templates`
  - Guidance for adapting the same workflow to `nonlipid` and other module-specific analyses.

## What is intentionally excluded

- Failed CFA branches
- Failed compact/2-factor rescue attempts
- Temporary debug scripts such as `tmp_*`
- Early Windows-only `userGWAS` attempts that were superseded by the successful WSL-native route
- One-off smoke tests and environment patch experiments

## Reproduction scope

This package is designed to reproduce the successful `lipid_final8` analysis:

1. Start from the original `Main_Zgt4_nonproportion` trait universe.
2. Retain the `Z > 4` main-analysis universe.
3. Restrict to lipid-related candidates.
4. Use the preserved manual review tables to move from 94 candidates to 20 compact traits, then to 10 ultra-pure markers, then to the final 8 markers.
5. Run the final 3-factor ESEM.
6. Prepare WSL-native factor GWAS inputs.
7. Run final 3-factor `userGWAS` in WSL.
8. Merge outputs, export standard GWAS text files, and validate factor-level LDSC.

## Script order for lipid_final8

Run the scripts in `lipid_final8/scripts` in this order:

1. `01_step1_ldsc_main112_qc.R`
2. `02_prepare_lipid94_candidates.py`
3. `03_ldsc_compact20.R`
4. `04_efa_esem_compact20.R`
5. `05_ldsc_ultrapure10.R`
6. `06_efa_esem_ultrapure10.R`
7. `07_ldsc_final8.R`
8. `08_efa_esem_final8.R`
9. `09_prepare_wsl_factor_gwas_inputs.R`
10. `10_run_wsl_usergwas_final8.R`
11. `11_merge_usergwas_final8.R`
12. `12_export_standard_usergwas_final8.R`
13. `13_ldsc_validate_factor_outputs.R`
14. `14_compile_supplement_tables.R`
15. `15_build_supplement_workbook.py`

## Important note about manual selection

This pipeline cannot be reproduced from scripts alone because the successful route includes manual selection decisions.

Those decisions are preserved here as TSV inputs:

- `lipid_module_compact_review.tsv`
- `lipid_module_compact_kept.tsv`
- `lipid_module_ultrapure3_review.tsv`
- `lipid_module_ultrapure3_kept.tsv`
- `lipid_module_ultrapure3_final8_review.tsv`
- `lipid_module_ultrapure3_final8_kept.tsv`

These files are part of the reproducible analysis input and should be versioned together with the scripts.

## How to adapt this to nonlipid or another module

Use `lipid_final8` as the working example, but do not mechanically rename lipid files only. The correct adaptation path is:

1. Build a new candidate universe from the same `Main_Zgt4_nonproportion` trait space.
2. Create new module-specific review tables with the same schema as the lipid review tables.
3. Update the top-level `panel_path`, `output_dir`, and model text in the copied scripts.
4. For downstream factor GWAS, prefer the full multi-factor `userGWAS` route if any factor has only two indicators.

See `templates/nonlipid_adaptation_guide.md`.

## Path assumptions

The uploaded scripts still use the same absolute path conventions that were used in the successful run:

- GWAS root: `D:/metabolic/GWAS`
- Main analysis root: `D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion`
- LDSC reference: `D:/LDSC/ldsc-master/eur_w_ld_chr`
- WSL R library: `/home/shenjing/R/genomicsem_fix_lib`

If the upload target machine uses different locations, these paths must be edited consistently before rerunning.
