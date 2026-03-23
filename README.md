# Metabolic Clean Upload

This repository contains the successful clean GenomicSEM workflows for both the `lipid_final8` and `nonlipid_final8` modules, together with code for rebuilding supplementary tables and figures.

It is intentionally separated from the larger working directory because the original worktree also contained failed model searches, temporary tests, debug scripts, and one-off diagnostics that should not be part of the public reproducible pipeline.

## Repository layout

- `lipid_final8/`
  - Successful lipid-module pipeline from the shared Main112 universe to final factor GWAS, LDSC validation, supplementary tables, and manuscript assets.
- `nonlipid_final8/`
  - Successful nonlipid-module pipeline from the same Main112 universe to final factor GWAS, LDSC validation, and supplementary tables.
- `templates/`
  - Guidance for adapting the same workflow to additional module-specific analyses.
- `wsl_native_environment.md`
  - Working WSL-native GenomicSEM environment notes inherited from the original `genomicGEM` workflow.

## What is included

- Successful scripts only
- Manual review tables required to reproduce the retained staged-reduction paths
- Final manifests and factor structures for both modules
- Supplement-table compilation scripts and workbook builders
- Figure-generation code

## What is intentionally excluded

- Failed CFA branches
- Failed rescue attempts and abandoned model lines
- Temporary debug scripts such as `tmp_*`
- `__pycache__`, `.pyc`, and local environment artifacts
- Early Windows-only factor-GWAS attempts superseded by the successful WSL-native route
- One-off smoke tests and intermediate work directories
- Draft manuscripts, generated figures, and supplementary-result exports

## Shared starting universe

Both modules begin from the same `Main_Zgt4_nonproportion` trait universe and retain all non-proportion metabolic traits with `h2 Z > 4`.

The shared starting inputs are versioned in both module folders so that each module can be reproduced independently:

- `inputs/selection/main112_trait_manifest.tsv`
- `inputs/selection/Main_Zgt4_nonproportion_ldsc_summary.tsv`
- `inputs/selection/trait_qc_and_selection.tsv`

## Reproduction scope

### lipid_final8

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
16. `16_build_lipid_nonlipid_result_figures.py`
17. `17_build_factor_gwas_manhattan_qq.py`

### nonlipid_final8

Run the scripts in `nonlipid_final8/scripts` in this order:

1. `01_step1_ldsc_main112_qc.R`
2. `02_prepare_nonlipid18_candidates.py`
3. `03_prepare_nonlipid15_compact_review.py`
4. `04_ldsc_nonlipid_compact15.R`
5. `05_efa_esem_nonlipid_compact15.R`
6. `06_search_nonlipid_ultrapure3_panels.R`
7. `07_prepare_nonlipid_ultrapure8_review.py`
8. `08_ldsc_nonlipid_ultrapure8.R`
9. `09_efa_esem_nonlipid_ultrapure8.R`
10. `10_finalize_nonlipid8_model.py`
11. `11_prepare_wsl_factor_gwas_inputs_nonlipid8.R`
12. `12_run_wsl_usergwas_nonlipid8.R`
13. `13_launch_wsl_prepare_nonlipid8.ps1`
14. `14_check_wsl_prepare_nonlipid8.ps1`
15. `15_launch_wsl_usergwas_nonlipid8_3way.ps1`
16. `16_check_wsl_usergwas_nonlipid8_3way.ps1`
17. `17_merge_usergwas_nonlipid8.R`
18. `18_export_standard_usergwas_nonlipid8.R`
19. `19_ldsc_validate_factor_outputs_nonlipid8.R`
20. `20_compile_nonlipid8_supplement.R`
21. `21_build_nonlipid8_workbook.py`
22. `22_build_combined_lipid_nonlipid_workbook.py`

## Important note about manual selection

These pipelines cannot be reproduced from scripts alone because the successful routes include manual selection decisions.

Those decisions are preserved here as versioned TSV inputs and should be treated as part of the reproducible analysis input.

Examples include:

- `lipid_final8/inputs/selection/lipid_module_compact_review.tsv`
- `lipid_final8/inputs/selection/lipid_module_ultrapure3_review.tsv`
- `lipid_final8/inputs/selection/lipid_module_ultrapure3_final8_review.tsv`
- `nonlipid_final8/inputs/selection/nonlipid_module_compact_review.tsv`
- `nonlipid_final8/inputs/selection/nonlipid_module_ultrapure3_review.tsv`
- `nonlipid_final8/inputs/selection/nonlipid_group_inclusion_table.tsv`
- `nonlipid_final8/inputs/selection/nonlipid_module_grouping.tsv`

## Figure-generation code

The repository includes code for rebuilding the main lipid+nonlipid figures and the factor-GWAS Manhattan/QQ plots:

- `lipid_final8/scripts/16_build_lipid_nonlipid_result_figures.py`
- `lipid_final8/scripts/17_build_factor_gwas_manhattan_qq.py`

These scripts regenerate the combined lipid+nonlipid result figures, as well as factor-GWAS Manhattan and QQ plots for both modules.

## Path assumptions

The uploaded scripts still use the same absolute path conventions that were used in the successful run:

- GWAS root: `D:/metabolic/GWAS`
- Main analysis root: `D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion`
- LDSC reference: `D:/LDSC/ldsc-master/eur_w_ld_chr`
- WSL R library: `/home/shenjing/R/genomicsem_fix_lib`

If the target machine uses different locations, these paths must be edited consistently before rerunning.
