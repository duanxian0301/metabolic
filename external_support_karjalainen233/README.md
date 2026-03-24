# external_support_karjalainen233

This directory contains the clean external-validation code used after the `lipid_final8` and `nonlipid_final8` factor-GWAS analyses.

It starts from locally downloaded external GWAS summary statistics and carries the workflow through:

1. external trait cataloging
2. GWAS standardization and QC
3. univariate LDSC screening
4. external trait prioritization
5. round-1 matched-trait bivariate LDSC
6. round-2 same-domain support bivariate LDSC
7. combined external-validation workbook export

## Contents

- `inputs/`
  - versioned manifests used to reproduce the external-validation selection path
  - round-specific trait manifests and requested-pair tables
- `scripts/`
  - clean scripts for standardization, retries, univariate LDSC, pair prioritization, bivariate LDSC, result collection, and workbook export

## Scope

This directory is intentionally code-first and manifest-first.

It does not include:

- raw downloaded external GWAS files
- generated QC outputs
- LDSC result directories
- supplementary workbooks generated from runtime outputs

## Script flow

1. `01_build_external_ldsc_tables.py`
2. `02_build_external_processing_manifest.py`
3. `03_standardize_external_gwas_to_txt.py`
4. `04_launch_external233_standardization_3way.ps1`
5. `05_check_external233_standardization_3way.ps1`
6. `06_build_failed_retry_manifest.py`
7. `07_launch_failed_retry_batches.ps1`
8. `08_build_tail_retry_manifest.py`
9. `09_launch_tail_retry_batches.ps1`
10. `10_build_external233_ldsc_manifest.py`
11. `11_run_external233_univariate_ldsc.py`
12. `12_collect_external233_univariate_ldsc.py`
13. `13_launch_external233_univariate_ldsc_3way.ps1`
14. `14_check_external233_univariate_ldsc_3way.ps1`
15. `15_run_external233_univariate_ldsc_wsl.py`
16. `16_launch_external233_univariate_ldsc_wsl.ps1`
17. `17_build_external_validation_priority_tables.py`
18. `18_build_external_bivariate_round1_manifests.py`
19. `19_run_external_bivariate_lipid_round1.R`
20. `20_run_external_bivariate_nonlipid_round1.R`
21. `21_collect_external_bivariate_round1.py`
22. `22_build_external_bivariate_round2_manifests.py`
23. `23_run_external_bivariate_lipid_round2.R`
24. `24_run_external_bivariate_nonlipid_round2.R`
25. `25_collect_external_bivariate_round1_round2.py`
26. `26_build_external_validation_combined_workbook.py`
27. `27_build_external_validation_combined_workbook_with_neff.py`

## Notes

- The current scripts preserve the original absolute-path conventions used in the successful run.
- Factor-GWAS `N` labels were later replaced with factor-specific `Neff` estimates; the combined workbook can be rebuilt with the updated `N` values via script 27 after running the shared root-level `estimate_factor_neff_and_update_standard_txt.py`.
- The reusable, project-agnostic version of the same post-GWAS workflow is maintained separately in the `postGWAS` repository.
