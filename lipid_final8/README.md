# lipid_final8

This directory contains the successful clean lipid-module workflow derived from the shared `Main_Zgt4_nonproportion` `h2 Z > 4` universe.

The retained route is:

- `112 -> 94 -> 20 -> 10 -> 8`
- shared Main112 universe -> lipid candidate set -> compact20 -> ultra-pure10 -> final lipid final8 panel

The final 8-trait model resolves three interpretable factors:

- `F1_TG_rich_axis`
- `F2_HDL_core_axis`
- `F3_CE_structural_axis`

## Contents

- `inputs/selection/`
  - Shared Main112 inputs plus lipid-specific candidate, compact, ultrapure, and final review tables.
- `inputs/final_model/`
  - Final manifest, factor structure, and model summary.
- `inputs/validation/`
  - Trait manifest and requested pair list for internal bivariate LDSC validation against final indicators and same-domain support traits.
- `scripts/`
  - Successful lipid scripts from shared step1 through factor GWAS, LDSC validation, supplement-table generation, figure builders, and internal validation.

## Notes

- The successful downstream factor-GWAS route uses WSL-native `GenomicSEM::userGWAS`.
- Internal validation is defined through `inputs/validation/validation_traits.tsv` and `inputs/validation/requested_pairs.tsv`, then executed by `18_internal_bivariate_ldsc_lipid_vs_indicators.R`.
- Failed model lines, temporary diagnostics, and abandoned Windows-only downstream routes are intentionally excluded.
