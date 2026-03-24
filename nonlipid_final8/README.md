# nonlipid_final8

This directory contains the successful clean nonlipid-module workflow derived from the shared `Main_Zgt4_nonproportion` `h2 Z > 4` universe.

The retained route is:

- `112 -> 18 -> 15 -> 8`
- shared Main112 universe -> nonlipid candidate set -> compact15 -> final nonlipid final8 panel

The final 8-trait model resolves three interpretable factors:

- `F1_ketone_axis`
- `F2_amino_acid_axis`
- `F3_energy_bridge_axis`

## Contents

- `inputs/selection/`
  - Shared Main112 inputs plus nonlipid-specific candidate, grouping, compact, and ultrapure review tables.
- `inputs/final_model/`
  - Final manifest, factor structure, and model summary.
- `scripts/`
  - Successful nonlipid scripts from shared step1 through factor GWAS, LDSC validation, supplement-table generation, and combined workbook building.
- `inputs/validation/`
  - Trait manifest and requested pair list for internal bivariate LDSC validation against final indicators and support traits.

## Notes

- `01_step1_ldsc_main112_qc.R` is shared with the lipid workflow and is copied here so the nonlipid module can be reproduced independently.
- The successful downstream factor-GWAS route uses WSL-native `GenomicSEM::userGWAS`.
- Internal validation is defined through `inputs/validation/validation_traits.tsv` and `inputs/validation/requested_pairs.tsv`, then executed by `23_internal_bivariate_ldsc_nonlipid_vs_indicators.R`.
- Failed model lines, temporary diagnostics, and abandoned Windows-only downstream routes are intentionally excluded.
