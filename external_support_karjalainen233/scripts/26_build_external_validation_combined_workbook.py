import shutil
from pathlib import Path

import pandas as pd
from openpyxl import load_workbook


BASE_WORKBOOK = Path(
    r"D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion/supplement_combined_lipid_nonlipid_final8/"
    r"metabolic_lipid_nonlipid_final8_combined_supplement_workbook.xlsx"
)
OUT_DIR = BASE_WORKBOOK.parent
OUT_WORKBOOK = OUT_DIR / "metabolic_lipid_nonlipid_final8_combined_supplement_workbook_with_external_validation.xlsx"

UNIVARIATE = Path(r"D:/metabolic/233/ldsc_univariate/external233_univariate_ldsc_merged.tsv")
PRIORITY = Path(r"D:/metabolic/233/ldsc_univariate/external_validation_priority/lipid_nonlipid_external_validation_priority.tsv")
MATCHED = Path(r"D:/metabolic/233/ldsc_univariate/external_bivariate_final_summary/external_bivariate_matched_pairs_round1_round2.tsv")
ALL_PAIRS = Path(r"D:/metabolic/233/ldsc_univariate/external_bivariate_final_summary/external_bivariate_all_pairs_round1_round2.tsv")


def build_sheet_27():
    dt = pd.read_csv(UNIVARIATE, sep="\t")
    keep = [
        "study_accession",
        "trait",
        "raw_url",
        "sample_total_n",
        "discovery_sample_ancestry",
        "raw_position_build",
        "status",
        "ldsc_status",
        "n_raw_rows",
        "n_variant_id_mismatch",
        "n_invalid_rsid",
        "n_retained_pre_dedup",
        "n_duplicates_removed",
        "n_final",
        "n_sumstats_read",
        "n_after_ref_merge",
        "n_after_weight_merge",
        "h2",
        "h2_se",
        "h2_z",
        "intercept",
        "intercept_se",
        "ratio",
        "ratio_se",
        "lambda_gc",
        "mean_chi2",
        "analysis_finished",
        "time_elapsed",
    ]
    return dt[keep].copy()


def build_sheet_28():
    priority = pd.read_csv(PRIORITY, sep="\t")
    matched = pd.read_csv(MATCHED, sep="\t")

    round_map = (
        matched.groupby(["module", "study_accession"], as_index=False)["round"]
        .agg(lambda x: ",".join(sorted(set(x))))
        .rename(columns={"round": "selected_rounds"})
    )
    priority = priority.merge(round_map, on=["module", "study_accession"], how="left")
    priority["selected_for_external_bivariate"] = priority["selected_rounds"].notna()
    priority["selected_round1"] = priority["selected_rounds"].fillna("").str.contains("round1")
    priority["selected_round2"] = priority["selected_rounds"].fillna("").str.contains("round2")

    keep = [
        "module",
        "anchor_factor",
        "internal_trait",
        "external_trait",
        "study_accession",
        "subset",
        "match_quality",
        "validation_priority",
        "selected_for_external_bivariate",
        "selected_round1",
        "selected_round2",
        "selected_rounds",
        "recommendation",
        "sample_total_n",
        "discovery_sample_ancestry",
        "h2",
        "h2_se",
        "h2_z",
        "intercept",
        "ratio",
        "lambda_gc",
        "mean_chi2",
        "rationale",
    ]
    out = priority[keep].copy()
    out = out.sort_values(
        ["module", "anchor_factor", "selected_for_external_bivariate", "validation_priority", "study_accession"],
        ascending=[True, True, False, True, True],
    ).reset_index(drop=True)
    return out


def build_sheet_29():
    dt = pd.read_csv(ALL_PAIRS, sep="\t")
    keep = [
        "round",
        "module",
        "factor",
        "pair_set",
        "expected_alignment",
        "target_trait",
        "study_accession.x",
        "external_trait",
        "internal_trait.x",
        "validation_priority.x",
        "covariance",
        "covariance_se",
        "p_cov",
        "rg",
        "rg_se",
        "p_rg",
        "aligned_rg",
        "orientation_multiplier",
        "abs_rg",
    ]
    out = dt[keep].copy()
    out = out.rename(
        columns={
            "study_accession.x": "study_accession",
            "internal_trait.x": "internal_trait",
            "validation_priority.x": "validation_priority",
        }
    )
    out = out.sort_values(["round", "module", "factor", "pair_set", "study_accession"]).reset_index(drop=True)
    return out


def update_index(writer):
    index_df = pd.read_excel(BASE_WORKBOOK, sheet_name="00_Index")
    additions = pd.DataFrame(
        [
            {
                "sheet_name": "27_ExtData_QC_UniLDSC",
                "description": "External support source, standardization QC, and univariate LDSC summary merged into one sheet.",
                "Unnamed: 2": None,
                "notes": "Rows correspond to the 233 external-support metabolite GWAS from Karjalainen et al. 2024 and include data-source fields, QC retention metrics, and single-trait LDSC estimates.",
            },
            {
                "sheet_name": "28_ExtPairSelection",
                "description": "External validation pair-selection table with match tier, recommendation, and round-selection flags.",
                "Unnamed: 2": None,
                "notes": "This sheet records exact, proxy, same-domain-support, and exploratory candidates, plus which traits were selected for external bivariate LDSC in rounds 1 and 2.",
            },
            {
                "sheet_name": "29_ExtBivariateLDSC",
                "description": "External bivariate LDSC results across rounds 1 and 2 with raw and orientation-aligned genetic correlations.",
                "Unnamed: 2": None,
                "notes": "Primary matches, support matches, and cross-factor checks are combined so external validation strength and specificity can be reviewed in one place.",
            },
        ]
    )
    index_df = pd.concat([index_df, additions], ignore_index=True)
    index_df.to_excel(writer, sheet_name="00_Index", index=False)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(BASE_WORKBOOK, OUT_WORKBOOK)

    sheet27 = build_sheet_27()
    sheet28 = build_sheet_28()
    sheet29 = build_sheet_29()

    with pd.ExcelWriter(OUT_WORKBOOK, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        update_index(writer)
        sheet27.to_excel(writer, sheet_name="27_ExtData_QC_UniLDSC", index=False)
        sheet28.to_excel(writer, sheet_name="28_ExtPairSelection", index=False)
        sheet29.to_excel(writer, sheet_name="29_ExtBivariateLDSC", index=False)

    print(OUT_WORKBOOK)


if __name__ == "__main__":
    main()
