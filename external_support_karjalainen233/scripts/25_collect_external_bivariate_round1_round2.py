import pandas as pd
from pathlib import Path


ROOT = Path(r"D:/metabolic/233/ldsc_univariate")


def load_pairs(round_label, module_name):
    if round_label == "round1" and module_name == "lipid":
        path = ROOT / "external_bivariate_round1_lipid" / "external_bivariate_lipid_round1_requested_pairs.tsv"
    elif round_label == "round1" and module_name == "nonlipid":
        path = ROOT / "external_bivariate_round1_nonlipid" / "external_bivariate_nonlipid_round1_requested_pairs.tsv"
    elif round_label == "round2" and module_name == "lipid":
        path = ROOT / "external_bivariate_round2_lipid" / "external_bivariate_lipid_round2_requested_pairs.tsv"
    else:
        path = ROOT / "external_bivariate_round2_nonlipid" / "external_bivariate_nonlipid_round2_requested_pairs.tsv"
    dt = pd.read_csv(path, sep="\t")
    dt["module"] = module_name
    dt["round"] = round_label
    return dt


def factor_orientation(all_pairs):
    primary_like = all_pairs[all_pairs["pair_set"].isin(["primary_match", "support_match"])].copy()
    orientation = (
        primary_like.groupby(["module", "factor"], as_index=False)["rg"]
        .mean()
        .rename(columns={"rg": "reference_mean_rg"})
    )
    orientation["orientation_multiplier"] = orientation["reference_mean_rg"].apply(lambda x: -1 if x < 0 else 1)
    return orientation


def main():
    out_dir = ROOT / "external_bivariate_final_summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    all_pairs = pd.concat(
        [
            load_pairs("round1", "lipid"),
            load_pairs("round1", "nonlipid"),
            load_pairs("round2", "lipid"),
            load_pairs("round2", "nonlipid"),
        ],
        ignore_index=True,
    )

    orient = factor_orientation(all_pairs)
    all_pairs = all_pairs.merge(orient, on=["module", "factor"], how="left", validate="many_to_one")
    all_pairs["aligned_rg"] = all_pairs["rg"] * all_pairs["orientation_multiplier"]
    all_pairs["aligned_abs_rg"] = all_pairs["aligned_rg"].abs()

    final_matches = all_pairs[all_pairs["pair_set"].isin(["primary_match", "support_match"])].copy()
    final_matches = final_matches[
        [
            "round",
            "module",
            "factor",
            "pair_set",
            "target_trait",
            "study_accession.x",
            "external_trait",
            "internal_trait.x",
            "validation_priority.x",
            "rg",
            "rg_se",
            "p_rg",
            "aligned_rg",
            "orientation_multiplier",
        ]
    ].rename(
        columns={
            "study_accession.x": "study_accession",
            "internal_trait.x": "internal_trait",
            "validation_priority.x": "validation_priority",
        }
    )

    cross_factor = all_pairs[all_pairs["pair_set"] == "cross_factor"].copy()
    cross_factor = cross_factor[
        [
            "round",
            "module",
            "factor",
            "anchor_factor.x",
            "target_trait",
            "study_accession.x",
            "external_trait",
            "internal_trait.x",
            "rg",
            "rg_se",
            "p_rg",
            "aligned_rg",
        ]
    ].rename(
        columns={
            "anchor_factor.x": "anchor_factor",
            "study_accession.x": "study_accession",
            "internal_trait.x": "internal_trait",
        }
    )

    factor_summary = (
        all_pairs.groupby(["round", "module", "factor", "pair_set"], as_index=False)
        .agg(
            n_pairs=("target_trait", "size"),
            mean_rg=("rg", "mean"),
            mean_abs_rg=("abs_rg", "mean"),
            mean_aligned_rg=("aligned_rg", "mean"),
            median_aligned_rg=("aligned_rg", "median"),
            min_p_rg=("p_rg", "min"),
        )
    )

    workbook = out_dir / "external_bivariate_final_summary.xlsx"
    directory = pd.DataFrame(
        {
            "sheet": [
                "directory",
                "all_pairs",
                "matched_pairs",
                "cross_factor",
                "factor_summary",
                "factor_orientation",
            ],
            "description": [
                "Workbook directory",
                "All external bivariate LDSC pairs from rounds 1 and 2",
                "Primary and support matches with aligned rg",
                "Cross-factor specificity checks",
                "Factor-level summary across rounds and pair sets",
                "Orientation multipliers used for aligned interpretation",
            ],
        }
    )

    all_pairs.to_csv(out_dir / "external_bivariate_all_pairs_round1_round2.tsv", sep="\t", index=False)
    final_matches.to_csv(out_dir / "external_bivariate_matched_pairs_round1_round2.tsv", sep="\t", index=False)
    cross_factor.to_csv(out_dir / "external_bivariate_cross_factor_round1_round2.tsv", sep="\t", index=False)
    factor_summary.to_csv(out_dir / "external_bivariate_factor_summary_round1_round2.tsv", sep="\t", index=False)
    orient.to_csv(out_dir / "external_bivariate_factor_orientation_round1_round2.tsv", sep="\t", index=False)

    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        directory.to_excel(writer, sheet_name="directory", index=False)
        all_pairs.to_excel(writer, sheet_name="all_pairs", index=False)
        final_matches.to_excel(writer, sheet_name="matched_pairs", index=False)
        cross_factor.to_excel(writer, sheet_name="cross_factor", index=False)
        factor_summary.to_excel(writer, sheet_name="factor_summary", index=False)
        orient.to_excel(writer, sheet_name="factor_orientation", index=False)

    print(out_dir)


if __name__ == "__main__":
    main()
