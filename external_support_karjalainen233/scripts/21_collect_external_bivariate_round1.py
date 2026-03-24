import pandas as pd
from pathlib import Path


ROOT = Path(r"D:/metabolic/233/ldsc_univariate")


def load_pairs(module_name):
    if module_name == "lipid":
        path = ROOT / "external_bivariate_round1_lipid" / "external_bivariate_lipid_round1_requested_pairs.tsv"
    else:
        path = ROOT / "external_bivariate_round1_nonlipid" / "external_bivariate_nonlipid_round1_requested_pairs.tsv"
    dt = pd.read_csv(path, sep="\t")
    dt["module"] = module_name
    return dt


def load_h2(module_name):
    if module_name == "lipid":
        path = ROOT / "external_bivariate_round1_lipid" / "external_bivariate_lipid_round1_h2.tsv"
    else:
        path = ROOT / "external_bivariate_round1_nonlipid" / "external_bivariate_nonlipid_round1_h2.tsv"
    dt = pd.read_csv(path, sep="\t")
    dt["module"] = module_name
    return dt


def factor_orientation(pairs):
    primary = (
        pairs[pairs["pair_set"] == "primary_match"]
        .groupby(["module", "factor"], as_index=False)["rg"]
        .mean()
        .rename(columns={"rg": "primary_mean_rg"})
    )
    primary["orientation_multiplier"] = primary["primary_mean_rg"].apply(lambda x: -1 if x < 0 else 1)
    return primary


def main():
    out_dir = ROOT / "external_bivariate_round1_summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    pairs = pd.concat([load_pairs("lipid"), load_pairs("nonlipid")], ignore_index=True)
    h2 = pd.concat([load_h2("lipid"), load_h2("nonlipid")], ignore_index=True)
    orient = factor_orientation(pairs)

    pairs = pairs.merge(orient, on=["module", "factor"], how="left", validate="many_to_one")
    pairs["aligned_rg"] = pairs["rg"] * pairs["orientation_multiplier"]
    pairs["aligned_abs_rg"] = pairs["aligned_rg"].abs()

    factor_summary = (
        pairs.groupby(["module", "factor", "pair_set"], as_index=False)
        .agg(
            n_pairs=("target_trait", "size"),
            mean_rg=("rg", "mean"),
            mean_abs_rg=("abs_rg", "mean"),
            mean_aligned_rg=("aligned_rg", "mean"),
            median_aligned_rg=("aligned_rg", "median"),
            min_p_rg=("p_rg", "min"),
        )
    )

    primary_view = pairs[pairs["pair_set"] == "primary_match"].copy()
    primary_view = primary_view[
        [
            "module",
            "factor",
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

    cross_view = pairs[pairs["pair_set"] == "cross_factor"].copy()
    cross_view = cross_view[
        [
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

    workbook = out_dir / "external_bivariate_round1_summary.xlsx"
    directory = pd.DataFrame(
        {
            "sheet": [
                "directory",
                "pair_results_all",
                "primary_matches",
                "cross_factor",
                "factor_summary",
                "trait_h2",
                "factor_orientation",
            ],
            "description": [
                "Workbook directory",
                "All pairwise round1 external bivariate LDSC results",
                "Primary factor-trait matches with raw and aligned rg",
                "Cross-factor checks for specificity",
                "Factor-level summary across primary vs cross-factor pairs",
                "Trait-level h2/intercept results from round1 bivariate runs",
                "Orientation multipliers used to align factor signs for interpretation",
            ],
        }
    )

    pairs.to_csv(out_dir / "external_bivariate_round1_all_pairs.tsv", sep="\t", index=False)
    primary_view.to_csv(out_dir / "external_bivariate_round1_primary_matches.tsv", sep="\t", index=False)
    cross_view.to_csv(out_dir / "external_bivariate_round1_cross_factor.tsv", sep="\t", index=False)
    factor_summary.to_csv(out_dir / "external_bivariate_round1_factor_summary.tsv", sep="\t", index=False)
    h2.to_csv(out_dir / "external_bivariate_round1_h2.tsv", sep="\t", index=False)
    orient.to_csv(out_dir / "external_bivariate_round1_factor_orientation.tsv", sep="\t", index=False)

    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        directory.to_excel(writer, sheet_name="directory", index=False)
        pairs.to_excel(writer, sheet_name="pair_results_all", index=False)
        primary_view.to_excel(writer, sheet_name="primary_matches", index=False)
        cross_view.to_excel(writer, sheet_name="cross_factor", index=False)
        factor_summary.to_excel(writer, sheet_name="factor_summary", index=False)
        h2.to_excel(writer, sheet_name="trait_h2", index=False)
        orient.to_excel(writer, sheet_name="factor_orientation", index=False)

    print(out_dir)


if __name__ == "__main__":
    main()
