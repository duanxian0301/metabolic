import pandas as pd
from pathlib import Path


ROOT = Path(r"D:/metabolic/233/ldsc_univariate")
INPUT_ROOT = Path(r"D:/codex/metabolic_repo/external_support_karjalainen233/inputs")


def load_external_ldsc():
    return pd.read_csv(ROOT / "external233_univariate_ldsc_merged.tsv", sep="\t")


def classify_row(row):
    subset = row.get("subset", "")
    match_quality = str(row.get("match_quality", ""))
    h2_z = row.get("h2_z")
    ratio = row.get("ratio")
    intercept = row.get("intercept")

    if pd.isna(h2_z) or pd.isna(ratio) or pd.isna(intercept):
        return "incomplete_metrics"

    if subset == "external_exploratory_bridge":
        return "exploratory_only"

    if subset == "external_primary_proxy":
        if h2_z >= 4 and ratio <= 0.2:
            return "proxy_bivariate_candidate"
        return "proxy_support_only"

    if subset == "external_same_domain_support":
        if h2_z >= 4 and ratio <= 0.2 and intercept <= 1.05:
            return "same_domain_bivariate_candidate"
        return "same_domain_support_only"

    if subset == "external_primary_exact":
        if h2_z >= 4 and ratio <= 0.2 and intercept <= 1.05:
            return "primary_bivariate_candidate"
        if h2_z >= 4:
            return "primary_bivariate_candidate_caution"
        return "primary_support_only"

    if "exact" in match_quality and h2_z >= 4 and ratio <= 0.2 and intercept <= 1.05:
        return "primary_bivariate_candidate"
    return "support_only"


def recommendation_text(row):
    rec = row["validation_priority"]
    if rec == "primary_bivariate_candidate":
        return "Recommended for external bivariate LDSC as a primary matched trait."
    if rec == "primary_bivariate_candidate_caution":
        return "Usable for bivariate LDSC, but report caution because intercept or attenuation ratio is elevated."
    if rec == "same_domain_bivariate_candidate":
        return "Suitable for same-domain external bivariate LDSC support analysis."
    if rec == "proxy_bivariate_candidate":
        return "Usable as a proxy-trait bivariate LDSC analysis, but do not overstate as exact replication."
    if rec in {"primary_support_only", "same_domain_support_only", "proxy_support_only", "support_only"}:
        return "Better treated as external support rather than a headline bivariate LDSC validation target."
    if rec == "exploratory_only":
        return "Exploratory bridge trait only; keep outside the primary external validation set."
    return "Review manually."


def build_one(module_name, path, ldsc):
    dt = pd.read_csv(path, sep="\t")
    ext = dt[dt["category"] == "external_trait"].copy()
    ext = ext.merge(
        ldsc[
            [
                "study_accession",
                "trait",
                "sample_total_n",
                "discovery_sample_ancestry",
                "h2",
                "h2_se",
                "h2_z",
                "intercept",
                "intercept_se",
                "ratio",
                "ratio_se",
                "lambda_gc",
                "mean_chi2",
                "status",
                "ldsc_status",
            ]
        ].rename(columns={"trait": "external_trait_ldsc"}),
        on="study_accession",
        how="left",
        validate="many_to_one",
    )
    ext["module"] = module_name
    ext["validation_priority"] = ext.apply(classify_row, axis=1)
    ext["recommendation"] = ext.apply(recommendation_text, axis=1)
    return ext


def main():
    ldsc = load_external_ldsc()
    lipid = build_one("lipid", INPUT_ROOT / "lipid_final8_external_validation_traits.tsv", ldsc)
    nonlipid = build_one("nonlipid", INPUT_ROOT / "nonlipid_final8_external_validation_traits.tsv", ldsc)
    combined = pd.concat([lipid, nonlipid], ignore_index=True)

    priority_order = [
        "primary_bivariate_candidate",
        "primary_bivariate_candidate_caution",
        "same_domain_bivariate_candidate",
        "proxy_bivariate_candidate",
        "primary_support_only",
        "same_domain_support_only",
        "proxy_support_only",
        "exploratory_only",
        "support_only",
        "incomplete_metrics",
    ]
    combined["validation_priority"] = pd.Categorical(
        combined["validation_priority"],
        categories=priority_order,
        ordered=True,
    )
    combined = combined.sort_values(
        ["module", "anchor_factor", "validation_priority", "subset", "study_accession"]
    ).reset_index(drop=True)

    out_dir = ROOT / "external_validation_priority"
    out_dir.mkdir(parents=True, exist_ok=True)

    lipid.to_csv(out_dir / "lipid_final8_external_validation_priority.tsv", sep="\t", index=False)
    nonlipid.to_csv(out_dir / "nonlipid_final8_external_validation_priority.tsv", sep="\t", index=False)
    combined.to_csv(out_dir / "lipid_nonlipid_external_validation_priority.tsv", sep="\t", index=False)

    summary = (
        combined.groupby(["module", "anchor_factor", "validation_priority"], dropna=False, observed=True)
        .size()
        .reset_index(name="n_traits")
        .query("n_traits > 0")
        .sort_values(["module", "anchor_factor", "validation_priority"])
    )
    summary.to_csv(out_dir / "external_validation_priority_summary.tsv", sep="\t", index=False)

    directory = pd.DataFrame(
        {
            "sheet": [
                "directory",
                "combined_priority",
                "lipid_priority",
                "nonlipid_priority",
                "priority_summary",
            ],
            "description": [
                "Workbook directory",
                "Combined lipid and nonlipid external validation priority table",
                "Lipid-only external validation priority table",
                "Nonlipid-only external validation priority table",
                "Counts by module, factor, and priority tier",
            ],
        }
    )
    workbook = out_dir / "lipid_nonlipid_external_validation_priority.xlsx"
    with pd.ExcelWriter(workbook, engine="openpyxl") as writer:
        directory.to_excel(writer, sheet_name="directory", index=False)
        combined.to_excel(writer, sheet_name="combined_priority", index=False)
        lipid.to_excel(writer, sheet_name="lipid_priority", index=False)
        nonlipid.to_excel(writer, sheet_name="nonlipid_priority", index=False)
        summary.to_excel(writer, sheet_name="priority_summary", index=False)

    print(out_dir)


if __name__ == "__main__":
    main()
