from pathlib import Path

import pandas as pd


REPO_ROOT = Path(r"D:/codex/metabolic_repo")
EXCEL_PATH = Path(r"D:/metabolic/Metabolic-data.xlsx")
OUTPUT_DIR = REPO_ROOT / "external_support_karjalainen233" / "inputs"
RAW_BASE = Path(r"D:/metabolic")
MUNGED_BASE = Path(r"D:/metabolic/external_support_karjalainen233/sumstats")
DATASET = "karjalainen_nature_2024_233_metabolites"
DATASET_NOTE = (
    "Public GWAS Catalog summary statistics from Karjalainen et al. Nature 2024; "
    "trait-specific files are predominantly multi-ancestry meta-analysis rather than "
    "European-only releases."
)


LIPID_MAPPINGS = [
    {
        "anchor_factor": "lipid8_F1",
        "anchor_factor_label": "F1_TG_rich_axis",
        "internal_trait": "M_HDL_TG",
        "internal_biomarker_name": "Triglycerides in medium HDL",
        "study_accession": "GCST90302041",
        "external_trait": "Triglycerides in medium HDL",
        "subset": "external_primary_exact",
        "match_quality": "exact",
        "expected_alignment": "primary",
        "rationale": "Direct matched external trait for medium-HDL triglycerides.",
    },
    {
        "anchor_factor": "lipid8_F1",
        "anchor_factor_label": "F1_TG_rich_axis",
        "internal_trait": "VLDL_size",
        "internal_biomarker_name": "Average diameter for VLDL particles",
        "study_accession": "GCST90302124",
        "external_trait": "Mean diameter of VLDL particles",
        "subset": "external_primary_exact",
        "match_quality": "exact_terminology_variant",
        "expected_alignment": "primary",
        "rationale": "Same VLDL particle-size construct with wording differences only.",
    },
    {
        "anchor_factor": "lipid8_F1",
        "anchor_factor_label": "F1_TG_rich_axis",
        "internal_trait": "MUFA",
        "internal_biomarker_name": "Monounsaturated fatty acids",
        "study_accession": "GCST90302055",
        "external_trait": "Monounsaturated fatty acids (16:1, 18:1) levels",
        "subset": "external_primary_exact",
        "match_quality": "exact_platform_variant",
        "expected_alignment": "primary",
        "rationale": "Direct MUFA trait on the same Nightingale biomarker platform.",
    },
    {
        "anchor_factor": "lipid8_F1",
        "anchor_factor_label": "F1_TG_rich_axis",
        "internal_trait": "S_VLDL_TG",
        "internal_biomarker_name": "Triglycerides in small VLDL",
        "study_accession": "GCST90302114",
        "external_trait": "Triglycerides in small VLDL",
        "subset": "external_primary_exact",
        "match_quality": "exact",
        "expected_alignment": "primary",
        "rationale": "Direct matched external trait for small-VLDL triglycerides.",
    },
    {
        "anchor_factor": "lipid8_F2",
        "anchor_factor_label": "F2_HDL_core_axis",
        "internal_trait": "ApoA1",
        "internal_biomarker_name": "Apolipoprotein A1",
        "study_accession": "GCST90301945",
        "external_trait": "Apolipoprotein A-I levels",
        "subset": "external_primary_exact",
        "match_quality": "exact_terminology_variant",
        "expected_alignment": "primary",
        "rationale": "Direct matched ApoA1 biomarker with naming normalization only.",
    },
    {
        "anchor_factor": "lipid8_F2",
        "anchor_factor_label": "F2_HDL_core_axis",
        "internal_trait": "HDL_CE",
        "internal_biomarker_name": "Cholesteryl esters in HDL",
        "study_accession": "GCST90301997",
        "external_trait": "Cholesterol esters in large HDL",
        "subset": "external_same_domain_support",
        "match_quality": "approximate_size_specific",
        "expected_alignment": "same_domain_support",
        "rationale": "Overall HDL-CE trait is unavailable; large-HDL CE used as size-specific support.",
    },
    {
        "anchor_factor": "lipid8_F2",
        "anchor_factor_label": "F2_HDL_core_axis",
        "internal_trait": "HDL_CE",
        "internal_biomarker_name": "Cholesteryl esters in HDL",
        "study_accession": "GCST90302033",
        "external_trait": "Cholesterol esters in medium HDL",
        "subset": "external_same_domain_support",
        "match_quality": "approximate_size_specific",
        "expected_alignment": "same_domain_support",
        "rationale": "Overall HDL-CE trait is unavailable; medium-HDL CE used as size-specific support.",
    },
    {
        "anchor_factor": "lipid8_F2",
        "anchor_factor_label": "F2_HDL_core_axis",
        "internal_trait": "HDL_CE",
        "internal_biomarker_name": "Cholesteryl esters in HDL",
        "study_accession": "GCST90302081",
        "external_trait": "Cholesterol esters in small HDL",
        "subset": "external_same_domain_support",
        "match_quality": "approximate_size_specific",
        "expected_alignment": "same_domain_support",
        "rationale": "Overall HDL-CE trait is unavailable; small-HDL CE used as size-specific support.",
    },
    {
        "anchor_factor": "lipid8_F2",
        "anchor_factor_label": "F2_HDL_core_axis",
        "internal_trait": "HDL_CE",
        "internal_biomarker_name": "Cholesteryl esters in HDL",
        "study_accession": "GCST90302128",
        "external_trait": "Cholesterol esters in very large HDL",
        "subset": "external_same_domain_support",
        "match_quality": "approximate_size_specific",
        "expected_alignment": "same_domain_support",
        "rationale": "Overall HDL-CE trait is unavailable; very-large-HDL CE used as size-specific support.",
    },
    {
        "anchor_factor": "lipid8_F3",
        "anchor_factor_label": "F3_CE_structural_axis",
        "internal_trait": "XS_VLDL_FC",
        "internal_biomarker_name": "Free cholesterol in very small VLDL",
        "study_accession": "GCST90302154",
        "external_trait": "Free cholesterol in very small VLDL",
        "subset": "external_primary_exact",
        "match_quality": "exact",
        "expected_alignment": "primary",
        "rationale": "Direct matched external trait for very-small-VLDL free cholesterol.",
    },
    {
        "anchor_factor": "lipid8_F3",
        "anchor_factor_label": "F3_CE_structural_axis",
        "internal_trait": "VLDL_CE",
        "internal_biomarker_name": "Cholesteryl esters in VLDL",
        "study_accession": "GCST90302021",
        "external_trait": "Cholesterol esters in large VLDL",
        "subset": "external_same_domain_support",
        "match_quality": "approximate_size_specific",
        "expected_alignment": "same_domain_support",
        "rationale": "Overall VLDL-CE trait is unavailable; large-VLDL CE used as size-specific support.",
    },
    {
        "anchor_factor": "lipid8_F3",
        "anchor_factor_label": "F3_CE_structural_axis",
        "internal_trait": "VLDL_CE",
        "internal_biomarker_name": "Cholesteryl esters in VLDL",
        "study_accession": "GCST90302059",
        "external_trait": "Cholesterol esters in medium VLDL",
        "subset": "external_same_domain_support",
        "match_quality": "approximate_size_specific",
        "expected_alignment": "same_domain_support",
        "rationale": "Overall VLDL-CE trait is unavailable; medium-VLDL CE used as size-specific support.",
    },
    {
        "anchor_factor": "lipid8_F3",
        "anchor_factor_label": "F3_CE_structural_axis",
        "internal_trait": "VLDL_CE",
        "internal_biomarker_name": "Cholesteryl esters in VLDL",
        "study_accession": "GCST90302106",
        "external_trait": "Cholesterol esters in small VLDL",
        "subset": "external_same_domain_support",
        "match_quality": "approximate_size_specific",
        "expected_alignment": "same_domain_support",
        "rationale": "Overall VLDL-CE trait is unavailable; small-VLDL CE used as size-specific support.",
    },
    {
        "anchor_factor": "lipid8_F3",
        "anchor_factor_label": "F3_CE_structural_axis",
        "internal_trait": "VLDL_CE",
        "internal_biomarker_name": "Cholesteryl esters in VLDL",
        "study_accession": "GCST90302140",
        "external_trait": "Cholesterol esters in very large VLDL",
        "subset": "external_same_domain_support",
        "match_quality": "approximate_size_specific",
        "expected_alignment": "same_domain_support",
        "rationale": "Overall VLDL-CE trait is unavailable; very-large-VLDL CE used as size-specific support.",
    },
    {
        "anchor_factor": "lipid8_F3",
        "anchor_factor_label": "F3_CE_structural_axis",
        "internal_trait": "VLDL_CE",
        "internal_biomarker_name": "Cholesteryl esters in VLDL",
        "study_accession": "GCST90302152",
        "external_trait": "Cholesterol esters in very small VLDL",
        "subset": "external_same_domain_support",
        "match_quality": "approximate_size_specific",
        "expected_alignment": "same_domain_support",
        "rationale": "Overall VLDL-CE trait is unavailable; very-small-VLDL CE used as size-specific support.",
    },
]


NONLIPID_MAPPINGS = [
    {
        "anchor_factor": "nonlipid8_F1",
        "anchor_factor_label": "F1_ketone_axis",
        "internal_trait": "bOHbutyrate",
        "internal_biomarker_name": "3-hydroxybutyrate",
        "study_accession": "GCST90301948",
        "external_trait": "3-Hydroxybutyrate levels",
        "subset": "external_primary_exact",
        "match_quality": "exact_terminology_variant",
        "expected_alignment": "primary",
        "rationale": "Direct matched ketone-body biomarker with terminology normalization only.",
    },
    {
        "anchor_factor": "nonlipid8_F1",
        "anchor_factor_label": "F1_ketone_axis",
        "internal_trait": "Acetoacetate",
        "internal_biomarker_name": "Acetoacetate",
        "study_accession": "GCST90301942",
        "external_trait": "Acetone levels",
        "subset": "external_primary_proxy",
        "match_quality": "proxy_same_ketone_pathway",
        "expected_alignment": "primary_proxy",
        "rationale": "Acetoacetate is not available; acetone serves as a ketone-body proxy in the same pathway.",
    },
    {
        "anchor_factor": "nonlipid8_F2",
        "anchor_factor_label": "F2_amino_acid_axis",
        "internal_trait": "Leu",
        "internal_biomarker_name": "Leucine",
        "study_accession": "GCST90301994",
        "external_trait": "Leucine levels",
        "subset": "external_primary_exact",
        "match_quality": "exact_terminology_variant",
        "expected_alignment": "primary",
        "rationale": "Direct matched amino-acid biomarker with spelling normalization only.",
    },
    {
        "anchor_factor": "nonlipid8_F2",
        "anchor_factor_label": "F2_amino_acid_axis",
        "internal_trait": "Phe",
        "internal_biomarker_name": "Phenylalanine",
        "study_accession": "GCST90302070",
        "external_trait": "Phenylalanine levels",
        "subset": "external_primary_exact",
        "match_quality": "exact",
        "expected_alignment": "primary",
        "rationale": "Direct matched aromatic-amino-acid biomarker.",
    },
    {
        "anchor_factor": "nonlipid8_F2",
        "anchor_factor_label": "F2_amino_acid_axis",
        "internal_trait": "Val",
        "internal_biomarker_name": "Valine",
        "study_accession": "GCST90302122",
        "external_trait": "Valine levels",
        "subset": "external_primary_exact",
        "match_quality": "exact",
        "expected_alignment": "primary",
        "rationale": "Direct matched branched-chain-amino-acid biomarker.",
    },
    {
        "anchor_factor": "nonlipid8_F2",
        "anchor_factor_label": "F2_amino_acid_axis",
        "internal_trait": "Ala",
        "internal_biomarker_name": "Alanine",
        "study_accession": "GCST90301943",
        "external_trait": "Alanine levels",
        "subset": "external_same_domain_support",
        "match_quality": "same_domain_support",
        "expected_alignment": "same_domain_support",
        "rationale": "Additional amino-acid support trait from the same nonlipid domain.",
    },
    {
        "anchor_factor": "nonlipid8_F2",
        "anchor_factor_label": "F2_amino_acid_axis",
        "internal_trait": "Ile",
        "internal_biomarker_name": "Isoleucine",
        "study_accession": "GCST90301987",
        "external_trait": "Isoleucine levels",
        "subset": "external_same_domain_support",
        "match_quality": "same_domain_support",
        "expected_alignment": "same_domain_support",
        "rationale": "Additional BCAA support trait from the same amino-acid domain.",
    },
    {
        "anchor_factor": "nonlipid8_F3",
        "anchor_factor_label": "F3_energy_bridge_axis",
        "internal_trait": "Acetate",
        "internal_biomarker_name": "Acetate",
        "study_accession": "GCST90301941",
        "external_trait": "Acetate levels",
        "subset": "external_primary_exact",
        "match_quality": "exact",
        "expected_alignment": "primary",
        "rationale": "Direct matched energy-metabolism biomarker.",
    },
    {
        "anchor_factor": "nonlipid8_F3",
        "anchor_factor_label": "F3_energy_bridge_axis",
        "internal_trait": "Glucose",
        "internal_biomarker_name": "Glucose",
        "study_accession": "GCST90301964",
        "external_trait": "Glucose levels",
        "subset": "external_primary_exact",
        "match_quality": "exact",
        "expected_alignment": "primary",
        "rationale": "Direct matched glycolysis-related biomarker.",
    },
    {
        "anchor_factor": "nonlipid8_F3",
        "anchor_factor_label": "F3_energy_bridge_axis",
        "internal_trait": "Lactate",
        "internal_biomarker_name": "Lactate",
        "study_accession": "GCST90301990",
        "external_trait": "Lactate levels",
        "subset": "external_primary_exact",
        "match_quality": "exact",
        "expected_alignment": "primary",
        "rationale": "Direct matched glycolysis-related biomarker.",
    },
    {
        "anchor_factor": "nonlipid8_F3",
        "anchor_factor_label": "F3_energy_bridge_axis",
        "internal_trait": "Citrate",
        "internal_biomarker_name": "Citrate",
        "study_accession": "GCST90301949",
        "external_trait": "Citrate levels",
        "subset": "external_same_domain_support",
        "match_quality": "same_domain_support",
        "expected_alignment": "same_domain_support",
        "rationale": "Additional TCA-cycle support trait from the same energy-metabolism domain.",
    },
    {
        "anchor_factor": "nonlipid8_F3",
        "anchor_factor_label": "F3_energy_bridge_axis",
        "internal_trait": "Albumin",
        "internal_biomarker_name": "Albumin",
        "study_accession": "GCST90301944",
        "external_trait": "Albumin levels",
        "subset": "external_exploratory_bridge",
        "match_quality": "exploratory_bridge",
        "expected_alignment": "exploratory_bridge",
        "rationale": "Bridge trait relevant to fluid-balance interpretation rather than core factor definition.",
    },
    {
        "anchor_factor": "nonlipid8_F3",
        "anchor_factor_label": "F3_energy_bridge_axis",
        "internal_trait": "Creatinine",
        "internal_biomarker_name": "Creatinine",
        "study_accession": "GCST90301952",
        "external_trait": "Creatinine levels",
        "subset": "external_exploratory_bridge",
        "match_quality": "exploratory_bridge",
        "expected_alignment": "exploratory_bridge",
        "rationale": "Bridge trait relevant to renal/fluid-balance interpretation rather than core factor definition.",
    },
]


def build_catalog():
    df = pd.read_excel(EXCEL_PATH, sheet_name="external support").copy()
    df = df.rename(
        columns={
            "trait": "external_trait",
            "Unnamed: 2": "source_url",
        }
    )
    df = df[["study_accession", "external_trait", "source_url"]].copy()
    df["local_raw_file"] = df["study_accession"].map(lambda x: str(RAW_BASE / f"{x}.tsv").replace("\\", "/"))
    df["munged_sumstats_file"] = df["study_accession"].map(
        lambda x: str(MUNGED_BASE / f"{x}.sumstats.gz").replace("\\", "/")
    )
    df["raw_file_exists"] = df["study_accession"].map(lambda x: (RAW_BASE / f"{x}.tsv").exists())
    df["dataset"] = DATASET
    df["dataset_note"] = DATASET_NOTE
    return df


def factor_rows(module):
    if module == "lipid":
        return [
            {
                "trait": "lipid8_F1",
                "category": "factor",
                "subset": "final_factor",
                "anchor_factor": "lipid8_F1",
                "anchor_factor_label": "F1_TG_rich_axis",
                "internal_trait": "lipid8_F1",
                "internal_biomarker_name": "lipid factor 1",
            },
            {
                "trait": "lipid8_F2",
                "category": "factor",
                "subset": "final_factor",
                "anchor_factor": "lipid8_F2",
                "anchor_factor_label": "F2_HDL_core_axis",
                "internal_trait": "lipid8_F2",
                "internal_biomarker_name": "lipid factor 2",
            },
            {
                "trait": "lipid8_F3",
                "category": "factor",
                "subset": "final_factor",
                "anchor_factor": "lipid8_F3",
                "anchor_factor_label": "F3_CE_structural_axis",
                "internal_trait": "lipid8_F3",
                "internal_biomarker_name": "lipid factor 3",
            },
        ]
    return [
        {
            "trait": "nonlipid8_F1",
            "category": "factor",
            "subset": "final_factor",
            "anchor_factor": "nonlipid8_F1",
            "anchor_factor_label": "F1_ketone_axis",
            "internal_trait": "nonlipid8_F1",
            "internal_biomarker_name": "nonlipid factor 1",
        },
        {
            "trait": "nonlipid8_F2",
            "category": "factor",
            "subset": "final_factor",
            "anchor_factor": "nonlipid8_F2",
            "anchor_factor_label": "F2_amino_acid_axis",
            "internal_trait": "nonlipid8_F2",
            "internal_biomarker_name": "nonlipid factor 2",
        },
        {
            "trait": "nonlipid8_F3",
            "category": "factor",
            "subset": "final_factor",
            "anchor_factor": "nonlipid8_F3",
            "anchor_factor_label": "F3_energy_bridge_axis",
            "internal_trait": "nonlipid8_F3",
            "internal_biomarker_name": "nonlipid factor 3",
        },
    ]


def build_module_tables(module_name, mappings, catalog):
    catalog_map = catalog.set_index("study_accession")
    rows = factor_rows(module_name)
    pair_rows = []

    for item in mappings:
        cat = catalog_map.loc[item["study_accession"]]
        trait_id = f"ext_{item['internal_trait']}_{item['study_accession']}"
        row = {
            "trait": trait_id,
            "category": "external_trait",
            "subset": item["subset"],
            "anchor_factor": item["anchor_factor"],
            "anchor_factor_label": item["anchor_factor_label"],
            "internal_trait": item["internal_trait"],
            "internal_biomarker_name": item["internal_biomarker_name"],
            "external_trait": item["external_trait"],
            "study_accession": item["study_accession"],
            "source_url": cat["source_url"],
            "raw_sumstats_file": cat["local_raw_file"],
            "munged_sumstats_file": cat["munged_sumstats_file"],
            "raw_file_exists": cat["raw_file_exists"],
            "match_quality": item["match_quality"],
            "expected_alignment": item["expected_alignment"],
            "dataset": DATASET,
            "dataset_note": DATASET_NOTE,
            "rationale": item["rationale"],
        }
        rows.append(row)
        pair_rows.append(
            {
                "factor": item["anchor_factor"],
                "trait": trait_id,
                "pair_set": item["subset"],
                "expected_alignment": item["expected_alignment"],
                "internal_trait": item["internal_trait"],
                "internal_biomarker_name": item["internal_biomarker_name"],
                "external_trait": item["external_trait"],
                "study_accession": item["study_accession"],
                "match_quality": item["match_quality"],
                "dataset": DATASET,
                "rationale": item["rationale"],
            }
        )

    traits_df = pd.DataFrame(rows)
    pair_df = pd.DataFrame(pair_rows)
    return traits_df, pair_df


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    MUNGED_BASE.mkdir(parents=True, exist_ok=True)

    catalog = build_catalog().sort_values("study_accession").reset_index(drop=True)
    catalog.to_csv(
        OUTPUT_DIR / "karjalainen233_external_support_catalog.tsv",
        sep="\t",
        index=False,
    )

    lipid_traits, lipid_pairs = build_module_tables("lipid", LIPID_MAPPINGS, catalog)
    lipid_traits.to_csv(
        OUTPUT_DIR / "lipid_final8_external_validation_traits.tsv",
        sep="\t",
        index=False,
    )
    lipid_pairs.to_csv(
        OUTPUT_DIR / "lipid_final8_external_requested_pairs.tsv",
        sep="\t",
        index=False,
    )

    nonlipid_traits, nonlipid_pairs = build_module_tables("nonlipid", NONLIPID_MAPPINGS, catalog)
    nonlipid_traits.to_csv(
        OUTPUT_DIR / "nonlipid_final8_external_validation_traits.tsv",
        sep="\t",
        index=False,
    )
    nonlipid_pairs.to_csv(
        OUTPUT_DIR / "nonlipid_final8_external_requested_pairs.tsv",
        sep="\t",
        index=False,
    )

    summary = pd.DataFrame(
        [
            {
                "table_name": "karjalainen233_external_support_catalog.tsv",
                "n_rows": len(catalog),
                "notes": "All 233 downloaded external support traits from the Excel sheet.",
            },
            {
                "table_name": "lipid_final8_external_validation_traits.tsv",
                "n_rows": len(lipid_traits),
                "notes": "Factors plus exact/approximate matched external support traits for lipid final8.",
            },
            {
                "table_name": "lipid_final8_external_requested_pairs.tsv",
                "n_rows": len(lipid_pairs),
                "notes": "Requested factor-trait external validation pairs for lipid final8.",
            },
            {
                "table_name": "nonlipid_final8_external_validation_traits.tsv",
                "n_rows": len(nonlipid_traits),
                "notes": "Factors plus exact/proxy/support external traits for nonlipid final8.",
            },
            {
                "table_name": "nonlipid_final8_external_requested_pairs.tsv",
                "n_rows": len(nonlipid_pairs),
                "notes": "Requested factor-trait external validation pairs for nonlipid final8.",
            },
        ]
    )
    summary.to_csv(
        OUTPUT_DIR / "external_ldsc_table_manifest.tsv",
        sep="\t",
        index=False,
    )

    print("Wrote external LDSC tables to:", OUTPUT_DIR)


if __name__ == "__main__":
    main()
