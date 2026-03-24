import pandas as pd
from pathlib import Path


PRIORITY_PATH = Path(r"D:/metabolic/233/ldsc_univariate/external_validation_priority/lipid_nonlipid_external_validation_priority.tsv")
OUT_DIR = Path(r"D:/codex/metabolic_repo/external_support_karjalainen233/inputs/external_bivariate_round2")


LIPID_ACCESSIONS = [
    "GCST90302081",  # HDL CE small HDL
    "GCST90302021",  # VLDL CE large VLDL
    "GCST90302140",  # VLDL CE very large VLDL
]

NONLIPID_ACCESSIONS = [
    "GCST90301943",  # Alanine
    "GCST90301987",  # Isoleucine
    "GCST90301949",  # Citrate
]

LIPID_FACTORS = ["lipid8_F1", "lipid8_F2", "lipid8_F3"]
NONLIPID_FACTORS = ["nonlipid8_F1", "nonlipid8_F2", "nonlipid8_F3"]


def build_module(module_name, accessions, factors):
    dt = pd.read_csv(PRIORITY_PATH, sep="\t")
    dt = dt[(dt["module"] == module_name) & (dt["study_accession"].isin(accessions))].copy()
    dt["sumstats_file"] = dt["study_accession"].map(lambda x: f"D:/metabolic/233/sumstats/{x}.sumstats.gz")
    dt["target_trait"] = dt["study_accession"].map(lambda x: f"ext_{x}")
    dt["category"] = "external_trait"

    trait_manifest = dt[
        [
            "target_trait",
            "study_accession",
            "external_trait",
            "internal_trait",
            "anchor_factor",
            "anchor_factor_label",
            "subset",
            "match_quality",
            "validation_priority",
            "sample_total_n",
            "discovery_sample_ancestry",
            "h2",
            "h2_se",
            "h2_z",
            "intercept",
            "ratio",
            "lambda_gc",
            "mean_chi2",
            "sumstats_file",
        ]
    ].rename(
        columns={
            "external_trait": "biomarker_name",
            "sample_total_n": "N",
            "target_trait": "trait",
        }
    )

    pair_rows = []
    for _, row in dt.iterrows():
        for factor in factors:
            pair_rows.append(
                {
                    "factor": factor,
                    "target_trait": row["target_trait"],
                    "anchor_factor": row["anchor_factor"],
                    "internal_trait": row["internal_trait"],
                    "external_trait": row["external_trait"],
                    "study_accession": row["study_accession"],
                    "validation_priority": row["validation_priority"],
                    "pair_set": "support_match" if factor == row["anchor_factor"] else "cross_factor",
                    "expected_alignment": "support" if factor == row["anchor_factor"] else "nonprimary",
                }
            )
    pairs = pd.DataFrame(pair_rows)

    module_dir = OUT_DIR / module_name
    module_dir.mkdir(parents=True, exist_ok=True)
    trait_manifest.to_csv(module_dir / "external_traits_manifest.tsv", sep="\t", index=False)
    pairs.to_csv(module_dir / "requested_pairs.tsv", sep="\t", index=False)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    build_module("lipid", LIPID_ACCESSIONS, LIPID_FACTORS)
    build_module("nonlipid", NONLIPID_ACCESSIONS, NONLIPID_FACTORS)
    print(OUT_DIR)


if __name__ == "__main__":
    main()
