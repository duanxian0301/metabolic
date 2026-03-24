import pandas as pd
from pathlib import Path


EXCEL_PATH = Path(r"D:/metabolic/Metabolic-data.xlsx")
EXCEL_SHEET = "external support"
PROCESSING_MANIFEST = Path(r"D:/metabolic/233/manifests/external233_processing_manifest.tsv")
GWAS_DIR = Path(r"D:/metabolic/233/GWAS")
QC_DIR = Path(r"D:/metabolic/233/QC")
SUMSTATS_DIR = Path(r"D:/metabolic/233/sumstats")
LDSC_DIR = Path(r"D:/metabolic/233/ldsc_univariate")
OUT_MANIFEST = LDSC_DIR / "external233_ldsc_manifest.tsv"


def main():
    LDSC_DIR.mkdir(parents=True, exist_ok=True)
    SUMSTATS_DIR.mkdir(parents=True, exist_ok=True)

    ext = pd.read_excel(EXCEL_PATH, sheet_name=EXCEL_SHEET)
    ext = ext.rename(columns={"Unnamed: 2": "raw_url"})
    ext["study_accession"] = ext["study_accession"].astype(str)

    proc = pd.read_csv(PROCESSING_MANIFEST, sep="\t")
    proc["study_accession"] = proc["study_accession"].astype(str)

    merged = ext.merge(
        proc[
            [
                "study_accession",
                "reported_trait",
                "raw_url",
                "sample_total_n",
                "discovery_sample_ancestry",
                "raw_file",
                "yaml_file",
                "output_txt",
                "qc_summary_file",
                "raw_position_build",
                "harmonised_url",
            ]
        ],
        on=["study_accession", "raw_url"],
        how="left",
        validate="one_to_one",
    )

    merged["txt_exists"] = merged["output_txt"].map(lambda x: Path(x).exists())
    merged["qc_exists"] = merged["qc_summary_file"].map(lambda x: Path(x).exists())
    merged["sumstats_file"] = merged["study_accession"].map(lambda x: str(SUMSTATS_DIR / f"{x}.sumstats.gz"))
    merged["munge_log"] = merged["study_accession"].map(lambda x: str(LDSC_DIR / "logs" / f"{x}.munge.log"))
    merged["h2_prefix"] = merged["study_accession"].map(lambda x: str(LDSC_DIR / "results" / f"{x}_h2"))
    merged["h2_log"] = merged["study_accession"].map(lambda x: str(LDSC_DIR / "results" / f"{x}_h2.log"))

    ordered = [
        "study_accession",
        "trait",
        "reported_trait",
        "raw_url",
        "sample_total_n",
        "discovery_sample_ancestry",
        "raw_position_build",
        "raw_file",
        "yaml_file",
        "output_txt",
        "qc_summary_file",
        "txt_exists",
        "qc_exists",
        "harmonised_url",
        "sumstats_file",
        "munge_log",
        "h2_prefix",
        "h2_log",
    ]
    merged = merged[ordered]
    merged.to_csv(OUT_MANIFEST, sep="\t", index=False)
    print(f"Wrote manifest with {len(merged)} traits to {OUT_MANIFEST}")


if __name__ == "__main__":
    main()
