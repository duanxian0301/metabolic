from pathlib import Path
import json
import re
import urllib.request

import pandas as pd
import yaml


EXCEL_PATH = Path(r"D:/metabolic/Metabolic-data.xlsx")
OUTPUT_DIR = Path(r"D:/metabolic/233/manifests")
RAW_BASE = Path(r"D:/metabolic")
YAML_BASE = Path(r"D:/metabolic/233/yaml")
GWAS_OUT_DIR = Path(r"D:/metabolic/233/GWAS")
QC_OUT_DIR = Path(r"D:/metabolic/233/QC")
PMID = "38448586"
API_URL = f"https://www.ebi.ac.uk/gwas/api/v2/publications/{PMID}/studies?fullPvalueSet=true&size=500"


def parse_total_n(sample_list):
    total = 0
    for item in sample_list or []:
        m = re.match(r"^\s*([\d,]+)\s+", str(item))
        if m:
            total += int(m.group(1).replace(",", ""))
    return total if total > 0 else None


def parse_total_n_from_yaml(yaml_path):
    with open(yaml_path, "r", encoding="utf-8") as fh:
        meta = yaml.safe_load(fh)
    samples = meta.get("samples", []) or []
    total = 0
    for sample in samples:
        n = sample.get("sample_size")
        if n is not None:
            total += int(n)
    assembly = meta.get("genome_assembly")
    return total if total > 0 else None, assembly


def harmonised_url_from_raw(raw_url, accession):
    raw_url = str(raw_url).strip()
    if not raw_url:
        return ""
    base = raw_url.rsplit("/", 1)[0]
    return f"{base}/harmonised/{accession}.h.tsv.gz"


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    GWAS_OUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_OUT_DIR.mkdir(parents=True, exist_ok=True)

    ext_df = pd.read_excel(EXCEL_PATH, sheet_name="external support").rename(
        columns={"trait": "reported_trait", "Unnamed: 2": "raw_url"}
    )
    ext_df = ext_df[["study_accession", "reported_trait", "raw_url"]].copy()

    with urllib.request.urlopen(API_URL) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    api_df = pd.DataFrame(payload["_embedded"]["studies"])[
        ["accessionId", "reportedTrait", "initialSampleDescription", "discoverySampleAncestry"]
    ].rename(
        columns={
            "accessionId": "study_accession",
            "reportedTrait": "api_reported_trait",
            "initialSampleDescription": "initial_sample_description",
            "discoverySampleAncestry": "discovery_sample_ancestry",
        }
    )
    api_df["sample_total_n"] = api_df["initial_sample_description"].apply(parse_total_n)

    df = ext_df.merge(api_df, on="study_accession", how="left", validate="one_to_one")
    df["raw_file"] = df["study_accession"].map(lambda x: str(RAW_BASE / f"{x}.tsv").replace("\\", "/"))
    df["yaml_file"] = df["study_accession"].map(lambda x: str(YAML_BASE / f"{x}.tsv-meta.yaml").replace("\\", "/"))
    df["raw_file_exists"] = df["study_accession"].map(lambda x: (RAW_BASE / f"{x}.tsv").exists())
    df["yaml_file_exists"] = df["study_accession"].map(lambda x: (YAML_BASE / f"{x}.tsv-meta.yaml").exists())
    df["harmonised_url"] = df.apply(lambda r: harmonised_url_from_raw(r["raw_url"], r["study_accession"]), axis=1)
    df["output_txt"] = df["study_accession"].map(lambda x: str(GWAS_OUT_DIR / f"{x}.txt").replace("\\", "/"))
    df["qc_summary_file"] = df["study_accession"].map(lambda x: str(QC_OUT_DIR / f"{x}.qc_summary.tsv").replace("\\", "/"))
    parsed = df["study_accession"].map(
        lambda x: parse_total_n_from_yaml(YAML_BASE / f"{x}.tsv-meta.yaml")
    )
    df["sample_total_n_yaml"] = [x[0] for x in parsed]
    df["raw_position_build"] = [x[1] for x in parsed]
    df["sample_total_n"] = df["sample_total_n_yaml"].fillna(df["sample_total_n"])
    df["build_decision"] = "Use raw chromosome/base_pair_location as GRCh37 coordinates; do not liftover."
    df["rsid_source"] = "GWAS Catalog harmonised file"
    df["standard_output_columns"] = "SNP CHR BP A1 A2 BETA SE P N"
    df["p_floor"] = 1e-300
    df["notes"] = (
        "N and assembly come from local YAML metadata; rsID is taken from the harmonised GWAS Catalog file; "
        "coordinates and effect estimates come from the raw GWAS file."
    )

    df = df.sort_values("study_accession").reset_index(drop=True)

    df.to_csv(OUTPUT_DIR / "external233_processing_manifest.tsv", sep="\t", index=False)

    build_note = pd.DataFrame(
        [
            {
                "field": "publication",
                "value": "Karjalainen MK et al. Nature 2024 (PMID 38448586)",
            },
            {
                "field": "raw_position_build",
                "value": "GRCh37",
            },
            {
                "field": "harmonised_usage",
                "value": "rsID lookup only; harmonised coordinates are not used for final output.",
            },
            {
                "field": "qc_rules",
                "value": "autosomes 1-22, valid rsID, A/C/G/T alleles, SE>0, 0<P<=1 with floor 1e-300, duplicate SNPs removed by lowest P.",
            },
        ]
    )
    build_note.to_csv(OUTPUT_DIR / "external233_build_and_qc_notes.tsv", sep="\t", index=False)
    print("Wrote manifest to", OUTPUT_DIR / "external233_processing_manifest.tsv")


if __name__ == "__main__":
    main()
