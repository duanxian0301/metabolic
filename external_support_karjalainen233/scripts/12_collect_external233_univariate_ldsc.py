import re
from pathlib import Path

import pandas as pd


EXCEL_PATH = Path(r"D:/metabolic/Metabolic-data.xlsx")
EXCEL_SHEET = "external support"
PROCESSING_MANIFEST = Path(r"D:/metabolic/233/manifests/external233_processing_manifest.tsv")
LDSC_MANIFEST = Path(r"D:/metabolic/233/ldsc_univariate/external233_ldsc_manifest.tsv")
QC_DIR = Path(r"D:/metabolic/233/QC")
RESULTS_DIR = Path(r"D:/metabolic/233/ldsc_univariate/results")
OUT_DIR = Path(r"D:/metabolic/233/ldsc_univariate")

PATTERNS = {
    "n_sumstats_read": re.compile(r"Read summary statistics for ([\d,]+) SNPs\."),
    "n_ref_ld_read": re.compile(r"Read reference panel LD Scores for ([\d,]+) SNPs\."),
    "n_after_ref_merge": re.compile(r"After merging with reference panel LD, ([\d,]+) SNPs remain\."),
    "n_after_weight_merge": re.compile(r"After merging with regression SNP LD, ([\d,]+) SNPs remain\."),
    "h2": re.compile(r"Total Observed scale h2:\s+([-\d\.eE]+)\s+\(([-\d\.eE]+)\)"),
    "lambda_gc": re.compile(r"Lambda GC:\s+([-\d\.eE]+)"),
    "mean_chi2": re.compile(r"Mean Chi\^2:\s+([-\d\.eE]+)"),
    "intercept": re.compile(r"Intercept:\s+([-\d\.eE]+)\s+\(([-\d\.eE]+)\)"),
    "ratio": re.compile(r"Ratio:\s+([-\d\.eE]+)\s+\(([-\d\.eE]+)\)"),
    "ratio_lt0": re.compile(r"Ratio < 0"),
    "analysis_finished": re.compile(r"Analysis finished at (.+)"),
    "time_elapsed": re.compile(r"Total time elapsed:\s+([-\d\.eE]+)s"),
}


def parse_qc(qc_path: Path):
    if not qc_path.exists():
        return {}
    dt = pd.read_csv(qc_path, sep="\t")
    return dict(zip(dt["metric"], dt["value"]))


def parse_log(log_path: Path):
    out = {"ldsc_log_exists": log_path.exists(), "ldsc_status": "missing_log"}
    if not log_path.exists():
        return out
    txt = log_path.read_text(encoding="utf-8", errors="ignore")
    out["ldsc_status"] = "completed" if "Total Observed scale h2" in txt else "incomplete_log"
    for key, pattern in PATTERNS.items():
        if key == "h2":
            m = pattern.search(txt)
            if m:
                out["h2"] = float(m.group(1))
                out["h2_se"] = float(m.group(2))
        elif key == "intercept":
            m = pattern.search(txt)
            if m:
                out["intercept"] = float(m.group(1))
                out["intercept_se"] = float(m.group(2))
        elif key == "ratio":
            m = pattern.search(txt)
            if m:
                out["ratio"] = float(m.group(1))
                out["ratio_se"] = float(m.group(2))
        elif key == "ratio_lt0":
            out["ratio_note"] = "Ratio < 0 (usually indicates GC correction)." if pattern.search(txt) else None
        else:
            m = pattern.search(txt)
            if m:
                val = m.group(1)
                if key.startswith("n_"):
                    out[key] = int(val.replace(",", ""))
                elif key in {"lambda_gc", "mean_chi2", "time_elapsed"}:
                    out[key] = float(val)
                else:
                    out[key] = val
    if "h2" in out and "h2_se" in out and out["h2_se"] not in [0, None]:
        out["h2_z"] = out["h2"] / out["h2_se"]
    return out


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    ext = pd.read_excel(EXCEL_PATH, sheet_name=EXCEL_SHEET).rename(columns={"Unnamed: 2": "raw_url"})
    proc = pd.read_csv(PROCESSING_MANIFEST, sep="\t")
    ldsc_manifest = pd.read_csv(LDSC_MANIFEST, sep="\t")

    rows = []
    for row in ldsc_manifest.itertuples(index=False):
        qc = parse_qc(Path(row.qc_summary_file))
        ldsc = parse_log(Path(row.h2_log))
        rows.append(
            {
                "study_accession": row.study_accession,
                "trait": row.trait,
                "reported_trait": row.reported_trait,
                "raw_url": row.raw_url,
                "sample_total_n": row.sample_total_n,
                "discovery_sample_ancestry": row.discovery_sample_ancestry,
                "raw_position_build": row.raw_position_build,
                "output_txt": row.output_txt,
                "sumstats_file": row.sumstats_file,
                "h2_log": row.h2_log,
                **qc,
                **ldsc,
            }
        )

    out = pd.DataFrame(rows)
    out = ext.merge(out, on=["study_accession", "trait", "raw_url"], how="left", validate="one_to_one")
    out = out.merge(
        proc[["study_accession", "harmonised_url", "yaml_file", "raw_file"]],
        on="study_accession",
        how="left",
        validate="one_to_one",
    )

    ordered_front = [
        "study_accession",
        "trait",
        "reported_trait",
        "raw_url",
        "harmonised_url",
        "raw_file",
        "yaml_file",
        "output_txt",
        "sumstats_file",
        "h2_log",
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
        "n_ref_ld_read",
        "n_after_ref_merge",
        "n_after_weight_merge",
        "h2",
        "h2_se",
        "h2_z",
        "intercept",
        "intercept_se",
        "ratio",
        "ratio_se",
        "ratio_note",
        "lambda_gc",
        "mean_chi2",
        "analysis_finished",
        "time_elapsed",
    ]
    ordered_existing = [c for c in ordered_front if c in out.columns]
    remaining = [c for c in out.columns if c not in ordered_existing]
    out = out[ordered_existing + remaining]

    out_path = OUT_DIR / "external233_univariate_ldsc_merged.tsv"
    out.to_csv(out_path, sep="\t", index=False)

    short = out[
        [
            "study_accession",
            "trait",
            "sample_total_n",
            "discovery_sample_ancestry",
            "raw_position_build",
            "status",
            "ldsc_status",
            "n_final",
            "n_sumstats_read",
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
        ]
    ].copy()
    short.to_csv(OUT_DIR / "external233_univariate_ldsc_summary.tsv", sep="\t", index=False)
    print(f"Wrote merged LDSC tables to {OUT_DIR}")


if __name__ == "__main__":
    main()
