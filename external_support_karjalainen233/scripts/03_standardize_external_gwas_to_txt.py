import argparse
import csv
import gzip
import io
import math
import os
import time
import urllib.request
from pathlib import Path

import pandas as pd


MANIFEST = Path(r"D:/metabolic/233/manifests/external233_processing_manifest.tsv")
GWAS_OUT_DIR = Path(r"D:/metabolic/233/GWAS")
QC_OUT_DIR = Path(r"D:/metabolic/233/QC")
VALID_ALLELES = {"A", "C", "G", "T"}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=1, help="1-based row start in manifest")
    parser.add_argument("--end", type=int, default=None, help="1-based row end in manifest")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--max-files", type=int, default=None)
    parser.add_argument("--accessions-file", type=str, default=None, help="Optional TSV/TXT with study_accession column or one accession per line")
    return parser.parse_args()


def open_harmonised_reader(url):
    last_error = None
    for _ in range(3):
        try:
            response = urllib.request.urlopen(url, timeout=180)
            gzip_stream = gzip.GzipFile(fileobj=response)
            text_stream = io.TextIOWrapper(gzip_stream, encoding="utf-8")
            reader = csv.DictReader(text_stream, delimiter="\t")
            return response, gzip_stream, text_stream, reader
        except Exception as exc:
            last_error = exc
            time.sleep(3)
    raise last_error


def load_harmonised_chrom(target_chrom, h_reader, h_pending_row):
    chrom_map = {}
    pending = h_pending_row

    if pending is None:
        try:
            pending = next(h_reader)
        except StopIteration:
            return chrom_map, None

    while pending is not None:
        chrom = coerce_int(pending.get("chromosome"))
        if chrom is None:
            try:
                pending = next(h_reader)
                continue
            except StopIteration:
                pending = None
                break
        if chrom < target_chrom:
            try:
                pending = next(h_reader)
                continue
            except StopIteration:
                pending = None
                break
        break

    while pending is not None:
        chrom = coerce_int(pending.get("chromosome"))
        if chrom != target_chrom:
            break
        variant_id = pending.get("variant_id")
        if variant_id is not None and variant_id not in chrom_map:
            chrom_map[variant_id] = pending.get("rsid", "")
        try:
            pending = next(h_reader)
        except StopIteration:
            pending = None
            break

    return chrom_map, pending


def is_valid_rsid(value):
    value = str(value).strip()
    return value.startswith("rs") and len(value) > 2


def coerce_int(value):
    try:
        return int(float(value))
    except Exception:
        return None


def coerce_float(value):
    try:
        x = float(value)
        if math.isnan(x):
            return None
        return x
    except Exception:
        return None


def standardize_one(row, overwrite=False):
    accession = row.study_accession
    raw_file = Path(row.raw_file)
    out_file = Path(row.output_txt)
    qc_file = Path(row.qc_summary_file)
    tmp_file = out_file.with_suffix(".tmp.txt")
    p_floor = float(row.p_floor)
    n_value = int(row.sample_total_n)

    if out_file.exists() and qc_file.exists() and not overwrite:
        return {"study_accession": accession, "status": "skipped_existing"}

    tmp_file.unlink(missing_ok=True)

    GWAS_OUT_DIR.mkdir(parents=True, exist_ok=True)
    QC_OUT_DIR.mkdir(parents=True, exist_ok=True)

    counters = {
        "n_raw_rows": 0,
        "n_variant_id_mismatch": 0,
        "n_missing_basic": 0,
        "n_invalid_rsid": 0,
        "n_non_autosome": 0,
        "n_invalid_bp": 0,
        "n_invalid_alleles": 0,
        "n_invalid_beta": 0,
        "n_invalid_se": 0,
        "n_invalid_p": 0,
        "n_p_capped": 0,
        "n_retained_pre_dedup": 0,
    }

    response, gzip_stream, text_stream, h_reader = open_harmonised_reader(row.harmonised_url)
    h_pending_row = None
    h_current_chrom = None
    h_current_map = {}
    try:
        with raw_file.open("r", encoding="utf-8", newline="") as raw_handle, tmp_file.open(
            "w", encoding="utf-8", newline=""
        ) as out_handle:
            raw_reader = csv.DictReader(raw_handle, delimiter="\t")
            writer = csv.writer(out_handle, delimiter="\t")
            writer.writerow(["SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "N"])

            for raw_row in raw_reader:
                counters["n_raw_rows"] += 1

                raw_variant = raw_row.get("variant_id")
                chrom = coerce_int(raw_row.get("chromosome"))
                if chrom is None:
                    counters["n_missing_basic"] += 1
                    continue

                if h_current_chrom != chrom:
                    h_current_map, h_pending_row = load_harmonised_chrom(chrom, h_reader, h_pending_row)
                    h_current_chrom = chrom

                rsid = h_current_map.get(raw_variant)
                if rsid is None:
                    counters["n_variant_id_mismatch"] += 1
                    continue

                rsid = rsid.strip()
                if not is_valid_rsid(rsid):
                    counters["n_invalid_rsid"] += 1
                    continue

                bp = coerce_int(raw_row.get("base_pair_location"))
                a1 = str(raw_row.get("effect_allele", "")).upper().strip()
                a2 = str(raw_row.get("other_allele", "")).upper().strip()
                beta = coerce_float(raw_row.get("beta"))
                se = coerce_float(raw_row.get("standard_error"))
                p = coerce_float(raw_row.get("p_value"))

                if chrom is None or bp is None or a1 == "" or a2 == "":
                    counters["n_missing_basic"] += 1
                    continue
                if chrom < 1 or chrom > 22:
                    counters["n_non_autosome"] += 1
                    continue
                if bp <= 0:
                    counters["n_invalid_bp"] += 1
                    continue
                if a1 not in VALID_ALLELES or a2 not in VALID_ALLELES or a1 == a2:
                    counters["n_invalid_alleles"] += 1
                    continue
                if beta is None:
                    counters["n_invalid_beta"] += 1
                    continue
                if se is None or se <= 0:
                    counters["n_invalid_se"] += 1
                    continue
                if p is None or p <= 0 or p > 1:
                    counters["n_invalid_p"] += 1
                    continue
                if p < p_floor:
                    p = p_floor
                    counters["n_p_capped"] += 1

                writer.writerow(
                    [
                        rsid,
                        chrom,
                        bp,
                        a1,
                        a2,
                        format(beta, ".10g"),
                        format(se, ".10g"),
                        format(p, ".10g"),
                        n_value,
                    ]
                )
                counters["n_retained_pre_dedup"] += 1
    finally:
        text_stream.close()
        gzip_stream.close()
        response.close()

    if not tmp_file.exists() or os.path.getsize(tmp_file) == 0:
        qc_dt = pd.DataFrame(
            [
                {"metric": "study_accession", "value": accession},
                {"metric": "reported_trait", "value": row.reported_trait},
                {"metric": "sample_total_n", "value": n_value},
                {"metric": "raw_position_build", "value": row.raw_position_build},
                {"metric": "n_raw_rows", "value": counters["n_raw_rows"]},
                {"metric": "n_variant_id_mismatch", "value": counters["n_variant_id_mismatch"]},
                {"metric": "n_missing_basic", "value": counters["n_missing_basic"]},
                {"metric": "n_invalid_rsid", "value": counters["n_invalid_rsid"]},
                {"metric": "n_non_autosome", "value": counters["n_non_autosome"]},
                {"metric": "n_invalid_bp", "value": counters["n_invalid_bp"]},
                {"metric": "n_invalid_alleles", "value": counters["n_invalid_alleles"]},
                {"metric": "n_invalid_beta", "value": counters["n_invalid_beta"]},
                {"metric": "n_invalid_se", "value": counters["n_invalid_se"]},
                {"metric": "n_invalid_p", "value": counters["n_invalid_p"]},
                {"metric": "n_p_capped", "value": counters["n_p_capped"]},
                {"metric": "n_retained_pre_dedup", "value": 0},
                {"metric": "n_duplicates_removed", "value": 0},
                {"metric": "n_final", "value": 0},
                {"metric": "status", "value": "failed_zero_rows"},
            ]
        )
        qc_dt.to_csv(qc_file, sep="\t", index=False)
        tmp_file.unlink(missing_ok=True)
        return {"study_accession": accession, "status": "failed_zero_rows", "n_final": 0}

    dedup_df = pd.read_csv(tmp_file, sep="\t")
    if dedup_df.empty:
        qc_dt = pd.DataFrame(
            [
                {"metric": "study_accession", "value": accession},
                {"metric": "reported_trait", "value": row.reported_trait},
                {"metric": "sample_total_n", "value": n_value},
                {"metric": "raw_position_build", "value": row.raw_position_build},
                {"metric": "n_raw_rows", "value": counters["n_raw_rows"]},
                {"metric": "n_variant_id_mismatch", "value": counters["n_variant_id_mismatch"]},
                {"metric": "n_missing_basic", "value": counters["n_missing_basic"]},
                {"metric": "n_invalid_rsid", "value": counters["n_invalid_rsid"]},
                {"metric": "n_non_autosome", "value": counters["n_non_autosome"]},
                {"metric": "n_invalid_bp", "value": counters["n_invalid_bp"]},
                {"metric": "n_invalid_alleles", "value": counters["n_invalid_alleles"]},
                {"metric": "n_invalid_beta", "value": counters["n_invalid_beta"]},
                {"metric": "n_invalid_se", "value": counters["n_invalid_se"]},
                {"metric": "n_invalid_p", "value": counters["n_invalid_p"]},
                {"metric": "n_p_capped", "value": counters["n_p_capped"]},
                {"metric": "n_retained_pre_dedup", "value": 0},
                {"metric": "n_duplicates_removed", "value": 0},
                {"metric": "n_final", "value": 0},
                {"metric": "status", "value": "failed_zero_rows"},
            ]
        )
        qc_dt.to_csv(qc_file, sep="\t", index=False)
        tmp_file.unlink(missing_ok=True)
        return {"study_accession": accession, "status": "failed_zero_rows", "n_final": 0}

    n_pre_dedup = len(dedup_df)
    dedup_df = dedup_df.sort_values(["P", "CHR", "BP", "SNP"], ascending=[True, True, True, True])
    dedup_df = dedup_df.drop_duplicates(subset=["SNP"], keep="first")
    n_final = len(dedup_df)
    n_duplicates_removed = n_pre_dedup - n_final
    dedup_df.to_csv(out_file, sep="\t", index=False)
    tmp_file.unlink(missing_ok=True)

    qc_dt = pd.DataFrame(
        [
            {"metric": "study_accession", "value": accession},
            {"metric": "reported_trait", "value": row.reported_trait},
            {"metric": "sample_total_n", "value": n_value},
            {"metric": "raw_position_build", "value": row.raw_position_build},
            {"metric": "n_raw_rows", "value": counters["n_raw_rows"]},
            {"metric": "n_variant_id_mismatch", "value": counters["n_variant_id_mismatch"]},
            {"metric": "n_missing_basic", "value": counters["n_missing_basic"]},
            {"metric": "n_invalid_rsid", "value": counters["n_invalid_rsid"]},
            {"metric": "n_non_autosome", "value": counters["n_non_autosome"]},
            {"metric": "n_invalid_bp", "value": counters["n_invalid_bp"]},
            {"metric": "n_invalid_alleles", "value": counters["n_invalid_alleles"]},
            {"metric": "n_invalid_beta", "value": counters["n_invalid_beta"]},
            {"metric": "n_invalid_se", "value": counters["n_invalid_se"]},
            {"metric": "n_invalid_p", "value": counters["n_invalid_p"]},
            {"metric": "n_p_capped", "value": counters["n_p_capped"]},
            {"metric": "n_retained_pre_dedup", "value": n_pre_dedup},
            {"metric": "n_duplicates_removed", "value": n_duplicates_removed},
            {"metric": "n_final", "value": n_final},
            {"metric": "status", "value": "completed"},
        ]
    )
    qc_dt.to_csv(qc_file, sep="\t", index=False)

    return {
        "study_accession": accession,
        "status": "completed",
        "n_final": n_final,
        "n_duplicates_removed": n_duplicates_removed,
        "n_p_capped": counters["n_p_capped"],
    }


def main():
    args = parse_args()
    manifest = pd.read_csv(MANIFEST, sep="\t")
    if args.accessions_file is not None:
        acc_path = Path(args.accessions_file)
        if acc_path.suffix.lower() == ".tsv":
            acc_df = pd.read_csv(acc_path, sep="\t")
            if "study_accession" in acc_df.columns:
                accessions = acc_df["study_accession"].astype(str).tolist()
            else:
                accessions = acc_df.iloc[:, 0].astype(str).tolist()
        else:
            accessions = [line.strip() for line in acc_path.read_text(encoding="utf-8").splitlines() if line.strip()]
        subset = manifest[manifest["study_accession"].astype(str).isin(accessions)].copy()
        start_idx = 1
        end_idx = len(subset)
    else:
        start_idx = max(1, args.start)
        end_idx = len(manifest) if args.end is None else min(args.end, len(manifest))
        subset = manifest.iloc[start_idx - 1 : end_idx].copy()
    if args.max_files is not None:
        subset = subset.head(args.max_files)

    status_rows = []
    for row in subset.itertuples(index=False):
        try:
            result = standardize_one(row, overwrite=args.overwrite)
        except Exception as exc:
            result = {"study_accession": row.study_accession, "status": "failed_exception", "error": str(exc)}
        status_rows.append(result)
        print(result, flush=True)

    status_df = pd.DataFrame(status_rows)
    status_name = f"processing_status_{start_idx}_{end_idx}.tsv"
    status_df.to_csv(QC_OUT_DIR / status_name, sep="\t", index=False)
    print("Wrote status file to", QC_OUT_DIR / status_name)


if __name__ == "__main__":
    main()
