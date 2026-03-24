#!/usr/bin/env python3
import argparse
import csv
import subprocess
from pathlib import Path


MANIFEST = Path("/mnt/d/metabolic/233/ldsc_univariate/external233_ldsc_manifest.tsv")
LDSC_ROOT = Path("/mnt/d/LDSC/ldsc-master")
HM3 = LDSC_ROOT / "eur_w_ld_chr" / "w_hm3.snplist"
LD = LDSC_ROOT / "eur_w_ld_chr"
PY2 = "/home/shenjing/miniconda3/envs/ldsc/bin/python"
STATUS_DIR = Path("/mnt/d/metabolic/233/ldsc_univariate/status")


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--end", type=int, default=None)
    parser.add_argument("--accessions-file", type=str, default=None)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def read_tsv(path):
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def win_to_wsl(path_str):
    path = str(path_str).replace("\\", "/")
    if len(path) >= 2 and path[1] == ":":
        return "/mnt/{0}{1}".format(path[0].lower(), path[2:])
    return path


def select_subset(manifest_rows, args):
    if args.accessions_file:
        acc_rows = read_tsv(args.accessions_file)
        accessions = [str(r.get("study_accession", "")).strip() for r in acc_rows if str(r.get("study_accession", "")).strip()]
        order = {acc: idx for idx, acc in enumerate(accessions)}
        subset = [r for r in manifest_rows if str(r["study_accession"]) in order]
        subset.sort(key=lambda r: order[str(r["study_accession"])])
        start_idx = max(1, args.start)
        end_idx = len(subset) if args.end is None else min(args.end, len(subset))
        subset = subset[start_idx - 1 : end_idx]
    else:
        start_idx = max(1, args.start)
        end_idx = len(manifest_rows) if args.end is None else min(args.end, len(manifest_rows))
        subset = manifest_rows[start_idx - 1 : end_idx]
    return subset, start_idx, end_idx


def run_cmd(cmd):
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return proc.returncode, proc.stdout, proc.stderr


def run_one(row, overwrite=False):
    accession = row["study_accession"]
    txt_file = Path(win_to_wsl(row["output_txt"]))
    sumstats_file = Path(win_to_wsl(row["sumstats_file"]))
    munge_log = Path(win_to_wsl(row["munge_log"]))
    h2_prefix = Path(win_to_wsl(row["h2_prefix"]))
    h2_log = Path(win_to_wsl(row["h2_log"]))

    if not txt_file.exists():
        return {"study_accession": accession, "status": "failed_missing_txt"}

    if sumstats_file.exists() and h2_log.exists() and not overwrite:
        return {"study_accession": accession, "status": "skipped_existing"}

    if overwrite:
        for path in [sumstats_file, h2_log, munge_log]:
            if path.exists():
                path.unlink()

    munge_cmd = [
        PY2,
        str(LDSC_ROOT / "munge_sumstats.py"),
        "--sumstats",
        str(txt_file),
        "--out",
        str(sumstats_file).replace(".sumstats.gz", ""),
        "--merge-alleles",
        str(HM3),
        "--chunksize",
        "500000",
        "--a1",
        "A1",
        "--a2",
        "A2",
        "--snp",
        "SNP",
        "--N-col",
        "N",
        "--signed-sumstats",
        "BETA,0",
        "--p",
        "P",
    ]
    rc, out, err = run_cmd(munge_cmd)
    munge_log.parent.mkdir(parents=True, exist_ok=True)
    munge_log.write_text(out + ("\n[stderr]\n" + err if err else ""), encoding="utf-8")
    if rc != 0 or not sumstats_file.exists():
        return {"study_accession": accession, "status": "failed_munge", "returncode": rc}

    h2_cmd = [
        PY2,
        str(LDSC_ROOT / "ldsc.py"),
        "--h2",
        str(sumstats_file),
        "--ref-ld-chr",
        str(LD) + "/",
        "--w-ld-chr",
        str(LD) + "/",
        "--out",
        str(h2_prefix),
    ]
    rc, out, err = run_cmd(h2_cmd)
    if err and h2_log.exists():
        h2_log.write_text(h2_log.read_text(encoding="utf-8", errors="ignore") + "\n[stderr]\n" + err, encoding="utf-8")
    if rc != 0 or not h2_log.exists():
        return {"study_accession": accession, "status": "failed_h2", "returncode": rc}

    return {"study_accession": accession, "status": "completed"}


def main():
    args = parse_args()
    manifest_rows = read_tsv(MANIFEST)
    subset, start_idx, end_idx = select_subset(manifest_rows, args)
    status_rows = []
    for row in subset:
        try:
            result = run_one(row, overwrite=args.overwrite)
        except Exception as exc:
            result = {"study_accession": row["study_accession"], "status": "failed_exception", "error": str(exc)}
        status_rows.append(result)
        print(result, flush=True)

    status_path = STATUS_DIR / "external233_ldsc_wsl_status_{0}_{1}.tsv".format(start_idx, end_idx)
    fieldnames = sorted({key for row in status_rows for key in row.keys()})
    write_tsv(status_path, status_rows, fieldnames)
    print("Wrote status file to {0}".format(status_path))


if __name__ == "__main__":
    main()
