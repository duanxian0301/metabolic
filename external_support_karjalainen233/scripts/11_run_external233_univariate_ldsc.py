import argparse
import shlex
import subprocess
from pathlib import Path

import pandas as pd


MANIFEST = Path(r"D:/metabolic/233/ldsc_univariate/external233_ldsc_manifest.tsv")
LDSC_ROOT = Path(r"D:/LDSC/ldsc-master")
HM3 = LDSC_ROOT / "eur_w_ld_chr" / "w_hm3.snplist"
LD = LDSC_ROOT / "eur_w_ld_chr"
WSL_PY2 = "/home/shenjing/miniconda3/envs/ldsc/bin/python"
SUMSTATS_DIR = Path(r"D:/metabolic/233/sumstats")
WORK_DIR = Path(r"D:/metabolic/233/ldsc_univariate")
RESULTS_DIR = WORK_DIR / "results"
LOGS_DIR = WORK_DIR / "logs"
STATUS_DIR = WORK_DIR / "status"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start", type=int, default=1)
    parser.add_argument("--end", type=int, default=None)
    parser.add_argument("--accessions-file", type=str, default=None)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def select_subset(manifest: pd.DataFrame, args):
    if args.accessions_file:
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
        subset["study_accession"] = subset["study_accession"].astype(str)
        subset["__order"] = subset["study_accession"].map({acc: i for i, acc in enumerate(accessions)})
        subset = subset.sort_values("__order").drop(columns="__order").reset_index(drop=True)
        start_idx = max(1, args.start)
        end_idx = len(subset) if args.end is None else min(args.end, len(subset))
        subset = subset.iloc[start_idx - 1 : end_idx].copy()
    else:
        start_idx = max(1, args.start)
        end_idx = len(manifest) if args.end is None else min(args.end, len(manifest))
        subset = manifest.iloc[start_idx - 1 : end_idx].copy()
    return subset, start_idx, end_idx


def run_cmd(cmd, cwd):
    proc = subprocess.run(cmd, cwd=str(cwd), capture_output=True, text=True)
    return proc.returncode, proc.stdout, proc.stderr


def run_wsl_bash(command):
    proc = subprocess.run(
        ["wsl", "-e", "bash", "-lc", command],
        capture_output=True,
    )
    stdout = proc.stdout.decode("utf-8", errors="ignore")
    stderr = proc.stderr.decode("utf-8", errors="ignore")
    return proc.returncode, stdout, stderr


def win_to_wsl(path):
    path = str(path).replace("\\", "/")
    if len(path) >= 2 and path[1] == ":":
        drive = path[0].lower()
        rest = path[2:]
        return "/mnt/{0}{1}".format(drive, rest)
    return path


def ensure_dirs():
    for path in [SUMSTATS_DIR, WORK_DIR, RESULTS_DIR, LOGS_DIR, STATUS_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def run_one(row, overwrite=False):
    accession = row.study_accession
    txt_file = Path(row.output_txt)
    sumstats_file = Path(row.sumstats_file)
    munge_log = Path(row.munge_log)
    h2_prefix = Path(row.h2_prefix)
    h2_log = Path(row.h2_log)

    if not txt_file.exists():
        return {"study_accession": accession, "status": "failed_missing_txt"}

    if sumstats_file.exists() and h2_log.exists() and not overwrite:
        return {"study_accession": accession, "status": "skipped_existing"}

    if overwrite:
        sumstats_file.unlink(missing_ok=True)
        h2_log.unlink(missing_ok=True)
        munge_log.unlink(missing_ok=True)

    munge_cmd = " ".join(
        [
            shlex.quote(WSL_PY2),
            shlex.quote(win_to_wsl(LDSC_ROOT / "munge_sumstats.py")),
            "--sumstats",
            shlex.quote(win_to_wsl(txt_file)),
            "--out",
            shlex.quote(win_to_wsl(str(sumstats_file).replace(".sumstats.gz", ""))),
            "--merge-alleles",
            shlex.quote(win_to_wsl(HM3)),
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
    )
    rc, out, err = run_wsl_bash(munge_cmd)
    munge_log.write_text(out + ("\n[stderr]\n" + err if err else ""), encoding="utf-8")
    if rc != 0 or not sumstats_file.exists():
        return {"study_accession": accession, "status": "failed_munge", "returncode": rc}

    h2_cmd = " ".join(
        [
            shlex.quote(WSL_PY2),
            shlex.quote(win_to_wsl(LDSC_ROOT / "ldsc.py")),
            "--h2",
            shlex.quote(win_to_wsl(sumstats_file)),
            "--ref-ld-chr",
            shlex.quote(win_to_wsl(LD) + "/"),
            "--w-ld-chr",
            shlex.quote(win_to_wsl(LD) + "/"),
            "--out",
            shlex.quote(win_to_wsl(h2_prefix)),
        ]
    )
    rc, out, err = run_wsl_bash(h2_cmd)
    if err and h2_log.exists():
        h2_log.write_text(h2_log.read_text(encoding="utf-8", errors="ignore") + "\n[stderr]\n" + err, encoding="utf-8")
    if rc != 0 or not h2_log.exists():
        return {"study_accession": accession, "status": "failed_h2", "returncode": rc}

    return {"study_accession": accession, "status": "completed"}


def main():
    ensure_dirs()
    manifest = pd.read_csv(MANIFEST, sep="\t")
    args = parse_args()
    subset, start_idx, end_idx = select_subset(manifest, args)

    status_rows = []
    for row in subset.itertuples(index=False):
        try:
            result = run_one(row, overwrite=args.overwrite)
        except Exception as exc:
            result = {"study_accession": row.study_accession, "status": "failed_exception", "error": str(exc)}
        status_rows.append(result)
        print(result, flush=True)

    status_df = pd.DataFrame(status_rows)
    status_path = STATUS_DIR / f"external233_ldsc_status_{start_idx}_{end_idx}.tsv"
    status_df.to_csv(status_path, sep="\t", index=False)
    print(f"Wrote status file to {status_path}")


if __name__ == "__main__":
    main()
