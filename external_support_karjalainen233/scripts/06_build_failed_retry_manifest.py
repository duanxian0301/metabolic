from pathlib import Path

import pandas as pd


MANIFEST = Path(r"D:/metabolic/233/manifests/external233_processing_manifest.tsv")
GWAS_DIR = Path(r"D:/metabolic/233/GWAS")
QC_DIR = Path(r"D:/metabolic/233/QC")
OUT_DIR = Path(r"D:/metabolic/233/retry")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(MANIFEST, sep="\t")

    rows = []
    for acc in manifest["study_accession"].astype(str):
        txt_path = GWAS_DIR / f"{acc}.txt"
        tmp_path = GWAS_DIR / f"{acc}.tmp.txt"
        qc_path = QC_DIR / f"{acc}.qc_summary.tsv"
        rows.append(
            {
                "study_accession": acc,
                "txt_exists": txt_path.exists(),
                "tmp_exists": tmp_path.exists(),
                "qc_exists": qc_path.exists(),
                "txt_path": str(txt_path).replace("\\", "/"),
                "tmp_path": str(tmp_path).replace("\\", "/"),
                "qc_path": str(qc_path).replace("\\", "/"),
            }
        )

    status = pd.DataFrame(rows)
    failed = status[(~status["txt_exists"]) | (~status["qc_exists"])].copy()
    failed = failed.sort_values("study_accession").reset_index(drop=True)
    failed.to_csv(OUT_DIR / "external233_failed_retry_manifest.tsv", sep="\t", index=False)

    accessions = failed[["study_accession"]].copy()
    accessions.to_csv(OUT_DIR / "external233_failed_accessions.tsv", sep="\t", index=False)

    # Split into small batches to reduce the impact of transient SSL failures.
    batch_size = 10
    for i in range(0, len(accessions), batch_size):
        batch = accessions.iloc[i : i + batch_size].copy()
        batch_id = i // batch_size + 1
        batch.to_csv(OUT_DIR / f"external233_failed_batch_{batch_id:02d}.tsv", sep="\t", index=False)

    print(f"failed_n={len(failed)}")
    print("wrote", OUT_DIR / "external233_failed_retry_manifest.tsv")


if __name__ == "__main__":
    main()
