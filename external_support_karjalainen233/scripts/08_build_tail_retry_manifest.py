from pathlib import Path

import pandas as pd


MANIFEST = Path(r"D:/metabolic/233/manifests/external233_processing_manifest.tsv")
GWAS_DIR = Path(r"D:/metabolic/233/GWAS")
QC_DIR = Path(r"D:/metabolic/233/QC")
OUT_DIR = Path(r"D:/metabolic/233/retry_tail")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(MANIFEST, sep="\t")
    rows = []
    for acc in manifest["study_accession"].astype(str):
        txt_exists = (GWAS_DIR / f"{acc}.txt").exists()
        tmp_exists = (GWAS_DIR / f"{acc}.tmp.txt").exists()
        qc_exists = (QC_DIR / f"{acc}.qc_summary.tsv").exists()
        if (not txt_exists) or (not qc_exists):
            rows.append({"study_accession": acc, "txt_exists": txt_exists, "tmp_exists": tmp_exists, "qc_exists": qc_exists})

    pending = pd.DataFrame(rows).sort_values("study_accession").reset_index(drop=True)
    pending.to_csv(OUT_DIR / "external233_tail_retry_manifest.tsv", sep="\t", index=False)
    pending[["study_accession"]].to_csv(OUT_DIR / "external233_tail_accessions.tsv", sep="\t", index=False)

    batch_size = 3
    for i in range(0, len(pending), batch_size):
        batch = pending.iloc[i : i + batch_size][["study_accession"]].copy()
        batch_id = i // batch_size + 1
        batch.to_csv(OUT_DIR / f"external233_tail_batch_{batch_id:02d}.tsv", sep="\t", index=False)

    print(f"tail_pending_n={len(pending)}")


if __name__ == "__main__":
    main()
