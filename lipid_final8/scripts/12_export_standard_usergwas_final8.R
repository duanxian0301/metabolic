library(data.table)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
native_dir <- file.path(root_dir, "step14_native_wsl_usergwas_final8_results")
merged_dir <- file.path(native_dir, "merged_lipid_final8")
std_dir <- file.path(merged_dir, "standard_txt")
sumstats_lookup_file <- file.path(root_dir, "step14_native_wsl_factor_gwas_inputs", "final8_factorGWAS_sumstats.tsv.gz")
bim_lookup_file <- "D:/SMR/g1000/g1000_eur.bim"

dir.create(merged_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(std_dir, showWarnings = FALSE, recursive = TRUE)

factor_files <- c(
  F1 = file.path(merged_dir, "lipid_final8_F1_userGWAS_merged.tsv.gz"),
  F2 = file.path(merged_dir, "lipid_final8_F2_userGWAS_merged.tsv.gz"),
  F3 = file.path(merged_dir, "lipid_final8_F3_userGWAS_merged.tsv.gz")
)

missing_files <- factor_files[!file.exists(factor_files)]
if (length(missing_files) > 0) {
  stop("Merged userGWAS files are required before export:\n", paste(missing_files, collapse = "\n"))
}
if (!file.exists(sumstats_lookup_file)) {
  stop("Missing factor GWAS SNP lookup file: ", sumstats_lookup_file)
}

lookup_dt <- fread(sumstats_lookup_file, select = c("SNP", "A1", "A2", "MAF"))
lookup_dt <- unique(lookup_dt, by = "SNP")

coord_dt <- NULL
if (file.exists(bim_lookup_file)) {
  coord_dt <- fread(
    bim_lookup_file,
    col.names = c("CHR", "SNP", "CM", "BP", "BIM_A1", "BIM_A2"),
    select = c("CHR", "SNP", "BP")
  )
  coord_dt <- unique(coord_dt, by = "SNP")
}

estimate_neff <- function(maf, se) {
  keep <- !is.na(maf) & !is.na(se) & maf >= 0.1 & maf <= 0.4 & se > 0
  if (!any(keep)) {
    return(NA_integer_)
  }
  as.integer(round(mean(1 / (2 * maf[keep] * (1 - maf[keep]) * (se[keep]^2)))))
}

export_one <- function(infile, outfile) {
  dt <- fread(infile)
  dt <- merge(lookup_dt, dt[, .(SNP, est, SE, Pval_Estimate)], by = "SNP", all.x = FALSE, all.y = TRUE)
  if (!is.null(coord_dt)) {
    dt <- merge(coord_dt, dt, by = "SNP", all.y = TRUE)
  }
  n_value <- estimate_neff(dt$MAF, dt$SE)
  out <- data.table(
    SNP = dt$SNP,
    CHR = dt$CHR,
    BP = dt$BP,
    A1 = dt$A1,
    A2 = dt$A2,
    FRQ = dt$MAF,
    BETA = dt$est,
    SE = dt$SE,
    P = dt$Pval_Estimate,
    N = as.integer(n_value)
  )
  fwrite(out, outfile, sep = "\t")
  n_value
}

f1_n <- export_one(
  factor_files[["F1"]],
  file.path(std_dir, "lipid_final8_F1_standard.txt")
)

f2_n <- export_one(
  factor_files[["F2"]],
  file.path(std_dir, "lipid_final8_F2_standard.txt")
)

f3_n <- export_one(
  factor_files[["F3"]],
  file.path(std_dir, "lipid_final8_F3_standard.txt")
)

meta <- data.table(
  factor = c("F1", "F2", "F3"),
  merged_file = unname(factor_files),
  standard_file = c(
    file.path(std_dir, "lipid_final8_F1_standard.txt"),
    file.path(std_dir, "lipid_final8_F2_standard.txt"),
    file.path(std_dir, "lipid_final8_F3_standard.txt")
  ),
  n_used = c(f1_n, f2_n, f3_n),
  notes = "Native WSL GenomicSEM userGWAS rerun on lipid final8 3-factor model; FREQ uses MAF; CHR/BP merged from D:/SMR/g1000/g1000_eur.bim when available; N replaced with factor-specific Neff estimated as mean(1/(2*MAF*(1-MAF)*SE^2)) across SNPs with 0.1<=MAF<=0.4."
)

fwrite(meta, file.path(std_dir, "lipid_final8_standard_txt_metadata.tsv"), sep = "\t")
