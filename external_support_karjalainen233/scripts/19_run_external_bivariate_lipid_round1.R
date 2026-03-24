library(GenomicSEM)
library(data.table)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
std_dir <- file.path(
  root_dir,
  "step14_native_wsl_usergwas_final8_results",
  "merged_lipid_final8",
  "standard_txt"
)
repo_dir <- "D:/codex/metabolic_repo/external_support_karjalainen233/inputs/external_bivariate_round1/lipid"
ascii_work_dir <- "D:/codex/GenomicSEM/metabolic/external_bivariate_lipid_round1_work"
out_dir <- "D:/metabolic/233/ldsc_univariate/external_bivariate_round1_lipid"

dir.create(ascii_work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(ascii_work_dir)

hm3 <- "D:/LDSC/ldsc-master/eur_w_ld_chr/w_hm3.snplist"
ld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"
wld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"

traits_manifest <- fread(file.path(repo_dir, "external_traits_manifest.tsv"))
pairs_dt <- fread(file.path(repo_dir, "requested_pairs.tsv"))

factor_meta <- data.table(
  trait = c("lipid8_F1", "lipid8_F2", "lipid8_F3"),
  source_file = c(
    file.path(std_dir, "lipid_final8_F1_standard.txt"),
    file.path(std_dir, "lipid_final8_F2_standard.txt"),
    file.path(std_dir, "lipid_final8_F3_standard.txt")
  ),
  N = c(599249, 599249, 599249)
)

files <- file.path(ascii_work_dir, paste0(factor_meta$trait, ".txt"))
for (i in seq_len(nrow(factor_meta))) {
  file.copy(factor_meta$source_file[i], files[i], overwrite = TRUE)
}

munge(
  files = files,
  hm3 = hm3,
  trait.names = factor_meta$trait,
  N = factor_meta$N
)

factor_traits <- data.table(
  trait = factor_meta$trait,
  sumstats = file.path(ascii_work_dir, paste0(factor_meta$trait, ".sumstats.gz")),
  category = "factor"
)

external_traits <- traits_manifest[, .(
  trait,
  sumstats = sumstats_file,
  category = "external_trait"
)]

all_traits <- rbindlist(list(factor_traits, external_traits), use.names = TRUE, fill = TRUE)

ldsc_out <- ldsc(
  traits = all_traits$sumstats,
  sample.prev = rep(NA, nrow(all_traits)),
  population.prev = rep(NA, nrow(all_traits)),
  ld = ld,
  wld = wld,
  trait.names = all_traits$trait,
  ldsc.log = file.path(out_dir, "external_bivariate_lipid_round1")
)

saveRDS(ldsc_out, file.path(out_dir, "external_bivariate_lipid_round1.rds"))

S <- as.matrix(ldsc_out$S)
rownames(S) <- colnames(S)
cov_ratio <- tcrossprod(1 / sqrt(diag(S)))
R <- S * cov_ratio
rownames(R) <- rownames(S)
colnames(R) <- colnames(S)

se_cov <- matrix(0, nrow(S), ncol(S), dimnames = dimnames(S))
se_cov[lower.tri(se_cov, diag = TRUE)] <- sqrt(diag(ldsc_out$V))
se_cov[upper.tri(se_cov)] <- t(se_cov)[upper.tri(se_cov)]

scale_o <- lower.tri(cov_ratio, diag = TRUE)
scale_vec <- cov_ratio[scale_o]
V_std <- ldsc_out$V * tcrossprod(scale_vec)
se_rg <- matrix(0, nrow(R), ncol(R), dimnames = dimnames(R))
se_rg[lower.tri(se_rg, diag = TRUE)] <- sqrt(diag(V_std))
se_rg[upper.tri(se_rg)] <- t(se_rg)[upper.tri(se_rg)]

diag_dt <- data.table(
  trait = colnames(S),
  h2 = diag(S),
  h2_se = diag(se_cov),
  h2_z = diag(S) / diag(se_cov),
  intercept = diag(ldsc_out$I)
)

trait_info <- traits_manifest[, .(
  trait,
  biomarker_name,
  internal_trait,
  anchor_factor,
  validation_priority,
  study_accession
)]

pairs_dt <- merge(
  pairs_dt,
  trait_info,
  by.x = "target_trait",
  by.y = "trait",
  all.x = TRUE,
  sort = FALSE
)

pairs_dt[, `:=`(
  covariance = mapply(function(x, y) S[x, y], factor, target_trait),
  covariance_se = mapply(function(x, y) se_cov[x, y], factor, target_trait),
  rg = mapply(function(x, y) R[x, y], factor, target_trait),
  rg_se = mapply(function(x, y) se_rg[x, y], factor, target_trait)
)]

pairs_dt[, z_cov := covariance / covariance_se]
pairs_dt[, p_cov := 2 * pnorm(abs(z_cov), lower.tail = FALSE)]
pairs_dt[, z_rg := rg / rg_se]
pairs_dt[, p_rg := 2 * pnorm(abs(z_rg), lower.tail = FALSE)]
pairs_dt[, abs_rg := abs(rg)]

summary_dt <- pairs_dt[, .(
  n_pairs = .N,
  mean_rg = mean(rg, na.rm = TRUE),
  mean_abs_rg = mean(abs_rg, na.rm = TRUE),
  median_rg = median(rg, na.rm = TRUE),
  median_abs_rg = median(abs_rg, na.rm = TRUE),
  min_p_rg = min(p_rg, na.rm = TRUE)
), by = .(factor, pair_set, expected_alignment)]

contrast_dt <- dcast(
  pairs_dt[, .(
    mean_rg = mean(rg, na.rm = TRUE),
    mean_abs_rg = mean(abs_rg, na.rm = TRUE)
  ), by = .(factor, pair_set)],
  factor ~ pair_set,
  value.var = c("mean_rg", "mean_abs_rg")
)

fwrite(traits_manifest, file.path(out_dir, "external_bivariate_lipid_round1_traits_manifest.tsv"), sep = "\t")
fwrite(diag_dt, file.path(out_dir, "external_bivariate_lipid_round1_h2.tsv"), sep = "\t")
fwrite(as.data.table(R, keep.rownames = "trait"), file.path(out_dir, "external_bivariate_lipid_round1_rg_matrix.tsv"), sep = "\t")
fwrite(pairs_dt, file.path(out_dir, "external_bivariate_lipid_round1_requested_pairs.tsv"), sep = "\t")
fwrite(summary_dt, file.path(out_dir, "external_bivariate_lipid_round1_pairset_summary.tsv"), sep = "\t")
fwrite(contrast_dt, file.path(out_dir, "external_bivariate_lipid_round1_factor_contrast.tsv"), sep = "\t")

message("External lipid round1 bivariate LDSC written to: ", out_dir)
