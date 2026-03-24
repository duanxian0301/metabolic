library(GenomicSEM)
library(data.table)

root_dir <- "D:/metabolic/GWAS/genomicgem_main_zgt4_nonproportion"
repo_dir <- "D:/codex/metabolic_repo/lipid_final8"
std_dir <- file.path(
  root_dir,
  "step14_native_wsl_usergwas_final8_results",
  "merged_lipid_final8",
  "standard_txt"
)
ascii_work_dir <- "D:/codex/GenomicSEM/metabolic/step24_internal_lipid_ldsc_work"
out_dir <- file.path(root_dir, "step24_internal_bivariate_lipid_final8")

traits_manifest_path <- file.path(repo_dir, "inputs", "validation", "validation_traits.tsv")
pairs_path <- file.path(repo_dir, "inputs", "validation", "requested_pairs.tsv")

dir.create(ascii_work_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
setwd(ascii_work_dir)

hm3 <- "D:/LDSC/ldsc-master/eur_w_ld_chr/w_hm3.snplist"
ld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"
wld <- "D:/LDSC/ldsc-master/eur_w_ld_chr/"

traits_manifest <- fread(traits_manifest_path)
pairs_dt <- fread(pairs_path)

factor_meta <- data.table(
  trait = c("lipid8_F1", "lipid8_F2", "lipid8_F3"),
  source_file = c(
    file.path(std_dir, "lipid_final8_F1_standard.txt"),
    file.path(std_dir, "lipid_final8_F2_standard.txt"),
    file.path(std_dir, "lipid_final8_F3_standard.txt")
  ),
  N = c(599249, 599249, 599249)
)

if (!all(file.exists(factor_meta$source_file))) {
  stop("Missing standard factor files in: ", std_dir)
}

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

support_traits <- traits_manifest[category != "factor", .(
  trait,
  sumstats = sumstats_file,
  category
)]

all_traits <- rbindlist(list(factor_traits, support_traits), use.names = TRUE, fill = TRUE)

if (!all(file.exists(all_traits$sumstats))) {
  missing_sumstats <- all_traits[!file.exists(sumstats), trait]
  stop("Missing .sumstats.gz for: ", paste(missing_sumstats, collapse = ", "))
}

ldsc_out <- ldsc(
  traits = all_traits$sumstats,
  sample.prev = rep(NA, nrow(all_traits)),
  population.prev = rep(NA, nrow(all_traits)),
  ld = ld,
  wld = wld,
  trait.names = all_traits$trait,
  ldsc.log = file.path(out_dir, "lipid_internal_bivariate_ldsc")
)

saveRDS(ldsc_out, file.path(out_dir, "lipid_internal_bivariate_ldsc.rds"))

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
  category,
  subset,
  anchor_factor,
  biomarker_name,
  group
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

factor_contrast_dt <- pairs_dt[, .(
  mean_rg = mean(rg, na.rm = TRUE),
  mean_abs_rg = mean(abs_rg, na.rm = TRUE)
), by = .(factor, pair_set)]
dcast_dt <- dcast(
  factor_contrast_dt,
  factor ~ pair_set,
  value.var = c("mean_rg", "mean_abs_rg")
)

fwrite(traits_manifest, file.path(out_dir, "lipid_internal_validation_traits_manifest.tsv"), sep = "\t")
fwrite(diag_dt, file.path(out_dir, "lipid_internal_bivariate_ldsc_h2.tsv"), sep = "\t")
fwrite(as.data.table(R, keep.rownames = "trait"), file.path(out_dir, "lipid_internal_bivariate_ldsc_rg_matrix.tsv"), sep = "\t")
fwrite(pairs_dt, file.path(out_dir, "lipid_internal_bivariate_ldsc_requested_pairs.tsv"), sep = "\t")
fwrite(summary_dt, file.path(out_dir, "lipid_internal_bivariate_ldsc_pairset_summary.tsv"), sep = "\t")
fwrite(dcast_dt, file.path(out_dir, "lipid_internal_bivariate_ldsc_factor_contrast.tsv"), sep = "\t")

message("Lipid internal validation written to: ", out_dir)
