from itertools import permutations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(r"D:\metabolic\GWAS\genomicgem_main_zgt4_nonproportion")
OUT_DIR = Path(r"D:\codex\metabolic_repo\lipid_final8\manuscript\figures_lipid_nonlipid_final8")
OUT_DIR.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="white", context="talk")
plt.rcParams["font.family"] = "DejaVu Sans"
plt.rcParams["axes.titleweight"] = "bold"
plt.rcParams["figure.dpi"] = 180

MODULE_COLORS = {
    "lipid": "#B24745",
    "nonlipid": "#2F6C8F",
}

FACTOR_COLORS = {
    "F1_TG_rich_axis": "#C85C3A",
    "F2_HDL_core_axis": "#5B8E7D",
    "F3_CE_structural_axis": "#7B6AA8",
    "F1_ketone_axis": "#B05A2B",
    "F2_amino_acid_axis": "#417A56",
    "F3_energy_bridge_axis": "#2E6E9E",
}


def save_fig(fig, stem):
    fig.savefig(OUT_DIR / f"{stem}.png", bbox_inches="tight", dpi=300)
    fig.savefig(OUT_DIR / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def best_factor_mapping(loadings_df, manifest_df):
    latent_factors = sorted(loadings_df["factor"].unique())
    biological_factors = list(dict.fromkeys(manifest_df["final_factor_name"]))
    trait_to_group = dict(zip(manifest_df["trait_code"], manifest_df["final_factor_name"]))
    score = pd.DataFrame(0.0, index=latent_factors, columns=biological_factors)

    for latent in latent_factors:
        sub = loadings_df[loadings_df["factor"] == latent].copy()
        sub["bio_factor"] = sub["trait"].map(trait_to_group)
        for bio in biological_factors:
            vals = sub.loc[sub["bio_factor"] == bio, "std.all"].abs()
            score.loc[latent, bio] = vals.mean() if len(vals) else -np.inf

    best_perm = None
    best_score = -np.inf
    for perm in permutations(biological_factors, len(latent_factors)):
        total = sum(score.loc[latent, bio] for latent, bio in zip(latent_factors, perm))
        if total > best_score:
            best_score = total
            best_perm = perm
    return dict(zip(latent_factors, best_perm))


def prepare_loading_matrix(manifest_path, loadings_path):
    manifest = pd.read_csv(manifest_path, sep="\t")
    loadings = pd.read_csv(loadings_path, sep="\t")
    mapping = best_factor_mapping(loadings, manifest)
    loadings["factor_label"] = loadings["factor"].map(mapping)
    matrix = loadings.pivot(index="trait", columns="factor_label", values="std.all")

    factor_order = list(dict.fromkeys(manifest["final_factor_name"]))
    order_df = manifest[["trait_code", "final_factor_name"]].drop_duplicates().copy()
    order_df["factor_rank"] = order_df["final_factor_name"].map({f: i for i, f in enumerate(factor_order)})
    primary_abs = (
        loadings.assign(factor_label=loadings["factor"].map(mapping))
        .loc[lambda d: d["factor_label"] == d["ultra_pure_factor"]]
        .set_index("trait")["std.all"]
        .abs()
    )
    order_df["primary_abs"] = order_df["trait_code"].map(primary_abs).fillna(0)
    order_df = order_df.sort_values(["factor_rank", "primary_abs"], ascending=[True, False])
    row_order = order_df["trait_code"].tolist()
    matrix = matrix.reindex(index=row_order, columns=factor_order)
    return matrix, manifest


def read_square_matrix(path, labels):
    df = pd.read_csv(path, index_col=0)
    df.index = labels
    df.columns = labels
    return df


def read_factor_rg(path, labels):
    df = pd.read_csv(path, sep="\t")
    df.index = labels
    df.columns = labels
    return df


def add_panel_label(ax, label):
    ax.text(-0.12, 1.08, label, transform=ax.transAxes, fontsize=18, fontweight="bold", va="bottom")


def make_stage_reduction_figure():
    stage_df = pd.DataFrame(
        [
            ("lipid", "Main112", 112),
            ("lipid", "Candidate", 94),
            ("lipid", "Compact", 20),
            ("lipid", "Ultrapure", 10),
            ("lipid", "Final", 8),
            ("nonlipid", "Main112", 112),
            ("nonlipid", "Candidate", 18),
            ("nonlipid", "Compact", 15),
            ("nonlipid", "Ultrapure", 8),
            ("nonlipid", "Final", 8),
        ],
        columns=["module", "stage", "n_traits"],
    )
    stage_order = ["Main112", "Candidate", "Compact", "Ultrapure", "Final"]
    x = np.arange(len(stage_order))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)
    for ax, module in zip(axes, ["lipid", "nonlipid"]):
        sub = stage_df[stage_df["module"] == module].set_index("stage").loc[stage_order]
        color = MODULE_COLORS[module]
        ax.plot(x, sub["n_traits"], marker="o", markersize=10, linewidth=3, color=color)
        ax.fill_between(x, sub["n_traits"], alpha=0.18, color=color)
        for xi, yi in zip(x, sub["n_traits"]):
            ax.text(xi, yi + 3, f"{int(yi)}", ha="center", va="bottom", fontsize=12, fontweight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(stage_order)
        ax.set_title(f"{module.capitalize()} staged reduction", color=color)
        ax.set_xlabel("")
        ax.set_ylabel("Number of traits")
        ax.set_ylim(0, 125)
        ax.grid(axis="y", alpha=0.2)
        sns.despine(ax=ax)
    add_panel_label(axes[0], "A")
    add_panel_label(axes[1], "B")
    fig.suptitle("Module construction from the shared Main112 universe", y=1.02, fontsize=20, fontweight="bold")
    save_fig(fig, "Figure1_staged_reduction")


def make_loading_figure():
    lipid_matrix, lipid_manifest = prepare_loading_matrix(
        ROOT / "final_model_lipid_module_final8" / "final8_trait_manifest.tsv",
        ROOT / "step13_efa_esem_lipid_module_final8" / "ALL_3factor_loadings.tsv",
    )
    nonlipid_matrix, nonlipid_manifest = prepare_loading_matrix(
        ROOT / "final_model_nonlipid_module_final8" / "final8_trait_manifest.tsv",
        ROOT / "step20_efa_esem_nonlipid_module_ultrapure8" / "ALL_3factor_loadings.tsv",
    )

    fig, axes = plt.subplots(1, 2, figsize=(19, 8), gridspec_kw={"wspace": 0.42})
    cmap = sns.diverging_palette(240, 10, as_cmap=True)

    for ax, matrix, title, module in [
        (axes[0], lipid_matrix, "Lipid final8 standardized loadings", "lipid"),
        (axes[1], nonlipid_matrix, "Nonlipid final8 standardized loadings", "nonlipid"),
    ]:
        sns.heatmap(
            matrix,
            ax=ax,
            cmap=cmap,
            center=0,
            vmin=-1,
            vmax=1,
            annot=True,
            fmt=".2f",
            linewidths=0.5,
            cbar=ax is axes[1],
            cbar_kws={"shrink": 0.8, "label": "Standardized loading", "pad": 0.04},
        )
        ax.set_title(title, color=MODULE_COLORS[module])
        ax.set_xlabel("Latent factor")
        ax.set_ylabel("Trait")
        ax.tick_params(axis="x", labelrotation=90, labelsize=12)
        ax.tick_params(axis="y", labelrotation=0, labelsize=11, pad=2)
    add_panel_label(axes[0], "A")
    add_panel_label(axes[1], "B")
    fig.suptitle("Final 3-factor loading structure", y=1.02, fontsize=20, fontweight="bold")
    save_fig(fig, "Figure2_final_loadings")


def make_trait_rg_figure():
    lipid_labels = ["M_HDL_TG", "VLDL_size", "MUFA", "S_VLDL_TG", "ApoA1", "HDL_CE", "XS_VLDL_FC", "VLDL_CE"]
    nonlipid_labels = ["Acetoacetate", "bOHbutyrate", "Leu", "Phe", "Val", "Acetate", "Glucose", "Lactate"]

    lipid_rg = read_square_matrix(ROOT / "step12_ldsc_lipid_module_final8" / "lipid_module_final8_rg_matrix.csv", lipid_labels)
    nonlipid_rg = read_square_matrix(ROOT / "step19_ldsc_nonlipid_module_ultrapure8" / "nonlipid_module_ultrapure8_rg_matrix.csv", nonlipid_labels)

    fig, axes = plt.subplots(1, 2, figsize=(19, 7), gridspec_kw={"wspace": 0.38})
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    for ax, matrix, title, module in [
        (axes[0], lipid_rg, "Lipid final8 trait-level rg", "lipid"),
        (axes[1], nonlipid_rg, "Nonlipid final8 trait-level rg", "nonlipid"),
    ]:
        sns.heatmap(
            matrix,
            ax=ax,
            cmap=cmap,
            center=0,
            vmin=-1,
            vmax=1,
            annot=True,
            fmt=".2f",
            square=True,
            linewidths=0.5,
            cbar=ax is axes[1],
            cbar_kws={"shrink": 0.8, "label": "Genetic correlation (rg)", "pad": 0.04},
        )
        ax.set_title(title, color=MODULE_COLORS[module])
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.tick_params(axis="x", labelrotation=90, labelsize=11)
        ax.tick_params(axis="y", labelrotation=0, labelsize=10, pad=2)
    add_panel_label(axes[0], "A")
    add_panel_label(axes[1], "B")
    fig.suptitle("Trait-level genetic correlation structure in the final panels", y=1.02, fontsize=20, fontweight="bold")
    save_fig(fig, "Figure3_final_trait_rg")


def make_factor_validation_figure():
    lipid_uni = pd.read_csv(ROOT / "supplement_lipid_final8" / "Supplementary_Table_S5_lipid_final8_univariate_ldsc_detailed.tsv", sep="\t")
    nonlipid_uni = pd.read_csv(ROOT / "supplement_nonlipid_final8" / "Supplementary_Table_S5_nonlipid_final8_univariate_ldsc_detailed.tsv", sep="\t")
    lipid_bi = pd.read_csv(ROOT / "supplement_lipid_final8" / "Supplementary_Table_S6_lipid_final8_bivariate_ldsc_detailed.tsv", sep="\t")
    nonlipid_bi = pd.read_csv(ROOT / "supplement_nonlipid_final8" / "Supplementary_Table_S6_nonlipid_final8_bivariate_ldsc_detailed.tsv", sep="\t")

    lipid_factors = ["lipid8_F1", "lipid8_F2", "lipid8_F3"]
    nonlipid_factors = ["nonlipid8_F1", "nonlipid8_F2", "nonlipid8_F3"]
    lipid_factor_rg = read_factor_rg(ROOT / "step15_ldsc_validation_lipid_final8" / "lipid_final8_factor_rg_matrix.tsv", lipid_factors)
    nonlipid_factor_rg = read_factor_rg(ROOT / "step23_ldsc_validation_nonlipid_final8" / "nonlipid_final8_factor_rg_matrix.tsv", nonlipid_factors)

    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.35, wspace=0.25)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    for ax, uni, module, title in [
        (ax1, lipid_uni, "lipid", "Lipid factor SNP heritability"),
        (ax2, nonlipid_uni, "nonlipid", "Nonlipid factor SNP heritability"),
    ]:
        x = np.arange(len(uni))
        ax.bar(x, uni["h2"], color=MODULE_COLORS[module], alpha=0.85)
        ax.errorbar(x, uni["h2"], yerr=uni["h2_se"], fmt="none", ecolor="black", elinewidth=1.5, capsize=4)
        upper = (uni["h2"] + uni["h2_se"]).max()
        ax.set_ylim(0, upper + 0.03)
        for xi, yi, zi in zip(x, uni["h2"], uni["h2_z"]):
            ax.text(
                xi,
                yi + max(uni["h2_se"]) + 0.002,
                f"h2={yi:.3f}\nZ={zi:.2f}",
                ha="center",
                va="bottom",
                fontsize=10,
            )
        ax.set_xticks(x)
        ax.set_xticklabels(uni["trait"], rotation=0)
        ax.set_ylabel("SNP heritability (h2)")
        ax.set_title(title, color=MODULE_COLORS[module])
        ax.grid(axis="y", alpha=0.2)
        sns.despine(ax=ax)

    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    for ax, matrix, title, module in [
        (ax3, lipid_factor_rg, "Lipid factor-level rg", "lipid"),
        (ax4, nonlipid_factor_rg, "Nonlipid factor-level rg", "nonlipid"),
    ]:
        sns.heatmap(
            matrix,
            ax=ax,
            cmap=cmap,
            center=0,
            vmin=-1,
            vmax=1,
            annot=True,
            fmt=".2f",
            square=True,
            linewidths=0.5,
            cbar=ax is ax4,
            cbar_kws={"shrink": 0.85, "label": "Genetic correlation (rg)"},
        )
        ax.set_title(title, color=MODULE_COLORS[module])
        ax.set_xlabel("")
        ax.set_ylabel("")

    add_panel_label(ax1, "A")
    add_panel_label(ax2, "B")
    add_panel_label(ax3, "C")
    add_panel_label(ax4, "D")
    fig.suptitle("Factor-level GWAS validation across modules", y=1.01, fontsize=20, fontweight="bold")
    save_fig(fig, "Figure4_factor_validation")


def write_overview():
    overview = """Figure outputs
- Figure1_staged_reduction: staged reduction from shared Main112 to final8 for lipid and nonlipid.
- Figure2_final_loadings: standardized loading heatmaps for the final 3-factor models.
- Figure3_final_trait_rg: trait-level genetic-correlation heatmaps for the final8 panels.
- Figure4_factor_validation: factor-level SNP heritability barplots and factor-level rg heatmaps.

Design notes
- Lipid is encoded in muted red tones.
- Nonlipid is encoded in blue tones.
- Heatmaps use a shared diverging scale centered at zero for direct comparison.
"""
    (OUT_DIR / "figure_overview.txt").write_text(overview, encoding="utf-8")


def main():
    make_stage_reduction_figure()
    make_loading_figure()
    make_trait_rg_figure()
    make_factor_validation_figure()
    write_overview()
    print(f"Wrote figures to: {OUT_DIR}")


if __name__ == "__main__":
    main()
