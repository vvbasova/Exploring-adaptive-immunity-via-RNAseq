import os
import glob
import itertools

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from statsmodels.stats.multitest import multipletests

CORRELATION_THRESHOLD = 0.65
P_VALUE_THRESHOLD = 0.05


def correlation_table_maker_cross(
    df1: pd.DataFrame, df2: pd.DataFrame, localization: str, output_dir: str
) -> None:
    """
    Calculates pairwise Pearson correlations between columns of two dataframes
    and writes the result to a CSV file with multiple testing correction.

    Parameters:
        df1 (pd.DataFrame): First input dataframe.
        df2 (pd.DataFrame): Second input dataframe.
        localization (str): Label used for naming the output file.
        output_dir (str): Path to directory where results will be saved.
    """
    cols1 = df1.columns
    cols2 = df2.columns

    combinations = list(itertools.product(cols1, cols2))
    results = []

    for gene1, gene2 in combinations:
        try:
            corr, pval = pearsonr(df1[gene1], df2[gene2])
        except Exception:
            corr, pval = np.nan, np.nan
        results.append((gene1, gene2, corr, pval))

    results_df = pd.DataFrame(
        results, columns=["Var1", "Var2", "correlation", "p_value"]
    )
    results_df["p_adj"] = multipletests(results_df["p_value"], method="fdr_bh")[1]
    results_df["-log10(p_adj)"] = -np.log10(results_df["p_adj"])
    results_df["pair"] = results_df["Var1"] + " / " + results_df["Var2"]

    output_dir = os.path.expanduser(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    results_df.to_csv(os.path.join(output_dir, f"{localization}_cor.csv"), index=False)


def process_localizations(
    df1: pd.DataFrame, df2: pd.DataFrame, output_dir: str
) -> None:
    """
    Splits data by localization ('project_id') and computes correlation
    tables for each localization.

    Parameters:
        df1 (pd.DataFrame): First dataframe with 'project_id' and 'sample' columns.
        df2 (pd.DataFrame): Second dataframe with 'project_id' and 'sample' columns.
        output_dir (str): Directory to save correlation CSVs.
    """
    localizations = df1["project_id"].unique()

    for loco in localizations:
        df_loco1 = df1[df1["project_id"] == loco].copy()
        df_loco2 = df2[df2["project_id"] == loco].copy()

        drop_cols = ["project_id", "sample"]
        df_loco1 = df_loco1.drop(columns=[col for col in drop_cols if col in df_loco1])
        df_loco2 = df_loco2.drop(columns=[col for col in drop_cols if col in df_loco2])

        df_loco1 = df_loco1.loc[:, (df_loco1 != 0).any(axis=0)]
        df_loco2 = df_loco2.loc[:, (df_loco2 != 0).any(axis=0)]

        correlation_table_maker_cross(df_loco1, df_loco2, loco, output_dir)


def make_clustermaps(correlation_dir: str) -> None:
    """
    Generates cluster maps from correlation CSVs and saves them as PNGs.

    Parameters:
        correlation_dir (str): Directory containing correlation CSV files.
    """
    correlation_dir = os.path.expanduser(correlation_dir)
    clustermap_dir = os.path.join(correlation_dir, "clustermaps")
    os.makedirs(clustermap_dir, exist_ok=True)

    csv_files = glob.glob(os.path.join(correlation_dir, "*.csv"))

    for file_path in csv_files:
        key = os.path.splitext(os.path.basename(file_path))[0]
        df = pd.read_csv(file_path)

        pivot_table = df.pivot(index="Var1", columns="Var2", values="correlation")
        p_adj_table = df.pivot(index="Var1", columns="Var2", values="p_adj")

        significant_rows = p_adj_table.index[(p_adj_table < 0.05).any(axis=1)]
        significant_cols = p_adj_table.columns[(p_adj_table < 0.05).any(axis=0)]

        filtered_pivot_table = pivot_table.loc[significant_rows, significant_cols]
        filtered_p_adj_table = p_adj_table.loc[significant_rows, significant_cols]

        if filtered_pivot_table.shape[0] < 2 or filtered_pivot_table.shape[1] < 2:
            continue

        mask = filtered_p_adj_table >= 0.05

        cmap = sns.diverging_palette(240, 10, as_cmap=True)
        cmap.set_bad(color="lightgrey")

        clustergrid = sns.clustermap(
            filtered_pivot_table, cmap=cmap, mask=mask, linewidths=0.5, vmin=-1, vmax=1
        )

        clustergrid.fig.suptitle(f"{key}", fontsize=16)

        output_path = os.path.join(clustermap_dir, f"{key}.png")
        clustergrid.fig.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(clustergrid.fig)


def make_heatmaps(correlation_dir: str) -> None:
    """
    Creates heatmaps of correlation matrices with significant values only.

    Parameters:
        correlation_dir (str): Directory containing *_cor.csv files.
    """
    correlation_dir = os.path.expanduser(correlation_dir)
    heatmap_dir = os.path.join(correlation_dir, "heatmaps")
    os.makedirs(heatmap_dir, exist_ok=True)

    csv_files = glob.glob(os.path.join(correlation_dir, "*_cor.csv"))

    all_var1_genes = set()
    all_var2_genes = set()

    for file_path in csv_files:
        df = pd.read_csv(file_path)
        all_var1_genes.update(df["Var1"].dropna().unique())
        all_var2_genes.update(df["Var2"].dropna().unique())

    var1_order = sorted(all_var1_genes)
    var2_order = sorted(all_var2_genes)

    for file_path in csv_files:
        filename = os.path.splitext(os.path.basename(file_path))[0]
        project_id = filename.replace("_cor", "")

        df = pd.read_csv(file_path)
        pivot_table = df.pivot(index="Var1", columns="Var2", values="correlation")
        p_adj_table = df.pivot(index="Var1", columns="Var2", values="p_adj")

        pivot_table = pivot_table.reindex(index=var1_order, columns=var2_order)
        p_adj_table = p_adj_table.reindex(index=var1_order, columns=var2_order)

        significant_rows = p_adj_table.index[(p_adj_table < 0.05).any(axis=1)]
        significant_cols = p_adj_table.columns[(p_adj_table < 0.05).any(axis=0)]

        if len(significant_rows) < 2 or len(significant_cols) < 2:
            continue

        pivot_table = pivot_table.loc[significant_rows, significant_cols]
        p_adj_table = p_adj_table.loc[significant_rows, significant_cols]

        mask = p_adj_table >= 0.05
        pivot_table = pivot_table.fillna(0)

        cmap = sns.color_palette("coolwarm", as_cmap=True)
        cmap.set_bad(color="lightgrey")

        plt.figure(figsize=(12, 6))
        ax = sns.heatmap(
            pivot_table,
            cmap=cmap,
            mask=mask,
            linewidths=0.3,
            linecolor="white",
            cbar_kws={"label": "Correlation"},
            vmin=-1,
            vmax=1,
        )

        ax.set_title(f"{project_id}", fontsize=14)
        ax.set_xlabel("Var2 genes")
        ax.set_ylabel("Var1 genes")
        plt.xticks(rotation=90)
        plt.yticks(rotation=0)
        plt.tight_layout()

        output_path = os.path.join(heatmap_dir, f"{project_id}_heatmap.png")
        plt.savefig(output_path, dpi=150)
        plt.close()


def load_significant_pairs(file_path: str) -> set:
    """
    Load and filter significant correlation pairs from CSV file.
    """
    df = pd.read_csv(file_path, sep=",")
    filtered = df[
        (df["correlation"].abs() > CORRELATION_THRESHOLD)
        & (df["p_adj"] < P_VALUE_THRESHOLD)
    ]
    return set(zip(filtered["Var1"], filtered["Var2"]))


def plot_similarity_clustermap(group_dir: str, title: str) -> None:
    """
    Display a clustermap showing similarity between diagnoses based on shared
    significant correlation pairs.

    Parameters:
        group_dir (str): Path to the folder with *_cor.csv files.
        title (str): Title for the plot.
    """
    group_dir = os.path.expanduser(group_dir)
    diagnosis_files = [f for f in os.listdir(group_dir) if f.endswith("_cor.csv")]
    diagnoses = [f.replace("_cor.csv", "") for f in diagnosis_files]

    diag_to_pairs = {
        diag: load_significant_pairs(os.path.join(group_dir, file))
        for diag, file in zip(diagnoses, diagnosis_files)
    }

    similarity_matrix = pd.DataFrame(0.0, index=diagnoses, columns=diagnoses)

    for diag1, diag2 in itertools.combinations(diagnoses, 2):
        common_pairs = diag_to_pairs[diag1] & diag_to_pairs[diag2]
        total_significant_pairs = len(diag_to_pairs[diag1]) + len(diag_to_pairs[diag2])
        geometric_mean_filtered = (
            np.sqrt(total_significant_pairs / 2) if total_significant_pairs > 0 else 1
        )
        normalized_value = (
            len(common_pairs) / geometric_mean_filtered
            if geometric_mean_filtered > 0
            else 0
        )
        similarity_matrix.loc[diag1, diag2] = float(normalized_value)
        similarity_matrix.loc[diag2, diag1] = float(normalized_value)

    np.fill_diagonal(similarity_matrix.values, 0.0)

    cmap = sns.clustermap(
        similarity_matrix.astype(float),
        annot=True,
        fmt=".1f",
        cmap="coolwarm",
        linewidths=0.5,
        vmin=0,
        vmax=1,
        annot_kws={"size": 6},
        figsize=(10, 8),
    )
    cmap.ax_heatmap.set_title(title, fontsize=14, pad=20)
    plt.show()
