from collections import Counter
from typing import List

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import fisher_exact, hypergeom
from statsmodels.stats.multitest import multipletests
from sklearn.linear_model import LinearRegression


def fisher_by_sample(df: pd.DataFrame, target_gene: str) -> pd.DataFrame:
    """
    Perform Fisher's exact test per project to evaluate enrichment
    of a target antigen gene across samples.

    Parameters:
        df (pd.DataFrame): Input DataFrame containing columns
            ['vdj_antigen.gene', 'tcr_sample', 'tcr_project_id'].
        target_gene (str): The gene of interest to test for enrichment.

    Returns:
        pd.DataFrame: DataFrame with enriched projects (adjusted p < 0.05),
            including odds ratios and p-values.
    """
    df = df.dropna(subset=["vdj_antigen.gene", "tcr_sample", "tcr_project_id"])
    df["vdj_antigen.gene"] = df["vdj_antigen.gene"].str.strip().str.upper()
    target_gene = target_gene.strip().upper()

    sample_gene_presence = (
        df[df["vdj_antigen.gene"] == target_gene]
        .groupby("tcr_sample")
        .size()
        .gt(0)
        .astype(int)
        .rename("has_target_gene")
    )

    sample_info = (
        df[["tcr_sample", "tcr_project_id"]].drop_duplicates().set_index("tcr_sample")
    )
    sample_info = sample_info.join(sample_gene_presence).fillna(0)
    sample_info["has_target_gene"] = sample_info["has_target_gene"].astype(int)

    results = []
    for project in sample_info["tcr_project_id"].unique():
        in_project = sample_info[sample_info["tcr_project_id"] == project]
        out_project = sample_info[sample_info["tcr_project_id"] != project]

        a = in_project["has_target_gene"].sum()
        b = len(in_project) - a
        c = out_project["has_target_gene"].sum()
        d = len(out_project) - c

        table = [[a, b], [c, d]]
        try:
            odds_ratio, p_value = fisher_exact(table)
        except Exception:
            odds_ratio, p_value = float("nan"), 1.0

        results.append(
            {
                "project_id": project,
                "gene": target_gene,
                "odds_ratio": odds_ratio,
                "p_value": p_value,
                "a": a,
                "b": b,
                "c": c,
                "d": d,
            }
        )

    result_df = pd.DataFrame(results)
    result_df["p_value_adj"] = multipletests(result_df["p_value"], method="fdr_bh")[1]

    return result_df[result_df["p_value_adj"] < 0.05]


def analyse_locus(
    df_locus: pd.DataFrame,
    locus_label: str,
    min_cdr3_per_epitope: int = 30,
    epitope_size_range: tuple = (100, 150),
    top_labels: int = 10,
    plot: bool = True,
) -> pd.DataFrame:
    """
    Analyze project–epitope overlap for a given locus using weighted clone counts
    and visualize expected vs observed overlap with residual enrichment.

    Parameters:
        df_locus (pd.DataFrame): DataFrame with TCR sequencing data, including
            'tcr_sample', 'tcr_cdr3', 'tcr_project_id', 'tcr_duplicate_count',
            'vdj_antigen.epitope', 'vdj_antigen.gene', 'vdj_antigen.species'.
        locus_label (str): Label for the plot title.
        min_cdr3_per_epitope (int): Minimum number of unique CDR3s per epitope to include.
        epitope_size_range (tuple): Acceptable range for clone counts per epitope.
        top_labels (int): Number of top enriched/depleted projects to label on the plot.
        plot (bool): Whether to generate and show the scatter plot.

    Returns:
        pd.DataFrame: DataFrame with enrichment statistics, residuals and predictions.
    """
    df_locus = df_locus.copy()

    df_locus["n_epitopes"] = (
        df_locus.groupby(["tcr_sample", "tcr_cdr3"])["vdj_antigen.epitope"]
        .transform("nunique")
        .fillna(1)
    )
    df_locus["weighted_count"] = (
        df_locus["tcr_duplicate_count"] / df_locus["n_epitopes"]
    )

    all_clones = set(df_locus["tcr_cdr3"])
    M_total = len(all_clones)

    proj_stats = df_locus.groupby("tcr_project_id").agg(
        n_samples=("tcr_sample", "nunique"),
        dup_total=("weighted_count", "sum"),
        clones_unique=("tcr_cdr3", "nunique"),
    )
    proj_stats["dup_norm"] = proj_stats["dup_total"] / proj_stats["n_samples"]

    epi_stats = df_locus.groupby("vdj_antigen.epitope").agg(
        n_epitope_unique=("tcr_cdr3", "nunique"),
        dup_epitope=("weighted_count", "sum"),
        antigen_species=(
            "vdj_antigen.species",
            lambda x: Counter(x).most_common(1)[0][0],
        ),
        antigen_gene=("vdj_antigen.gene", lambda x: Counter(x).most_common(1)[0][0]),
    )
    epi_stats = epi_stats[
        (epi_stats["n_epitope_unique"] >= min_cdr3_per_epitope)
        & (epi_stats["n_epitope_unique"].between(*epitope_size_range))
    ]

    epi2clones = df_locus.groupby("vdj_antigen.epitope")["tcr_cdr3"].apply(set)
    epi2clones = epi2clones.loc[epi_stats.index]

    rows = []
    for epitope, epi_set in epi2clones.items():
        n_epi = len(epi_set)
        dup_epi = epi_stats.loc[epitope, "dup_epitope"]

        for proj, pstats in proj_stats.iterrows():
            proj_set = set(df_locus[df_locus["tcr_project_id"] == proj]["tcr_cdr3"])
            overlap_set = epi_set & proj_set
            k = len(overlap_set)

            overlap_dup = df_locus.loc[
                (df_locus["tcr_project_id"] == proj)
                & (df_locus["tcr_cdr3"].isin(overlap_set)),
                "weighted_count",
            ].sum()

            x_val = pstats["dup_norm"] * dup_epi
            y_val = overlap_dup

            M, n, N = M_total, n_epi, pstats["clones_unique"]
            p_raw = 1.0 - hypergeom.cdf(k - 1, M, n, N)

            rows.append(
                (
                    proj,
                    epitope,
                    epi_stats.loc[epitope, "antigen_species"],
                    epi_stats.loc[epitope, "antigen_gene"],
                    pstats["dup_norm"],
                    dup_epi,
                    k,
                    overlap_dup,
                    x_val,
                    y_val,
                    p_raw,
                )
            )

    res = pd.DataFrame(
        rows,
        columns=[
            "project_id",
            "epitope",
            "antigen_species",
            "antigen_gene",
            "dup_project_norm",
            "dup_epitope",
            "overlap_unique",
            "overlap_dup",
            "x",
            "y",
            "p_raw",
        ],
    )

    model = LinearRegression(fit_intercept=False)
    model.fit(res[["x"]], res["y"])
    beta = model.coef_[0]

    res["y_pred"] = beta * res["x"]
    res["residual"] = np.log10((res["y"] + 1e-9) / (res["y_pred"] + 1e-9))
    res["p_adj"] = multipletests(res["p_raw"], method="fdr_bh")[1]
    res["project_id"] = res["project_id"].str.replace("TCGA-", "", regex=False)

    if plot:
        sns.set(style="whitegrid", font_scale=1.1)
        plt.figure(figsize=(10, 10))

        sns.scatterplot(
            data=res,
            x="x",
            y="y",
            hue="project_id",
            palette="tab20",
            s=60,
            alpha=0.75,
            edgecolor="none",
            legend=False,
        )

        xx = np.linspace(res["x"].min(), res["x"].max(), 200)
        plt.plot(xx, beta * xx, "--", color="black", label=f"Fit β={beta:.2e}")

        top_pos = res.sort_values("residual", ascending=False).head(top_labels)
        top_neg = res[res["y"] > 0].sort_values("residual").head(top_labels)

        for _, r in top_pos.iterrows():
            plt.annotate(
                f"{r['project_id']} | {r['antigen_gene']}",
                (r["x"], r["y"]),
                textcoords="offset points",
                xytext=(6, 6),
                fontsize=8,
                weight="bold",
                color="darkgreen",
            )

        for _, r in top_neg.iterrows():
            plt.annotate(
                f"{r['project_id']} | {r['antigen_gene']}",
                (r["x"], r["y"]),
                textcoords="offset points",
                xytext=(-30, -8),
                fontsize=8,
                weight="bold",
                color="darkred",
            )

        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("Expected overlap")
        plt.ylabel("Observed duplicate overlap")
        plt.title(f"{locus_label}: project × epitope overlap")
        plt.tight_layout()
        plt.show()

    return res
