import os
from typing import List

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests


def analyze_gene_set_vs_cancer(
    df: pd.DataFrame,
    output_csv: str = "cancer_regression.csv",
    heatmap_file: str = "gene_cancer_heatmap.png",
) -> pd.DataFrame:
    """
    Perform linear regression for each gene against cancer type (project_id) using one-hot encoding,
    and generate a heatmap of regression coefficients.

    Parameters:
        df (pd.DataFrame): Input dataframe with 'sample', 'project_id', and gene columns.
        output_csv (str): Path to save the regression result table as CSV.
        heatmap_file (str): Path to save the heatmap of regression coefficients.

    Returns:
        pd.DataFrame: DataFrame containing significant gene-cancer associations with adjusted p-values.
    """
    results = []

    project_dummies = pd.get_dummies(df["project_id"], drop_first=True).astype(int)
    merged_df = pd.concat(
        [df.drop(columns=["sample", "project_id"]), project_dummies], axis=1
    )

    gene_columns = [col for col in df.columns if col not in ["sample", "project_id"]]
    cancer_types = project_dummies.columns

    coef_matrix = pd.DataFrame(index=gene_columns, columns=cancer_types, dtype=float)

    for gene in gene_columns:
        y = np.log1p(df[gene])
        X = sm.add_constant(project_dummies)
        model = sm.OLS(y, X).fit()

        p_values = model.pvalues[1:]
        coefficients = model.params[1:]
        genes = list(X.columns[1:])

        reject, p_corrected, _, _ = multipletests(p_values, method="fdr_bh")

        for cancer_type, coef, p_corr, sig in zip(
            genes, coefficients, p_corrected, reject
        ):
            if sig and not pd.isna(p_corr):
                results.append(
                    [gene, cancer_type, coef, p_corr, -np.log10(p_corr + 1e-12)]
                )
                coef_matrix.loc[gene, cancer_type] = coef
            else:
                coef_matrix.loc[gene, cancer_type] = 0.0

    results_df = pd.DataFrame(
        results, columns=["Gene", "Cancer", "coef", "adj_p", "-log10(p)"]
    )
    results_df.to_csv(output_csv, index=False)

    coef_matrix = coef_matrix.loc[~(coef_matrix == 0).all(axis=1)]

    plt.figure(figsize=(14, max(6, len(coef_matrix) * 0.3)))
    sns.heatmap(
        coef_matrix,
        cmap="vlag",
        center=0,
        annot=True,
        fmt=".1f",
        linewidths=0.5,
        cbar_kws={"label": "coefficient"},
    )
    plt.title("Gene ~ Cancer Type Regression Coefficient Heatmap")
    plt.xlabel("Cancer type")
    plt.tight_layout()
    plt.savefig(heatmap_file, dpi=300)
    plt.show()

    return results_df


def analyze_bcr_tcr(
    bcr_df: pd.DataFrame,
    tcr_df: pd.DataFrame,
    output_csv: str = "regression_results.csv",
) -> pd.DataFrame:
    """
    Perform regression analysis to test interaction effects between TCR expression and cancer types
    on BCR expression levels.

    Parameters:
        bcr_df (pd.DataFrame): DataFrame with BCR genes and columns ['sample', 'project_id'].
        tcr_df (pd.DataFrame): DataFrame with TCR genes and columns ['sample', 'project_id'].
        output_csv (str): Path to save the resulting regression table as CSV.

    Returns:
        pd.DataFrame: DataFrame of significant interaction terms with adjusted p-values.
    """
    merged_df = pd.merge(
        bcr_df, tcr_df, on=["sample", "project_id"], suffixes=("_bcr", "_tcr")
    )
    project_dummies = pd.get_dummies(merged_df["project_id"]).astype(int)
    merged_df = pd.concat(
        [merged_df.drop(columns=["project_id"]), project_dummies], axis=1
    )

    results = []

    bcr_genes = [col for col in bcr_df.columns if col not in ["sample", "project_id"]]
    tcr_genes = [col for col in tcr_df.columns if col not in ["sample", "project_id"]]

    for bcr_gene in bcr_genes:
        for tcr_gene in tcr_genes:
            interactions = {
                f"{tcr_gene}_x_{cancer}": merged_df[tcr_gene] * merged_df[cancer]
                for cancer in project_dummies.columns
            }
            interaction_df = pd.DataFrame(interactions)

            X = sm.add_constant(interaction_df)
            y = merged_df[bcr_gene]

            model = sm.OLS(y, X).fit()
            p_values = model.pvalues[1:]
            coefficients = model.params[1:]
            terms = list(X.columns[1:])

            reject, p_corrected, _, _ = multipletests(p_values, method="fdr_bh")

            for term, coef, p_corr, sig in zip(
                terms, coefficients, p_corrected, reject
            ):
                if sig and not pd.isna(p_corr):
                    tcr_pred, cancer_pred = term.split("_x_")
                    results.append(
                        [
                            bcr_gene,
                            tcr_pred,
                            cancer_pred,
                            coef,
                            p_corr,
                            -np.log10(p_corr + 1e-12),
                        ]
                    )

    results_df = pd.DataFrame(
        results,
        columns=[
            "dependent_var",
            "gene_predictor",
            "TCGA_predictor",
            "coef",
            "adj_p",
            "-log10(p-value)",
        ],
    )
    results_df.to_csv(output_csv, index=False)

    return results_df
