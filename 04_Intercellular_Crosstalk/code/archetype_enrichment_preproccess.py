import marimo

__generated_with = "0.14.13"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    return mo, pd


@app.cell
def _():
    ENDOTHELIAL_FILE = (
        "data/raw/enrichment/ParTi_Lite_Out_Endo_6_continuous_significant.csv"
    )
    FIBROBLAST_FILE = "data/raw/enrichment/fibroblast_enriched_genes.csv"
    MACROPHAGE_FILE = "data/raw/enrichment/Mac_5_gene_enrichment.csv"
    return ENDOTHELIAL_FILE, FIBROBLAST_FILE, MACROPHAGE_FILE


@app.cell
def _(pd):
    def read_enriched_genes(file, name, cell_type):
        df = (
            pd.read_csv(file)
            .rename(
                columns={
                    "archetype #": "archetype",
                    "Feature Name": "enriched_gene",
                    "P value (Mann-Whitney)": "p_val",
                    "Median Difference": "median_diff",
                    "Mean Difference": "mean_diff",
                }
            )
            .drop(
                columns=[
                    "Significant after Benjamini-Hochberg correction?",
                    "Is first bin maximal?",
                ]
            )
        )
        df["archetype"] = (name + df["archetype"].astype(str)).astype("category")
        df['cell_type'] = cell_type
        return df
    return (read_enriched_genes,)


@app.cell
def _(
    ENDOTHELIAL_FILE,
    FIBROBLAST_FILE,
    MACROPHAGE_FILE,
    pd,
    read_enriched_genes,
):
    arc_gene_df = pd.concat(
        [
            read_enriched_genes(
                ENDOTHELIAL_FILE,
                "AE",
                'endothelial'
            ),
            read_enriched_genes(
                FIBROBLAST_FILE,
                "AF",
                'fibroblast'
            ),
            read_enriched_genes(
                MACROPHAGE_FILE,
                "AM",
                'macrophage'
            ),
        ]
    ).sort_index()
    arc_gene_df["archetype"] = arc_gene_df["archetype"].astype("category")
    arc_gene_df = arc_gene_df.reset_index(drop=True)
    return (arc_gene_df,)


@app.cell
def _(arc_gene_df, mo):
    mo.ui.table(arc_gene_df)
    return


@app.cell
def _(arc_gene_df):
    #save
    arc_gene_df.to_csv('data/processed/all_enrichment.csv', index=False)
    return


if __name__ == "__main__":
    app.run()
