import marimo

__generated_with = "0.14.13"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    import json
    return np, pd, plt, sns


@app.cell
def _():
    ENDOTHELIAL_FILE = 'data/raw/tissues/endo_archetype_summary_thresh_0.30.csv'
    FIBROBLAST_FILE = 'data/raw/tissues/fib_archetype_summary_thresh_0.30.csv'
    MACROPHAGE_FILE = 'data/raw/tissues/mac_archetype_summary_thresh_0.30.csv'
    return ENDOTHELIAL_FILE, FIBROBLAST_FILE, MACROPHAGE_FILE


@app.cell
def _(pd):
    def tiss_enrich(name, file):
        df = pd.read_csv(file)
        df["Adjusted_Dominant"] = (
            name + df["Adjusted_Dominant"].astype(str)
        ).astype("category")
        grouped_dicts = df.groupby("Adjusted_Dominant", observed=False).apply(
            lambda g: dict(zip(g["Tissue"], g["Adjusted_Count"]))
        )
        return grouped_dicts
    return (tiss_enrich,)


@app.cell
def _(ENDOTHELIAL_FILE, tiss_enrich):
    endothelial_tissue_enrichment = tiss_enrich("AE", ENDOTHELIAL_FILE)
    return (endothelial_tissue_enrichment,)


@app.cell
def _(FIBROBLAST_FILE, tiss_enrich):
    fibroblast_tissue_enrichment = tiss_enrich("AF", FIBROBLAST_FILE)
    return (fibroblast_tissue_enrichment,)


@app.cell
def _(MACROPHAGE_FILE, tiss_enrich):
    macrophage_tissue_enrichment = tiss_enrich("AM", MACROPHAGE_FILE)
    return (macrophage_tissue_enrichment,)


@app.cell
def _(
    endothelial_tissue_enrichment,
    fibroblast_tissue_enrichment,
    macrophage_tissue_enrichment,
    pd,
):
    tiss_df = pd.concat(
        [
            endothelial_tissue_enrichment,
            fibroblast_tissue_enrichment,
            macrophage_tissue_enrichment
        ]
    ).sort_index()
    unwanted = ["AF0", "AE0", "AM0"]
    tiss_df = tiss_df[~tiss_df.index.isin(unwanted)]
    tiss_df.index.name = 'index'
    tiss_df
    return (tiss_df,)


@app.cell
def _(tiss_df):
    tiss_df.index
    return


@app.cell
def _(np, pd, plt, sns, tiss_df):
    # First compute dot products
    archetype_dicts = tiss_df  # assuming tiss_df is a Series with dicts as values
    archetypes = archetype_dicts.index
    similarity_df = pd.DataFrame(index=archetypes, columns=archetypes, dtype=float)
    for i in archetypes:
        for j in archetypes:
            dict_i = archetype_dicts[i]
            dict_j = archetype_dicts[j]
            shared_keys = set(dict_i) & set(dict_j)
            dot_product = sum(dict_i[k] * dict_j[k] for k in shared_keys)
            similarity_df.at[i, j] = dot_product

    # Now normalize using cosine similarity
    for i in archetypes:
        for j in archetypes:
            norm_i = np.sqrt(sum(v**2 for v in archetype_dicts[i].values()))
            norm_j = np.sqrt(sum(v**2 for v in archetype_dicts[j].values()))
            if norm_i > 0 and norm_j > 0:
                similarity_df.at[i, j] /= (norm_i * norm_j)
            else:
                similarity_df.at[i, j] = 0

    # # similarity_df
    sns.heatmap(similarity_df, cmap="viridis")
    plt.title("cosine distance map")
    return (similarity_df,)


@app.cell
def _(similarity_df):
    similarity_df
    return


@app.cell
def _(similarity_df):
    similarity_df.stack().to_dict()

    return


@app.cell
def _(similarity_df):
    # saves the similarity_df to csv
    similarity_df.to_csv("data/processed/archetypes_tissue_correlation.csv")
    return


if __name__ == "__main__":
    app.run()
