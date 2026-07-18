import marimo

__generated_with = "0.14.13"
app = marimo.App(width="medium")

with app.setup:
    import marimo as mo
    import pandas as pd
    import json
    import scipy.stats as ss
    from statsmodels.stats.multitest import multipletests


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    # making a lignad receptor data set easy to use

    [ligand receptor data base](https://doi.org/10.1093/bib/bbaa269)

    slicing and converting to jason, saving in *proccessed_data*
    """
    )
    return


@app.cell
def _():
    # lr data set
    lrds = pd.read_csv(
        "data/preprocessed/ligand_receptors/human_lr_pair.txt", sep=r"\s+"
    )[["lr_pair", "ligand_gene_symbol", "receptor_gene_symbol"]]
    lrds
    return (lrds,)


@app.cell
def _(lrds):
    ligands = lrds.ligand_gene_symbol.unique()
    receptors = lrds.receptor_gene_symbol.unique()

    lr_dict = dict()
    rl_dict = dict()
    for _, row in lrds.iterrows():   # _ = index, row = Series
        ligand = row['ligand_gene_symbol']
        receptor = row['receptor_gene_symbol']

        if ligand not in lr_dict:
            lr_dict[ligand] = [receptor]
        else:
            lr_dict[ligand].append(receptor)

        if receptor not in rl_dict:
            rl_dict[receptor] = [ligand]
        else:
            rl_dict[receptor].append(ligand)

    lr_dict
    ligands
    return ligands, lr_dict, receptors


@app.cell
def _(ligands, lr_dict, receptors):
    # Save
    with open("data/processed/lrds.json", "w") as f:
        json.dump(lr_dict, f)

    with open("data/processed/ligand_receptor_lists.json", "w") as j:
        json.dump({'ligands':list(ligands), 'receptors':list(receptors)}, j)
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    # Finding the interactions

    gets archtype gene enrichment, finds all ligand receptor enriched genes, and pairs them up
    """
    )
    return


@app.cell
def _():
    enrichment_df = pd.read_csv("data/processed/all_enrichment.csv")
    mo.ui.table(enrichment_df)
    return (enrichment_df,)


@app.cell
def _():
    mo.md(r"""making two tables one of enriched ligands, and one of enriched receptors""")
    return


@app.cell
def _(enrichment_df):
    arch_ligand = enrichment_df.query('enriched_gene in @ligands')
    arch_receptor = enrichment_df.query('enriched_gene in @receptors')

    return arch_ligand, arch_receptor


@app.cell
def _():
    mo.md(
        r"""
    merging the tables where the ligand and receptor match.
    checking significance using fischers method for combining pvalues $\chi_{2k}^{2} = -2 \sum_{i=1}^{k} \text{ln} \  p_{i}$ and binyamini-hochberg for muliple comparisons

    each ligand receptor pair gets a **pair weight** that is the product of the median differance of expression of the ligand and receptor
    """
    )
    return


@app.cell
def _(arch_ligand, arch_receptor, lr_dict):
    # Step  1: Cross-merge all combinations of ligand and receptor rows
    merged = arch_ligand.merge(
        arch_receptor, how="cross", suffixes=("_ligand", "_receptor")
    )

    # Step 2: Filter to keep only valid pairs from lr_dict
    lr_df = merged[
        merged.apply(
            lambda row: row["enriched_gene_ligand"] in lr_dict
            and row["enriched_gene_receptor"]
            in lr_dict[row["enriched_gene_ligand"]],
            axis=1,
        )
    ].reset_index(drop=True)

    # pair weight
    lr_df["pair_weight"] = lr_df.median_diff_ligand * lr_df.median_diff_receptor

    # statistics
    lr_df[["fisher_stat", "combined_pval"]] = lr_df.apply(
        lambda row: pd.Series(
            ss.combine_pvalues(
                [row["p_val_ligand"], row["p_val_receptor"]], method="fisher"
            )
        ),
        axis=1,
    )
    rejected, corrected_pvals, _, _ = multipletests(
        lr_df["combined_pval"], method="fdr_bh"
    )
    lr_df["pval_bh"] = corrected_pvals
    lr_df["rejected"] = rejected

    lr_df
    return (lr_df,)


@app.cell
def _(lr_df):
    lr_df.to_csv("data/processed/archetype_lr_pairs.csv", index=False)
    return


if __name__ == "__main__":
    app.run()
