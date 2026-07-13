import marimo



__generated_with = "0.14.17"
app = marimo.App(width="medium")


@app.cell
def _():
    import sys
    import os
    import pathlib
    input_root = pathlib.Path(os.environ["PUBLICATION_INPUT_ROOT"]).expanduser().resolve()
    output_root = pathlib.Path(os.environ["PUBLICATION_OUTPUT_ROOT"]).expanduser().resolve()
    current_dir = pathlib.Path(__file__).resolve().parent

    src_path = os.path.join(current_dir, "src")
    if src_path not in sys.path:
        sys.path.append(src_path)
    
    data_dir = input_root / "04_Crosstalk"
    results_dir = output_root / "Figure_4_Crosstalk"
    results_dir.mkdir(parents=True, exist_ok=True)

    required_inputs = [
        data_dir / "archetype_lr_pairs.csv",
        data_dir / "archetypes_tissue_correlation.csv",
    ]
    missing = [str(path) for path in required_inputs if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing Figure 4 input(s): {missing}")

    import marimo as mo
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns
    from lr_archetype.triangle import trig
    from lr_archetype.interactions import interactions
    return data_dir, interactions, mo, np, pd, plt, results_dir, sns


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # All archetype interactions

    looking numericly - statisticaly at all interactions
    """
    )
    return


@app.cell
def _(mo):
    mo.md(r"""## importing data and making interaction table""")
    return


@app.cell
def _(mo):
    mo.md(r"""importing ligand receptor archetpe pair table""")
    return


@app.cell
def _(data_dir, np, pd):
    # importing interaction df
    lr_df = pd.read_csv(data_dir / "archetype_lr_pairs.csv")
    # Replace inf with nan
    lr_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    # Fill nan in pair_weight with 0
    lr_df['pair_weight'] = lr_df['pair_weight'].fillna(0)
    lr_df
    return (lr_df,)


@app.cell
def _(mo):
    mo.md(r"""makeing a pairwise interaction table""")
    return


@app.cell
def _(lr_df):
    lr_pairs = (
        lr_df.groupby(by=["archetype_ligand", "archetype_receptor"])
        .agg({"enriched_gene_ligand": "count", "pair_weight": "sum"})
        .reset_index()
        .pivot(
            index="archetype_ligand",
            columns="archetype_receptor",
            values="pair_weight",
        )
    ).fillna(0)
    lr_pairs
    return (lr_pairs,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""importing tissue correlation data""")
    return


@app.cell
def _(data_dir, pd):
    tissue_correlation = pd.read_csv(
        data_dir / "archetypes_tissue_correlation.csv", index_col=0
    )
    tissue_correlation
    return (tissue_correlation,)


@app.cell
def _(results_dir, tissue_correlation):
    tissue_correlation.stack().to_csv(
        results_dir / "stacked_tissue_correlation.csv"
    )
    return


@app.cell
def _(mo):
    mo.md(r"""multipying the tissue and pairs to get correlation map""")
    return


@app.cell
def _(lr_pairs, tissue_correlation):
    lra_tissue = lr_pairs * tissue_correlation
    lra_tissue
    return (lra_tissue,)


@app.cell
def _(mo):
    mo.md(r"""# Heatmaps""")
    return


@app.cell
def _(lra_tissue, plt, sns):
    def _():
        fig = plt.figure()
        sns.heatmap(lra_tissue, cmap="viridis")
        plt.title("Number of ligand receptor pairs")
        plt.xlabel("receptor")
        plt.ylabel("ligand")
        separators = [6, 11]
        for sep in separators:
            # Horizontal line
            plt.axhline(sep, color="red", linewidth=2)
            # Vertical line
            plt.axvline(sep, color="red", linewidth=2)
        return fig


    _()
    return


@app.cell
def _(lra_tissue, np, pd, plt, sns):
    # add directional graph
    mask = pd.DataFrame(np.triu(np.ones((16, 16), dtype=int)), 
                        index=lra_tissue.index, 
                        columns=lra_tissue.columns)
    mask = mask.mask(mask == 0)
    directional = (lra_tissue - lra_tissue.T) * mask
    def _():
        fig = plt.figure()
        sns.heatmap(directional, cmap="seismic", center=0)
        plt.title("counts tissue corrected")
        plt.xlabel("receptor")
        plt.ylabel("ligand")
        separators = [6, 11]
        for sep in separators:
            # Horizontal line
            plt.axhline(sep, color='red', linewidth=2)
            # Vertical line
            plt.axvline(sep, color='red', linewidth=2)
        return fig
    _()
    return


@app.cell
def _(mo):
    mo.md(r"""# Triangle Graph""")
    return


@app.cell
def _(data_dir, pd):
    tissues = (
        pd.read_csv(
            data_dir / "archetypes_tissue_correlation.csv", index_col="index"
        )
        .stack()
        .to_dict()
    )
    return (tissues,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""## all interactions""")
    return


@app.cell
def _(interactions, lr_df, results_dir, tissues):
    # saveing the data for fig 6 b-c
    # this is old data TODO change it
    interactions(lr_df).tissue_correction(tissues).condense().df.to_csv(results_dir / "fig_b-c_data.csv")
    return


@app.cell
def _(interactions, lr_df, results_dir, tissues):
    interactions(lr_df).tissue_correction(tissues, cutoff=0.3).min_mean_diff(0.5).df.to_csv(results_dir / "fig_b-c_lr_pairs_data.csv")
    return


@app.cell
def _(interactions, lr_df, results_dir, tissues):
    #old fig 6 b 
    #fig_6_b = interactions(lr_df).tissue_correction(tissues).condence().set_min(8).log().plot(lw_scaler=0.7)
    fig_6_b = interactions(lr_df).tissue_correction(tissues, cutoff=0.3).min_mean_diff(0.5).condense().log().plot(lw_scaler=1, catnames=True)
    fig_6_b.savefig(results_dir / "fig_b_map_all.svg", format='svg')
    fig_6_b
    return


@app.cell
def _(mo):
    mo.md(r"""# Endothelial""")
    return


@app.cell
def _(interactions, lr_df, tissues):
    endo_l = interactions(lr_df).tissue_correction(tissues).min_mean_diff(0.3).condense().slice_cell_type({'ligands':['endothelial'], 'receptors':['endothelial', 'fibroblast', 'macrophage']}).log().plot(lw_scaler=0.7)
    endo_l
    return (endo_l,)


@app.cell
def _(endo_l, results_dir):
    endo_l.savefig(results_dir / "fig_c_endo_l.svg", format='svg')
    return


@app.cell
def _(interactions, lr_df, tissues):
    endo_r = interactions(lr_df).tissue_correction(tissues).min_mean_diff(0.3).condense().slice_cell_type({'receptors':['endothelial'], 'ligands':['endothelial', 'fibroblast', 'macrophage']}).log().plot(lw_scaler=0.7)
    endo_r
    return (endo_r,)


@app.cell
def _(endo_r, results_dir):
    endo_r.savefig(results_dir / "fig_c_endo_r.svg", format="svg")
    return


@app.cell
def _(mo):
    mo.md(r"""# Fibroblasts""")
    return


@app.cell
def _(interactions, lr_df, tissues):
    fib_l = interactions(lr_df).tissue_correction(tissues).min_mean_diff(0.3).condense().log().slice_cell_type({'ligands':['fibroblast'], 'receptors':['endothelial', 'fibroblast', 'macrophage']}).plot(lw_scaler=0.7)
    fib_l
    return (fib_l,)


@app.cell
def _(fib_l, results_dir):
    fib_l.savefig(results_dir / "fig_c_fib_l.svg", format="svg")
    return


@app.cell
def _(interactions, lr_df, tissues):
    fib_r = interactions(lr_df).tissue_correction(tissues).min_mean_diff(0.3).condense().log().slice_cell_type({'receptors':['fibroblast'], 'ligands':['endothelial', 'fibroblast', 'macrophage']}).plot(lw_scaler=0.7)
    fib_r
    return (fib_r,)


@app.cell
def _(fib_r, results_dir):
    fib_r.savefig(results_dir / "fig_c_fib_r.svg", format="svg")
    return


@app.cell
def _(mo):
    mo.md(r"""# Macrophages""")
    return


@app.cell
def _(interactions, lr_df, tissues):
    mac_l = interactions(lr_df).tissue_correction(tissues).min_mean_diff(0.3).condense().log().slice_cell_type({'ligands':['macrophage'], 'receptors':['endothelial', 'fibroblast', 'macrophage']}).plot(lw_scaler=0.7)
    mac_l
    return (mac_l,)


@app.cell
def _(mac_l, results_dir):
    mac_l.savefig(results_dir / "fig_c_mac_l.svg", format="svg")
    return


@app.cell
def _(interactions, lr_df, tissues):
    mac_r = interactions(lr_df).tissue_correction(tissues).min_mean_diff(0.3).condense().log().slice_cell_type({'receptors':['macrophage'], 'ligands':['endothelial', 'fibroblast', 'macrophage']}).plot(lw_scaler=0.7)
    mac_r
    return (mac_r,)


@app.cell
def _(mac_r, results_dir):
    mac_r.savefig(results_dir / "fig_c_mac_r.svg", format="svg")
    return


@app.cell
def _(interactions, lr_df, pd, tissues):
    ligand_by_cell_type = (
        interactions(lr_df)
        .tissue_correction(tissues)
        .condense()
        .df.groupby(by="cell_type_ligand")
        .sum("pair_weight")
        .rename(index={"cell_type_ligand": "cell_type"}, columns={"pair_weight": "ligand_sum"})
    )
    ligand_by_cell_type.index.name = 'cell_type'
    receptor_by_cell_type = (
        interactions(lr_df)
        .tissue_correction(tissues)
        .condense()
        .df.groupby(by="cell_type_receptor")
        .sum("pair_weight")
        .rename(columns={'pair_weight':'receptor_sum'})
    )
    receptor_by_cell_type.index.name = 'cell_type'
    cell_type_summery = pd.merge(ligand_by_cell_type, receptor_by_cell_type, how='inner', on='cell_type')
    cell_type_summery
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    # fig 6 a 

    finding a pair of ligand receptor that are unike to one pair of archetypes
    """
    )
    return


@app.cell
def _(interactions, lr_df, tissues):
    lig_count = interactions(lr_df).tissue_correction(tissues).df.groupby("enriched_gene_ligand").agg({"archetype_ligand":"nunique", "archetype_receptor":"nunique"})

    lig_count['total_count'] = lig_count["archetype_ligand"] + lig_count["archetype_receptor"]
    genes_of_interest = lig_count.query('total_count == 2').index.tolist()

    # Step 3: Filter original df based on those genes
    interactions(lr_df).tissue_correction(tissues).df.query("enriched_gene_ligand in @genes_of_interest and pair_weight != 0 and cell_type_ligand != cell_type_receptor")#[:,["archetype_ligand", "enriched_gene_ligand",]]
    return


@app.cell
def _(mo):
    mo.md(r"""# Autocrine vs Paracrine""")
    return


@app.cell
def _(interactions, lr_df, tissues):
    auto_para_df = interactions(lr_df).tissue_correction(tissues).min_mean_diff(0.3).condense().df
    auto_para_df["Autocrine"] = auto_para_df['cell_type_ligand'] == auto_para_df['cell_type_receptor']
    auto_para_df["auto_auto"] = auto_para_df.index.get_level_values(0) == auto_para_df.index.get_level_values(1)
    auto_para_df
    return (auto_para_df,)


@app.cell
def _(auto_para_df):
    ap_df = auto_para_df.groupby(["cell_type_ligand", "Autocrine"]).agg({"tiss_corl":"sum"}).reset_index()
    ap_df
    return (ap_df,)


@app.cell
def _(ap_df):
    size_dict = {
        'endothelial': 6,
        'fibroblast': 5,
        'macrophage': 5
    }

    def normlize_ap_df(row):
        if row["Autocrine"]:
            return row["tiss_corl"] / (size_dict[row['cell_type_ligand']] - 1) ** 2
        else:
            others = list(size_dict.keys())  # call keys()
            others.remove(row['cell_type_ligand'])  # remove current key
            return row['tiss_corl'] / ((size_dict[others[0]] + size_dict[others[1]]) * size_dict[row['cell_type_ligand']])

        # if row['Autocrine']:
        #     return row["count"]
        # else:
        #     return row['count'] / 2ap_df['normed_counts'] = ap_df.apply(normlize_ap_df, axis=1)


    ap_df['normed_counts'] = ap_df.apply(normlize_ap_df, axis=1)
    ap_df
    a_df = ap_df.query("Autocrine == True")
    p_df = ap_df.query("Autocrine == False")
    a_df
    return a_df, p_df, size_dict


@app.cell
def _(auto_para_df, size_dict):
    aa_df = auto_para_df.groupby(["cell_type_ligand", "auto_auto"]).agg({"tiss_corl":"sum"}).query("auto_auto == True").reset_index()
    aa_df['normed_counts'] = aa_df['tiss_corl'] / ((aa_df['cell_type_ligand'].map(size_dict) - 1) )
    aa_df
    return (aa_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
    graph of the mean number of lr pairs in each category, auto auto, autocrine (including auto_auto), and parancrine.

    nomelization

    - auto auto = counts / (num(cell type) - 1)
    - autocrine = counts / ((num(cell type) - 1) ** 2)
    - paracrine = counts / (num(other cell types) * num(cell type))
  

    as expected the highest value in each cell type is the auto-autocrine, with differant ratios for each cell type, endothelial has the largest ratio while the others have a much smaller ratio.

    the ratio of autocrine and paracrine is varsly differant in the differant cell types, while endothelial and fibroblast have a lower autocrine level macropheges have a much higher autocrine level
    """
    )
    return


@app.cell
def _(a_df, aa_df, np, p_df, plt, results_dir):
    def plot_grouped_bars():
        labels = a_df.cell_type_ligand.tolist()
        x = np.arange(len(labels))  # label locations

        width = 0.25  # smaller bar width

        fig, ax = plt.subplots()

        ax.bar(x - width, aa_df.normed_counts, width, label='Auto-Autocrine')
        ax.bar(x, a_df.normed_counts, width, label='Autocrine')
        ax.bar(x + width, p_df.normed_counts, width, label='Paracrine')

        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend()
        ax.set_ylabel("mean interactions")
        return fig

    auto_para_fig = plot_grouped_bars()
    auto_para_fig.savefig(results_dir / "autocrine_paracrine.svg", format='svg')
    auto_para_fig


    return


@app.cell
def _(a_df, aa_df, p_df, pd, results_dir):
    auto_para_table = pd.DataFrame(
        index=['endothelial', 'fibroblast', 'macrophage'],
        # columns=['auto_autocrine', 'autocrine', 'parancrine']
    )
    auto_para_table['auto_autocrine'] = auto_para_table.index.map(
        aa_df.set_index('cell_type_ligand')['normed_counts']
    )
    auto_para_table['autocrine'] = auto_para_table.index.map(
        a_df.set_index('cell_type_ligand')['normed_counts']
    )
    auto_para_table['paracrine'] = auto_para_table.index.map(
        p_df.set_index('cell_type_ligand')['normed_counts']
    )
    auto_para_table.to_csv(results_dir / "autocrine_paracrine.csv")
    auto_para_table
    return


@app.cell
def _(aa_df):
    aa_df.normed_counts
    return


@app.cell
def _(a_df, aa_df, size_dict):
    aaa_df = a_df.copy().reset_index()
    aaa_df['ap_counts'] = aaa_df['tiss_corl'] - aa_df['tiss_corl']
    aaa_df['ap_normed'] = aaa_df['ap_counts'] / (aaa_df['cell_type_ligand'].map(size_dict) - 1 )
    aaa_df['ratio'] = aaa_df['ap_normed'] / aaa_df['normed_counts']
    aaa_df
    return


if __name__ == "__main__":
    app.run()
