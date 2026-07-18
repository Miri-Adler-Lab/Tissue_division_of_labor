import marimo

__generated_with = "0.23.5"
app = marimo.App(
    width="columns",
    layout_file="layouts/all_interactions_update_presentation.slides.json",
)


@app.cell(column=0)
def _():
    import sys
    import os
    import pathlib
    # Add the src directory to sys.path so lr_archetype can be imported
    # current_dir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in locals() else os.getcwd()
    # if "Code" in current_dir:
    #     project_root_path = pathlib.Path(current_dir).parents[1]
    # else:
    #     project_root_path = pathlib.Path(current_dir)

    # src_path = os.path.join(current_dir, "src")
    # if src_path not in sys.path:
    #     sys.path.append(src_path)

    # data_dir = project_root_path / "Data" / "04_Crosstalk"
    # results_dir = project_root_path / "Results" / "Figure_4_Crosstalk"
    # results_dir.mkdir(parents=True, exist_ok=True)

    import marimo as mo
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    from scipy import stats
    from statsmodels.stats.multicomp import pairwise_tukeyhsd
    from lr_archetype.triangle import trig
    from lr_archetype.interactions import interactions

    import matplotlib.patches as patches

    import networkx as nx

    return interactions, mo, np, nx, pairwise_tukeyhsd, pd, plt, sns, stats


@app.cell
def _(np, pd):
    # importing interaction df
    lr_df = pd.read_csv("new_data/archetype_lr_pairs.csv")
    # Replace inf with nan
    lr_df.replace([np.inf, -np.inf], np.nan, inplace=True)
    # Fill nan in pair_weight with 0
    lr_df['pair_weight'] = lr_df['pair_weight'].fillna(0)
    lr_df
    return (lr_df,)


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
    return


@app.cell
def _(pd):
    tissues = (
        pd.read_csv(
            "new_data/archetypes_tissue_correlation.csv", index_col="index"
        )
        .stack()
        .to_dict()
    )
    tissues
    return (tissues,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## all interactions
    """)
    return


@app.cell
def _(interactions, lr_df, tissues):
    # saveing the data for fig 6 b-c
    # this is old data TODO change it
    interactions(lr_df).tissue_correction(tissues).condense().df.to_csv("results/fig4_b_data.csv")
    #this was 1
    return


@app.cell
def _(interactions, lr_df, tissues):
    interactions(lr_df).tissue_correction(tissues).condense().df
    return


@app.cell
def _(interactions, lr_df, tissues):
    interactions(lr_df).tissue_correction(tissues, cutoff=0.1).min_mean_diff(0.5).df
    return


@app.cell
def _(interactions, lr_df, tissues):
    interactions(lr_df).tissue_correction(tissues, cutoff=0.1).min_mean_diff(0.5).df.to_csv("results/fig4_b_pairs_data.csv")
    return


@app.cell
def _(interactions, lr_df, mo, tissues):
    #old fig 6 b 
    #fig_6_b = interactions(lr_df).tissue_correction(tissues).condence().set_min(8).log().plot(lw_scaler=0.7)
    fig_6_b = interactions(lr_df).tissue_correction(tissues, cutoff=0.2).min_mean_diff(0.3).condense().log().plot(lw_scaler=1, catnames=True)
    fig_6_b.savefig("results/fig4_b.svg", format='svg')
    mo.vstack([mo.md("#all interactions"),fig_6_b])
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Endothelial
    """)
    return


@app.cell
def _(interactions, lr_df, mo, tissues):
    endo_l = interactions(lr_df).tissue_correction(tissues,cutoff=0.1).min_mean_diff(0.3).condense().slice_cell_type({'ligands':['endothelial'], 'receptors':['endothelial', 'fibroblast', 'macrophage']}).log().plot(lw_scaler=0.7)
    endo_ln= interactions(lr_df).tissue_correction(tissues, cutoff=0.1).min_mean_diff(0.3).condense().slice_cell_type({'ligands':['endothelial'], 'receptors':['fibroblast', 'macrophage']}).log().plot(lw_scaler=0.7)

    mo.vstack([mo.md("# Endothealial - ligand"),mo.hstack([endo_l, endo_ln])])
    return


@app.cell
def _(interactions, lr_df, mo, tissues):
    endo_r = interactions(lr_df).tissue_correction(tissues, 0.1).min_mean_diff(0.3).condense().slice_cell_type({'receptors':['endothelial'], 'ligands':['endothelial', 'fibroblast', 'macrophage']}).log().plot(lw_scaler=0.7)
    endo_rn = interactions(lr_df).tissue_correction(tissues, 0.1).min_mean_diff(0.3).condense().slice_cell_type({'receptors':['endothelial'], 'ligands':['fibroblast', 'macrophage']}).log().plot(lw_scaler=0.7)

    mo.vstack([mo.md("# Endothealial - receptor"),mo.hstack([endo_r, endo_rn])])
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Fibroblasts
    """)
    return


@app.cell
def _(interactions, lr_df, mo, tissues):
    fib_l = interactions(lr_df).tissue_correction(tissues, 0.1).min_mean_diff(0.3).condense().log().slice_cell_type({'ligands':['fibroblast'], 'receptors':['endothelial', 'fibroblast', 'macrophage']}).plot(lw_scaler=0.7)
    fib_ln = interactions(lr_df).tissue_correction(tissues, 0.1).min_mean_diff(0.3).condense().log().slice_cell_type({'ligands':['fibroblast'], 'receptors':['endothelial','macrophage']}).plot(lw_scaler=0.7)
    mo.vstack([mo.md("# Fibroblasts - ligand"),mo.hstack([fib_l, fib_ln])])
    return (fib_l,)


@app.cell
def _(fib_l):
    fib_l.savefig("new_data/fig_c_fib_l.svg", format="svg")
    return


@app.cell
def _(interactions, lr_df, mo, tissues):
    fib_r = interactions(lr_df).tissue_correction(tissues,0.1).min_mean_diff(0.3).condense().log().slice_cell_type({'receptors':['fibroblast'], 'ligands':['endothelial', 'fibroblast', 'macrophage']}).plot(lw_scaler=0.7)
    fib_rn = interactions(lr_df).tissue_correction(tissues, 0.1).min_mean_diff(0.3).condense().log().slice_cell_type({'receptors':['fibroblast'], 'ligands':['endothelial', 'macrophage']}).plot(lw_scaler=0.7)
    mo.vstack([mo.md("# Fibroblasts - receptor"),mo.hstack([fib_r, fib_rn])])
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Macrophages
    """)
    return


@app.cell
def _(interactions, lr_df, mo, tissues):
    mac_l = interactions(lr_df).tissue_correction(tissues, 0.1).min_mean_diff(0.3).condense().log().slice_cell_type({'ligands':['macrophage'], 'receptors':['endothelial', 'fibroblast', 'macrophage']}).plot(lw_scaler=0.7)
    mac_ln = interactions(lr_df).tissue_correction(tissues,0.1).min_mean_diff(0.3).condense().log().slice_cell_type({'ligands':['macrophage'], 'receptors':['endothelial', 'fibroblast']}).plot(lw_scaler=0.7)
    mo.vstack([mo.md("# Macrophage - ligand"),mo.hstack([mac_l, mac_ln])])
    return


@app.cell
def _(interactions, lr_df, mo, tissues):
    mac_r = interactions(lr_df).tissue_correction(tissues, 0.1).min_mean_diff(0.3).condense().log().slice_cell_type({'receptors':['macrophage'], 'ligands':['endothelial', 'fibroblast', 'macrophage']}).plot(lw_scaler=0.7)
    mac_rn= interactions(lr_df).tissue_correction(tissues, 0.1).min_mean_diff(0.3).condense().log().slice_cell_type({'receptors':['macrophage'], 'ligands':['endothelial', 'fibroblast']}).plot(lw_scaler=0.7)
    mo.vstack([mo.md("# Macrophage - receptor"),mo.hstack([mac_r, mac_rn])])
    return


@app.cell
def _(mo):
    mo.md(r"""
    # Autocrine vs Paracrine
    """)
    return


@app.cell
def _(interactions, lr_df, tissues):
    auto_para_df = interactions(lr_df).tissue_correction(tissues).min_mean_diff(0.3).condense().df
    auto_para_df["Autocrine"] = auto_para_df['cell_type_ligand'] == auto_para_df['cell_type_receptor']
    auto_para_df["auto_auto"] = auto_para_df.index.get_level_values(0) == auto_para_df.index.get_level_values(1)
    # auto_para_df
    return (auto_para_df,)


@app.cell
def _(interactions, lr_df, tissues):
    ligand_by_archetype = (
        interactions(lr_df)
        .tissue_correction(tissues)
        .condense()
        .df.groupby(by=["cell_type_ligand","archetype_ligand"])["pair_weight"]
        .sum()
        .to_frame(name="ligand_sum")
    )
    ligand_by_archetype
    return


@app.cell
def _(auto_para_df):
    ap_df = auto_para_df.groupby(["cell_type_ligand", "Autocrine"]).agg({"tiss_corl":"sum"}).reset_index()
    # ap_df
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
    # a_df
    return a_df, p_df, size_dict


@app.cell
def _(auto_para_df, size_dict):
    aa_df = auto_para_df.groupby(["cell_type_ligand", "auto_auto"]).agg({"tiss_corl":"sum"}).query("auto_auto == True").reset_index()
    aa_df['normed_counts'] = aa_df['tiss_corl'] / ((aa_df['cell_type_ligand'].map(size_dict) - 1) )
    # aa_df
    return (aa_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    graph of the mean number of lr pairs in each category, auto auto, autocrine (including auto_auto), and parancrine.

    nomelization

    - auto auto = counts / (num(cell type) - 1)
    - autocrine = counts / ((num(cell type) - 1) ** 2)
    - paracrine = counts / (num(other cell types) * num(cell type))


    as expected the highest value in each cell type is the auto-autocrine, with differant ratios for each cell type, endothelial has the largest ratio while the others have a much smaller ratio.

    the ratio of autocrine and paracrine is varsly differant in the differant cell types, while endothelial and fibroblast have a lower autocrine level macropheges have a much higher autocrine level
    """)
    return


@app.cell
def _(a_df, aa_df, np, p_df, plt):
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
    # auto_para_fig.savefig(results_dir / "autocrine_paracrine.svg", format='svg')
    auto_para_fig
    return


@app.cell
def _(a_df, aa_df, p_df, pd):
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
    auto_para_table.to_csv("new_data/autocrine_paracrine.csv")
    auto_para_table
    return


@app.cell(column=1)
def _(lr_df):
    lr_df.groupby(["enriched_gene_ligand", "enriched_gene_receptor"]).count()
    return


@app.cell
def _(agged):
    agged.to_csv("results/fig4_d_data.csv")
    return


@app.cell
def _(lr_df):
    agged = lr_df.groupby(["cell_type_ligand", "cell_type_receptor"]).agg({"pair_weight":"sum"}).reset_index()
    edges_with_weights = list(zip(
        agged["cell_type_ligand"],
        agged["cell_type_receptor"],
        agged["pair_weight"]
    ))
    edges_with_weights
    return agged, edges_with_weights


@app.cell
def _(edges_with_weights, nx, plt):
    # 2. Create a Directed Graph
    G = nx.DiGraph()

    # 3. Add edges with their weight attributes
    for source, target, weight in edges_with_weights:
        G.add_edge(source, target, weight=weight)

    # 4. Set up the layout geometry (Circular keeps nodes separated and readable)
    pos = nx.circular_layout(G)

    # 5. Define aesthetic features
    node_colors = {
        "endothelial": "#4682B4",  # Steel Blue
        "fibroblast": "#E67E22",   # Orange
        "macrophage": "#2ECC71",   # Emerald Green
    }
    colors = [node_colors[_node] for _node in G.nodes()]

    # Extract weights for scaling arrow thicknesses dynamically
    weights = [G[u][v]["weight"] for u, v in G.edges()]
    # Normalize weights so they don't get comically thick or invisible
    max_weight = max(weights)
    edge_widths = [(w / max_weight) * 5 for w in weights]

    # 6. Plot the network
    fig, ax = plt.subplots(figsize=(7, 7))

    # Draw the cell type nodes
    nx.draw_networkx_nodes(
        G, pos, 
        node_size=2500, 
        node_color=colors, 
        alpha=0.9, 
        ax=ax
    )

    # Draw the node text labels
    nx.draw_networkx_labels(
        G, pos, 
        font_size=11, 
        font_weight="bold", 
        font_family="sans-serif", 
        ax=ax
    )

    # Draw the directional signaling arrows
    nx.draw_networkx_edges(
        G, pos,
        width=edge_widths,
        edge_color="#7F8C8D",
        arrowstyle="->",
        arrowsize=25,
        connectionstyle="arc3,rad=0.15",  # Curving the lines lets you see back-and-forth pairs clearly
        ax=ax
    )

    # 7. Polish the canvas
    ax.set_title("Cell-Type Ligand-Receptor Communication Network", fontsize=14, fontweight="bold", pad=20)
    plt.axis("off")
    plt.tight_layout()
    fig.savefig("results/fig4_d.svg", format="svg")
    fig
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


@app.cell
def _(cell_type_summery, sns):
    sns.violinplot(data=cell_type_summery, x="cell_type", y="lr_ratio")
    return


@app.cell(column=2)
def _():
    return


@app.cell
def _(interactions, lr_df, pd, tissues):
    ligand_by_cell_type = (
        interactions(lr_df)
        .tissue_correction(tissues)
        .condense()
        .df
        .groupby(by=["cell_type_ligand","archetype_ligand"])["pair_weight"]
        .sum()
        .to_frame(name="ligand_sum")
    )
    ligand_by_cell_type.reset_index(inplace=True)
    ligand_by_cell_type.columns = ["cell_type", "archetype", "ligand_sum"]
    ligand_by_cell_type.set_index = ['cell_type', 'archetype']
    # ligand_by_cell_type

    receptor_by_cell_type = (
        interactions(lr_df)
        .tissue_correction(tissues)
        .condense()
        .df.groupby(by=["cell_type_receptor","archetype_receptor"])["pair_weight"]
        .sum()
        .to_frame(name="receptor_sum")
    )
    receptor_by_cell_type.reset_index(inplace=True)
    receptor_by_cell_type.columns = ["cell_type", "archetype", "receptor_sum"]
    receptor_by_cell_type.set_index = ['cell_type', 'archetype']

    cell_type_summery = pd.merge(ligand_by_cell_type, receptor_by_cell_type, how='inner', on=['cell_type', 'archetype'])
    nn = {"endothelial": 6, "fibroblast": 5, "macrophage": 5}
    cell_type_summery["lr_ratio"] = cell_type_summery['ligand_sum'] / cell_type_summery['receptor_sum']
    cell_type_summery
    return (cell_type_summery,)


@app.cell
def _(cell_type_summery, sns):
    sns.scatterplot(data=cell_type_summery, x="archetype", y="receptor_sum", hue="cell_type")
    return


@app.cell
def _(cell_type_summery, sns):
    sns.scatterplot(data=cell_type_summery, x="archetype", y="ligand_sum", hue="cell_type")
    return


@app.cell
def _(cell_type_summery, sns):
    sns.scatterplot(data=cell_type_summery, x="archetype", y="lr_ratio", hue="cell_type")
    return


@app.cell
def _(cell_type_summery, ci_lower, ci_upper, np, plt, sample_mean, sns):
    cell_type_summery["Ligand + Receptor Sum"] = (
        cell_type_summery["ligand_sum"] + cell_type_summery["receptor_sum"]
    )
    _fig, _ax = plt.subplots()
    sns.scatterplot( x=cell_type_summery["archetype"], y=np.log2(cell_type_summery["lr_ratio"]), hue=cell_type_summery["cell_type"], size=cell_type_summery["Ligand + Receptor Sum"], sizes=(10, 250),ax=_ax)

    # _ax.axhline(y=0, color="red", label="0 -> equal ligand receptor", ls="-.", linewidth=1)
    _ax.axhline(y=sample_mean, color="purple", ls="--", label="mean lr ratio", alpha=0.3)
    _ax.axhspan(ci_lower, ci_upper, color="lightgray", alpha=0.4, zorder=0, label="95%ci")


    _ax.hlines(y=cell_type_summery.loc[cell_type_summery["cell_type"] == "endothelial", "log_lr"].mean(), xmin=0, xmax=5, ls=":", color="blue", label="mean endothelaial")
    _ax.hlines(y=cell_type_summery.loc[cell_type_summery["cell_type"] == "macrophage", "log_lr"].mean(), xmin=10, xmax=14, ls=":", color="green", label="mean macrophage")
    _ax.hlines(y=cell_type_summery.loc[cell_type_summery["cell_type"] == "fibroblast", "log_lr"].mean(), xmin=5.5, xmax=10, ls=":", color="orange", label="mean fibroblast")

    _ax.set_ylabel(r"$log_2 (\frac{ligand}{receptor})$")
    _ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    # _ax.grid(False)

    _ax.grid(True, axis='x', color='gainsboro', linestyle=':', linewidth=1)

    _fig.savefig('results/fig4_c.svg', format="svg")
    cell_type_summery.to_csv("results/fig4_c_data.csv")
    _ax
    return


@app.cell
def _(cell_type_summery, np, plt, sns):
    _fig, _ax = plt.subplots()
    sns.barplot(x=cell_type_summery["archetype"], y=np.log2(cell_type_summery["lr_ratio"]), hue=cell_type_summery["Ligand + Receptor Sum"], ax=_ax)
    _ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., title="Weighted sum \n ligand + receptor pairs")
    _ax.set_ylabel(r"$log_2 (\frac{ligand}{receptor})$")
    _ax.set_title("autocrine")
    _fig
    return


@app.cell
def _(cell_type_summery, np, sns):
    lr_ratios = np.log2(cell_type_summery["lr_ratio"])
    cell_type_summery["log_lr"] =  lr_ratios
    np.mean(lr_ratios)
    sns.kdeplot(lr_ratios, fill=True)
    return


@app.cell
def _(cell_type_summery, np, plt, sns, stats):
    # Generate some random, normally distributed data
    np.random.seed(42)
    # data = np.random.normal(loc=20, scale=5, size=100)
    data = cell_type_summery["log_lr"]

    # Set the seaborn style for a clean look
    sns.set_theme(style="whitegrid")

    # Create the Q-Q plot
    _fig, _ax = plt.subplots(figsize=(6, 6))
    stats.probplot(data, dist="norm", plot=_ax)

    # Customize titles/labels using matplotlib via seaborn styling
    _ax.set_title("Q-Q Plot for Normal Distribution")
    return (data,)


@app.cell
def _(data, np, stats):
    n = len(data)
    df = n - 1

    # Calculate sample mean and standard error of the mean (SEM)
    sample_mean = np.mean(data)
    sem = stats.sem(data)  # This computes s / sqrt(n)

    # Calculate 95% Confidence Interval
    confidence_level = 0.95
    ci_lower, ci_upper = stats.t.interval(confidence_level, df, loc=sample_mean, scale=sem)

    print(f"Sample Mean: {sample_mean:.2f}")
    print(f"95% Confidence Interval: ({ci_lower:.2f}, {ci_upper:.2f})")
    return ci_lower, ci_upper, sample_mean


@app.cell
def _():
    #todo - make the cell type in one graph colord by cell type and y ligand receptor
    #todo - make a graph of 3 nodes one for each cell type, with autocrine and parancrine
    return


@app.cell
def _(cell_type_summery, stats):
    endothelial_ligand = cell_type_summery[cell_type_summery['cell_type'] == "endothelial"]["ligand_sum"]
    fibroblast_ligand = cell_type_summery[cell_type_summery['cell_type'] == "fibroblast"]["ligand_sum"]
    macrophage_ligand = cell_type_summery[cell_type_summery['cell_type'] == "macrophage"]["ligand_sum"]

    f_stat_ligand, p_value_ligand = stats.f_oneway(endothelial_ligand, fibroblast_ligand, macrophage_ligand)

    print(f"F-Statistic: {f_stat_ligand:.4f}")
    print(f"P-Value: {p_value_ligand:.4f}")
    return


@app.cell
def _(cell_type_summery, stats):
    endothelial_receptor = cell_type_summery[cell_type_summery['cell_type'] == "endothelial"]["receptor_sum"]
    fibroblast_receptor = cell_type_summery[cell_type_summery['cell_type'] == "fibroblast"]["receptor_sum"]
    macrophage_receptor = cell_type_summery[cell_type_summery['cell_type'] == "macrophage"]["receptor_sum"]

    f_stat_receptor, p_value_receptor = stats.f_oneway(endothelial_receptor, fibroblast_receptor, macrophage_receptor)

    print(f"F-Statistic: {f_stat_receptor:.4f}")
    print(f"P-Value: {p_value_receptor:.4f}")
    return


@app.cell
def _(cell_type_summery, pairwise_tukeyhsd):
    _tukey = pairwise_tukeyhsd(endog=cell_type_summery['ligand_sum'], groups=cell_type_summery['cell_type'], alpha=0.05)

    print(_tukey)
    return


@app.cell
def _(cell_type_summery, pairwise_tukeyhsd):
    # Perform Tukey HSD
    _tukey = pairwise_tukeyhsd(endog=cell_type_summery['receptor_sum'], groups=cell_type_summery['cell_type'], alpha=0.05)

    print(_tukey)
    return


@app.cell
def _(cell_type_summery, plt, sns):
    _fig, _ax = plt.subplots(2, figsize=(7,10))


    sns.boxplot(x=cell_type_summery.cell_type, y=cell_type_summery.ligand_sum, ax=_ax[0])
    sns.boxplot(x=cell_type_summery.cell_type, y=cell_type_summery.receptor_sum, ax=_ax[1])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    for now boxplotted the ligands and receptors of eache cell type, using anova tuki-kremer no significant differacne between cell types.

    noticable archetypes outliers are AM1 with a low amount of ligands, AM4 with high amount of receptors, and AE1 with a low amount of receptors.
    """)
    return


@app.cell
def _(cell_type_summery, sns):
    sns.boxplot(data=cell_type_summery, x="cell_type", y="lr_ratio")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    looking at the ratios of the ligand receptor we see a two outliers, AF3 extreemly high and AM5 also high
    """)
    return


@app.cell
def _(cell_type_summery, plt, sns):
    _fig, _ax = plt.subplots(1,3, figsize=(15,3))
    sns.boxplot(cell_type_summery.ligand_sum, ax=_ax[0])
    sns.boxplot(cell_type_summery.receptor_sum, ax=_ax[1])
    sns.boxplot(cell_type_summery.lr_ratio, ax=_ax[2])
    return


@app.cell
def _(cell_type_summery, sns):
    sns.histplot(cell_type_summery.receptor_sum,binwidth=20)
    return


@app.cell
def _(cell_type_summery, sns):
    sns.histplot(cell_type_summery.lr_ratio,binwidth=0.5)
    return


@app.cell
def _(cell_type_summery, plt, sns):
    pdf_melted = cell_type_summery.melt(
        id_vars=["archetype", "cell_type"], 
        value_vars=["ligand_sum", "receptor_sum"],
        var_name="molecule_type", 
        value_name="expression_sum"
    )

    # Now plot in one command
    _fig, _ax = plt.subplots()
    sns.scatterplot(
        data=pdf_melted, 
        x="archetype", 
        y="expression_sum", 
        hue="cell_type", 
        style="molecule_type", # This gives Ligands and Receptors different shapes automatically
        s=100
    )
    _ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    return


@app.cell(column=3)
def _(cell_type_summery, interactions, lr_df, pd, tissues):
    ligand_by_cell_type_p = (
        interactions(lr_df)
        .tissue_correction(tissues)
        .condense()
        .df
        # .groupby(by=["cell_type_ligand","archetype_ligand"])["pair_weight"]
        # .sum()
        # .to_frame(name="ligand_sum")
    )
    ligand_by_cell_type_p.reset_index(inplace=True)

    ligand_by_cell_type_p = ligand_by_cell_type_p[ligand_by_cell_type_p["archetype_ligand"] != ligand_by_cell_type_p["archetype_receptor"]]
    ligand_by_cell_type_p = ligand_by_cell_type_p.groupby(by=["cell_type_ligand","archetype_ligand"])["pair_weight"].sum().to_frame(name="ligand_sum")
    ligand_by_cell_type_p.reset_index(inplace=True)
    ligand_by_cell_type_p.columns = ["cell_type", "archetype", "ligand_sum"]
    ligand_by_cell_type_p.set_index = ['cell_type', 'archetype']

    receptor_by_cell_type_p = (
        interactions(lr_df)
        .tissue_correction(tissues)
        .condense()
        .df
        # .groupby(by=["cell_type_receptor","archetype_receptor"])["pair_weight"]
        # .sum()
        # .to_frame(name="receptor_sum")
    )
    receptor_by_cell_type_p.reset_index(inplace=True)

    receptor_by_cell_type_p = receptor_by_cell_type_p[receptor_by_cell_type_p["archetype_ligand"] != receptor_by_cell_type_p["archetype_receptor"]]
    receptor_by_cell_type_p = receptor_by_cell_type_p.groupby(by=["cell_type_receptor","archetype_receptor"])["pair_weight"].sum().to_frame(name="receptor_sum")
    receptor_by_cell_type_p.reset_index(inplace=True)
    receptor_by_cell_type_p.columns = ["cell_type", "archetype", "receptor_sum"]
    receptor_by_cell_type_p.set_index = ['cell_type', 'archetype']

    cell_type_summery_p = pd.merge(ligand_by_cell_type_p, receptor_by_cell_type_p, how='inner', on=['cell_type', 'archetype'])
    # nn = {"endothelial": 6, "fibroblast": 5, "macrophage": 5}
    cell_type_summery_p["lr_ratio"] = cell_type_summery['ligand_sum'] / cell_type_summery['receptor_sum']
    cell_type_summery_p["Ligand + Receptor Sum"] = (
        cell_type_summery_p["ligand_sum"] + cell_type_summery_p["receptor_sum"]
    )
    cell_type_summery_p
    return (cell_type_summery_p,)


@app.cell
def _(cell_type_summery_p, np, plt, sns):
    _fig, _ax = plt.subplots()
    sns.barplot(x=cell_type_summery_p["archetype"], y=np.log2(cell_type_summery_p["lr_ratio"]), hue=cell_type_summery_p["Ligand + Receptor Sum"], ax=_ax)
    _ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., title="Weighted sum \n ligand + receptor pairs")
    _ax.set_ylabel(r"$log_2 (\frac{ligand}{receptor})$")
    _ax.set_title("no autocrine")
    no_autocrine = _ax
    _fig
    return


@app.cell
def _(cell_type_summery_p, np, sns):
    sns.scatterplot(x=cell_type_summery_p["archetype"], y=np.log2(cell_type_summery_p["receptor_sum"]),hue=cell_type_summery_p["cell_type"])
    return


@app.cell
def _(cell_type_summery_p, sns):
    sns.scatterplot(data=cell_type_summery_p, x="archetype", y="receptor_sum", hue="cell_type")
    return


@app.cell
def _(cell_type_summery_p, sns):
    sns.scatterplot(data=cell_type_summery_p, x="archetype", y="ligand_sum", hue="cell_type")
    return


@app.cell
def _(cell_type_summery_p, sns):
    sns.scatterplot(data=cell_type_summery_p, x="cell_type", y="lr_ratio", hue="cell_type")
    return


@app.cell
def _(cell_type_summery_p, plt, sns):
    _fig, _ax = plt.subplots()
    sns.scatterplot(data=cell_type_summery_p, x="archetype", y="ligand_sum", hue="cell_type", marker="s")
    sns.scatterplot(data=cell_type_summery_p, x="archetype", y="receptor_sum", hue="cell_type")
    _ax.legend()
    return


@app.cell
def _(cell_type_summery_p, plt, sns):
    df_melted = cell_type_summery_p.melt(
        id_vars=["archetype", "cell_type"], 
        value_vars=["ligand_sum", "receptor_sum"],
        var_name="molecule_type", 
        value_name="expression_sum"
    )

    # Now plot in one command
    _fig, _ax = plt.subplots()
    sns.scatterplot(
        data=df_melted, 
        x="archetype", 
        y="expression_sum", 
        hue="cell_type", 
        style="molecule_type", # This gives Ligands and Receptors different shapes automatically
        s=100
    )
    _ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    return


@app.cell(column=4)
def _(cell_type_summery, cell_type_summery_p, np, plt, sns):
    _fig, _ax = plt.subplots(1,2, figsize=(30,10))
    sns.barplot(x=cell_type_summery["archetype"], y=np.log2(cell_type_summery["lr_ratio"]), hue=cell_type_summery["Ligand + Receptor Sum"], ax=_ax[0])
    _ax[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., title="Weighted sum \n ligand + receptor pairs")
    _ax[0].set_ylabel(r"$log_2 (\frac{ligand}{receptor})$")
    _ax[0].set_title("autocrine")

    sns.barplot(x=cell_type_summery_p["archetype"], y=np.log2(cell_type_summery_p["lr_ratio"]), hue=cell_type_summery_p["Ligand + Receptor Sum"], ax=_ax[1])
    _ax[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0., title="Weighted sum \n ligand + receptor pairs")
    _ax[1].set_ylabel(r"$log_2 (\frac{ligand}{receptor})$")
    _ax[1].set_title("no autocrine")
    _fig
    return


if __name__ == "__main__":
    app.run()
