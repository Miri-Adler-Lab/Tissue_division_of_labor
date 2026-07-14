import marimo

__generated_with = "0.23.5"
app = marimo.App(width="medium")

with app.setup:
    # Initialization code that runs before all other cells
    import marimo as mo
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import json
    import gseapy as gp
    import warnings
    warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')
    go_bp = gp.get_library(name="GO_Biological_Process_2025", organism="Human")
    # import pathlib
    # import os
    # import sys

    # # Robustly find project root from script location
    # current_dir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in locals() else os.getcwd()
    # if "Code" in current_dir:
    #     project_root_path = pathlib.Path(current_dir).parents[1]
    # else:
    #     project_root_path = pathlib.Path(current_dir)

    # Add src to path if needed (assuming structure similar to all_interactions)
    # Add src to path if needed (assuming structure similar to all_interactions)
    # src_path = os.path.join(current_dir, "src")
    # if src_path not in sys.path:
    #     sys.path.append(src_path)

    # data_dir = project_root_path / "Data" / "04_Crosstalk"
    # local_data_dir = pathlib.Path(current_dir) / "new_data"
    # results_dir = project_root_path / "Results" / "Figure_4_Crosstalk"
    # results_dir.mkdir(parents=True, exist_ok=True)

    from lr_archetype.triangle import trig
    from lr_archetype.interactions import interactions


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # Pathway Analisis

    ## General idea
    mapping ligand receptor interactions related to different processes and pathways trying to see if they are specific to particular archetypes or generalists.
    specific archetypes will appear as distinctive arrows between a small amount of archetypes, but generalists can appear in two ways, either as a lot of small arrows between a wide range of archetypes without a hierarchy or order, or not appear at all making it more generalist thus no ligand or receptor is expressed significantly different between archetypes, and therefore is not in the analysis.
    """)
    return


@app.cell
def _():
    mo.md(r"""
    # Importing data
    ## interactions
    """)
    return


@app.cell
def _():
    lr_pairs = pd.read_csv( "new_data/archetype_lr_pairs.csv")
    # lr_pairs
    return (lr_pairs,)


@app.cell
def _():
    # imports tissue correlation as a dictionary
    tissues = (
        pd.read_csv(
            "new_data/archetypes_tissue_correlation.csv", index_col="index"
        )
        .stack()
        .to_dict()
    )
    return (tissues,)


@app.function
def get_lr_from_go(
    go_terms,
):
    """
    gets a list of ligand receptors that are related to the GO terms
    go_terms: list of full name go terms
    tip- make a multiselect of go_bp.keys()
    """
    # get a dict of lists for ligands and receptors
    with open("new_data/ligand_receptor_lists.json") as _f:
        lr_lists = json.load(_f)

    #### add go list functionality
    related_genes = list({_gene for go_term in go_terms for _gene in go_bp[go_term]})
    process_lr = dict()
    # filters the lr lists for related genes
    for key, gene_list in lr_lists.items():
        process_lr[key] = [
            _gene for _gene in gene_list if _gene in related_genes
        ]
    for k, gl in process_lr.items():
        process_lr[k] = list(set(process_lr[k]))

    return process_lr


@app.function
def get_lr_from_gene_list(
    gene_list,
):
    """
    for lr_dict from gene list if it is not in the go term finder
    """
    # get a dict of lists for ligands and receptors
    with open("new_data/ligand_receptor_lists.json") as _f:
        lr_lists = json.load(_f)

    #### add go list functionality
    process_lr = dict()
    # filters the lr lists for related genes
    for key, _list in lr_lists.items():
        process_temp = [
            _gene for _gene in _list if _gene in gene_list
        ]
        process_temp = list(set(process_temp)) if process_temp else lr_lists[key] # makes list and fills if one is empty
        process_lr[key] = process_temp

    if process_lr == lr_lists:
        raise ValueError('no genes fit')

    return process_lr


@app.cell
def _():
    mo.md(r"""
    # Interactive mapper
    """)
    return


@app.cell
def _():
    go_term_selector = mo.ui.multiselect(go_bp.keys())
    go_term_selector
    return (go_term_selector,)


@app.cell
def _(go_term_selector):
    go_term_list = go_term_selector.value
    go_term_list
    return (go_term_list,)


@app.cell
def _(go_term_list, lr_pairs, tissues):
    try:
        mo.output.append(interactions(lr_pairs).tissue_correction(tissues).slice_by_gene_list(get_lr_from_go(go_term_list)).condense().plot(log_line_width=False, lw=3))
    except Exception as e:
        print(f"Exception caught: {type(e).__name__}: {e}")
        if not go_term_list:
            mo.output.append(mo.callout("please select go terms from the list", "neutral"))
        else:
            mo.output.append(mo.callout("The GO terms selected have no ligand receptor pairs prominant in the data set", "danger"))
    return


@app.cell
def _():
    mo.md(r"""
    # specific processes
    """)
    return


@app.cell
def _(lr_pairs, tissues):
    go_list = go_bp.keys()

    def map_from_key_term(key_term):
        key_term_go_list = [
            _key for _key in go_list if key_term.lower() in _key.lower()
        ]
        for _go in key_term_go_list:
            try:
                mo.output.append(
                    interactions(lr_pairs)
                    .tissue_correction(tissues)
                    .slice_by_gene_list(get_lr_from_go([_go]))
                    .condense()
                    .plot(title=_go)
                )
                print(_go)
            except Exception as e:
                pass
        mo.output.append(
            interactions(lr_pairs)
            .tissue_correction(tissues)
            .slice_by_gene_list(get_lr_from_go(key_term_go_list))
            .condense()
            .plot(title="all terms")
        )

    return


@app.cell
def _():
    mo.md(r"""
    ## angiogenesis
    """)
    return


@app.cell
def _():
    # map_from_key_term('angiogenesis')
    return


@app.cell
def _():
    mo.md(r"""
    ### fig 6 D
    """)
    return


@app.cell
def _(lr_pairs, tissues):
    reg_of_angiogenesis = interactions(lr_pairs).tissue_correction(tissues).min_mean_diff(0.3).slice_by_gene_list(get_lr_from_go(["Regulation of Angiogenesis (GO:0045765)"]))
    reg_of_angiogenesis.df.to_csv("results/fig4_f_data.csv")
    reg_of_angiogenesis_plt = reg_of_angiogenesis.condense().plot(title="Regulation of Angiogenesis (GO:0045765)", lw_scaler=1)
    reg_of_angiogenesis_plt.savefig("results/fig4_f.svg", format='svg')
    reg_of_angiogenesis_plt
    return


@app.cell
def _():
    mo.md(r"""
    ## ECM
    """)
    return


@app.cell
def _():
    # map_from_key_term('extracellular matrix')
    return


@app.cell
def _(lr_pairs, tissues):
    ecm_organization = interactions(lr_pairs).tissue_correction(tissues).min_mean_diff(0.3).slice_by_gene_list(get_lr_from_go(["Extracellular Matrix Organization (GO:0030198)"]))
    ecm_organization.df.to_csv("results/fig4_h_data.csv")
    ecm_organization_plt = ecm_organization.condense().plot(title="Extracellular Matrix Organization (GO:0030198)", lw_scaler=1)
    ecm_organization_plt.savefig("results/fig4_h.svg", format='svg')
    ecm_organization_plt
    return


@app.cell
def _():
    mo.md(r"""
    ## inflamatory
    """)
    return


@app.cell
def _():
    # map_from_key_term('inflammatory')
    return


@app.cell
def _(lr_pairs, tissues):
    reg_inflamitory = interactions(lr_pairs).tissue_correction(tissues).min_mean_diff(0.3).slice_by_gene_list(get_lr_from_go(["Regulation of Inflammatory Response (GO:0050727)"]))
    reg_inflamitory.df.to_csv("results/fig4_g_data.csv")
    reg_inflamitory_plt = reg_inflamitory.condense().plot(title="Regulation of Inflammatory Response (GO:0050727)", lw_scalar=1)
    reg_inflamitory_plt.savefig("results/fig4_g.svg", format="svg")
    reg_inflamitory_plt
    return


@app.cell
def _():
    mo.md(r"""
    ## Growth factor
    """)
    return


@app.cell
def _():
    growth_factor_go_list = list(set(pd.read_excel("new_data/1GO_term_summary_20220615_085602.xlsx")['Symbol']))
    growth_factor_go_list = [x.upper() for x in growth_factor_go_list]
    # growth_factor_go_list
    return (growth_factor_go_list,)


@app.cell
def _(growth_factor_go_list, lr_pairs, tissues):
    interactions(lr_pairs).tissue_correction(tissues).min_mean_diff(0.3).slice_by_gene_list(get_lr_from_gene_list(growth_factor_go_list)).condense().plot(title="growth factor activity (GO:0008083)")
    return


@app.cell
def _(growth_factor_go_list, lr_pairs, tissues):
    growth_factors = interactions(lr_pairs).tissue_correction(tissues).min_mean_diff(0.3).slice_by_gene_list(get_lr_from_gene_list(growth_factor_go_list))
    growth_factors.df.to_csv("results/fig4_e_data.csv")
    growth_factors_plt = growth_factors.condense().plot(title="growth factor activity (GO:0008083)", lw_scalar=1)
    growth_factors_plt.savefig("results/fig4_e.svg", format="svg")
    growth_factors_plt
    return


if __name__ == "__main__":
    app.run()
