import marimo

__generated_with = "0.23.5"
app = marimo.App(width="medium")

with app.setup:
    # Initialization code that runs before all other cells
    import marimo as mo
    import pandas as pd
    import numpy as np
    from lr_archetype.triangle import trig


@app.class_definition
class interactions:
    def __init__(self, df, logd=False):
        self.df = df.copy(deep=True)
        self.logd = logd

    ### all logic needs to end with modifing self.df = modification
    def __repr__(self):
        return repr(self.df)

    def tissue_correction(self, tissues, cutoff=0.1):
        """
        multiply the pair_weight by the tissue correlation.

        Args:
            tissues: Dictionary mapping (ligand_archetype, receptor_archetype) to correlation values.
            cutoff: cutoff will map anything under it to zero and above to one, set to None for full range with no mapping.
        """
        self.df["tiss_corl"] = self.df.apply(
            lambda row: tissues[
                (row["archetype_ligand"], row["archetype_receptor"])
            ],
            axis=1,
        )
        if cutoff is not None:
            self.df["tiss_corl"] = self.df["tiss_corl"].map(
                lambda x: 1 if x >= cutoff else 0
            )
        else:
            self.df["tiss_corl"] *= 1
        self.df["pair_weight"] = self.df["pair_weight"] * self.df["tiss_corl"]
        self.df["pair_weight"] = self.df["pair_weight"].fillna(0)
        return self

    def min_mean_diff(self, min):
        """sets minimum mea diff"""
        self.df = self.df.query(
            "mean_diff_ligand >= @min and mean_diff_receptor >= @min"
        )
        return self

    def min_median_diff(self, min):
        """sets minimum mea diff"""
        self.df = self.df.query(
            "median_diff_ligand >= @min and median_diff_receptor >= @min"
        )
        return self

    def set_min(self, min):
        """sets minimum strength"""
        self.df = self.df.query("tiss_corl >= @min")
        return self

    def set_max(self, max):
        """sets maximum strength"""
        self.df = self.df.query("tiss_corl <= @max")
        return self

    def slice_cell_type(self, cell_types):
        """
        selects cell types.

        Args:
            cell_types: dict['lignads':[list of cell types], 'receptors': [list of cell types]]
        """
        ligands = cell_types['ligands']
        receptors = cell_types['receptors']
        self.df = self.df.query(
            "cell_type_ligand in @ligands and cell_type_receptor in @receptors"
        )
        return self

    def slice_archetypes(self, archetypes):
        """
        selects archetypes.

        Args:
            archetypes: dict['ligands':[list of archetypes], 'receptors':[list of receptors]]
        """
        ligands = archetypes['ligands']
        receptors = archetypes['receptors']
        self.df = self.df.query(
            "archetype_ligand in @ligands and archetype_receptor in @receptors"
        )
        return self

    def slice_by_gene_list(self, gene_list):
        """
        selects interactions by gene list.
        !!! must come before condense !!!

        Args:
            gene_list: dict['ligands':[list of genes], 'receptors':[list of genes]]
        """
        ligands = gene_list['ligands']
        receptors = gene_list['receptors']
        self.df = self.df.query(
            "enriched_gene_ligand in @ligands and enriched_gene_receptor in @receptors"
        )
        return self

    def condense(self):
        # to condense all the pairs
        """
        condenses df from single interactions to archetype pair interactions.
        must be used before plot.
        """
        self.df = (
            self.df.loc[
                :,
                [
                    "archetype_ligand",
                    "cell_type_ligand",
                    "archetype_receptor",
                    "cell_type_receptor",
                    "pair_weight",
                    "tiss_corl",
                ],
            ]
            .groupby(["archetype_ligand", "archetype_receptor"])
            .agg(
                {
                    "pair_weight": "sum",
                    "cell_type_ligand": "first",
                    "cell_type_receptor": "first",
                    "tiss_corl": "sum",
                }
            )
        )
        return self

    def log(self):
        """log2 of pair weight"""
        self.df["pair_weight"] = self.df["pair_weight"].map(
            lambda x: np.log2(x) if x != 0 else x
        )
        self.df["tiss_corl"] = self.df["tiss_corl"].map(
            lambda y: np.log2(y) if y != 0 else y
        )
        self.logd = True
        return self

    def to_trig_list(self, pair_weight=False):
        # to internaly for making trig list
        if not pair_weight:
            trig_list = [
                (i, j, v) for (i, j), v in self.df["tiss_corl"].items()
            ]
        else:
            trig_list = [
                (i, j, v) for (i, j), v in self.df["pair_weight"].items()
            ]
        return trig_list

    def plot(self, pair_weight=False, **kwargs):
        trig_list = self.to_trig_list(pair_weight)
        return trig(trig_list, logd=self.logd, **kwargs)


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # testing
    """)
    return


@app.cell
def _():
    pd.read_csv('data/processed/archetype_lr_pairs.csv')
    return


if __name__ == "__main__":
    app.run()
