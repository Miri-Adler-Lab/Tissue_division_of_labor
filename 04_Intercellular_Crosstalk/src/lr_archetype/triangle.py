import marimo

__generated_with = "0.23.5"
app = marimo.App(width="medium")

with app.setup:
    # Initialization code that runs before all other cells
    import marimo as mo
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.path import Path
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors

    ## setting up the triangle
    categories = {
        'endothelial': [f'AE{e}' for e in range(1, 7)],
        'fibroblast': [f'AF{f}' for f in range(1,5)],
        'macrophage': [f'AM{m}' for m in range(1,6)]
    }

    def get_colormap(interactions, log_color):
        # ----COLORMAP-----
        # Get all strengths for normalization
        if not log_color:
            all_strengths = [s for _lig, _rec, s in interactions]
        else:
            all_strengths = [np.log1p(s) for _lig, _rec, s in interactions]

        # Normalize strength values to [0, 1] for colormap
        if not all_strengths:
            norm = mcolors.Normalize(vmin=0, vmax=1)
        else:
            norm = mcolors.Normalize(vmin=min(all_strengths), vmax=max(all_strengths))
        cmap = cm.get_cmap("YlGnBu")  # choose any: 'viridis', 'coolwarm', etc.
        return norm, cmap

    def outward_normal(start, end, triangle_center):
        vec = end - start
        orthogonal = np.array([-vec[1], vec[0]])
        orth_unit = orthogonal / np.linalg.norm(orthogonal)
        midpoint = (start + end) / 2
        to_center = triangle_center - midpoint
        if np.dot(orth_unit, to_center) > 0:
            orth_unit *= -1  # Flip if it's pointing inward
        return orth_unit


    def trig(
        interactions,
        scale=3.0,
        marker_radius=0.12,
        markersize=40,
        fontsize=15,
        arrow_mutation=10,
        lw_scaler=1,
        log_line_width=False,
        log_color=False,
        colors="YlGnBu",
        categories=categories,
        title=None,
        logd=False,
        catnames=False,
        **kwargs,
    ):
        """
        TODO add dockstring

        """
        if log_color:
            logd = True

        # for reverse lookup
        instance_to_category = {
            instance: cat
            for cat, instances in categories.items()
            for instance in instances
        }
        corners = np.array([[0, 0], [1, 0], [0.5, np.sqrt(3) / 2]])
        category_edges = {
            list(categories.keys())[0]: (corners[0], corners[1]),
            list(categories.keys())[1]: (corners[1], corners[2]),
            list(categories.keys())[2]: (corners[2], corners[0]),
        }
        category_colors = {
            list(categories.keys())[0]: "blue",
            list(categories.keys())[1]: "green",
            list(categories.keys())[2]: "red",
        }

        norm, cmap = get_colormap(interactions, log_color)

        # ---- SETUP ----
        node_positions = {}
        for cat, nodes in categories.items():
            start, end = category_edges[cat]
            for i, node in enumerate(nodes):
                t = (i + 1) / (len(nodes) + 1)
                pos = start * (1 - t) + end * t
                node_positions[node] = pos

        fig, ax = plt.subplots(figsize=(10, 8))

        # Scale node positions outward from center
        center = corners.mean(axis=0)
        scaled_positions = {
            k: center + (np.array(v) - center) * scale
            for k, v in node_positions.items()
        }
        scaled_corners = np.array(
            [center + (pt - center) * scale for pt in corners]
        )
        category_edges_scaled = {
            cat: (scaled_corners[i], scaled_corners[j])
            for cat, (i, j) in {
                "endothelial": (0, 1),
                "fibroblast": (1, 2),
                "macrophage": (2, 0),
            }.items()
        }

        # ---- DRAW TRIANGLE EDGES ----
        for cat, (start, end) in category_edges_scaled.items():
            ax.plot(
                [start[0], end[0]],
                [start[1], end[1]],
                color=category_colors[cat],
                lw=2,
                alpha=0.2,
            )



        # ---- DRAW INTERACTIONS ----
        for src, tgt, strength in interactions:

            # get position and category
            src_pos = np.array(scaled_positions[src])
            tgt_pos = np.array(scaled_positions[tgt])
            src_cat = instance_to_category.get(src)
            tgt_cat = instance_to_category.get(tgt)
            if not log_line_width:
                line_strength = strength
            else:
                line_strength = np.log2(strength) if strength != 0 else strength
            # line_strength=np.log1p(strength) if log_line_width else strength

            # ---DIFFERANT_CELL_TYPES----
            if src_cat != tgt_cat:
                # Arrow across triangle
                direction = tgt_pos - src_pos
                unit_vec = direction / np.linalg.norm(direction)
                start = src_pos + unit_vec * marker_radius
                end = tgt_pos - unit_vec * marker_radius

                ax.annotate(
                    "",
                    xy=end,
                    xycoords="data",
                    xytext=start,
                    textcoords="data",
                    arrowprops=dict(
                        arrowstyle="-|>",
                        lw=lw_scaler * line_strength,
                        color=cmap(norm(strength)),
                        connectionstyle="arc3,rad=0.1",
                        mutation_scale=arrow_mutation,
                        shrinkA=0,
                        shrinkB=0,
                    ),
                )
            # ---SAME_ARCHETYPE----
            elif src == tgt:
                # Self loop using cubic Bezier
                cat_edge = category_edges_scaled.get(src_cat)
                direction = cat_edge[1] - cat_edge[0]
                angle = np.degrees(np.arctan2(direction[1], direction[0]))
                off_vec = outward_normal(cat_edge[0], cat_edge[1], center)

                start = src_pos + marker_radius * np.array(
                    [
                        np.cos(np.radians(angle - 150)),
                        np.sin(np.radians(angle - 150)),
                    ]
                )
                end = tgt_pos + marker_radius * np.array(
                    [
                        np.cos(np.radians(angle - 30)),
                        np.sin(np.radians(angle - 30)),
                    ]
                )

                control1 = (
                    start
                    + off_vec * 1.8 * marker_radius
                    - 0.08 * direction / np.linalg.norm(direction)
                )
                control2 = (
                    end
                    + off_vec * 1.8 * marker_radius
                    + 0.08 * direction / np.linalg.norm(direction)
                )

                path = Path(
                    [start, control1, control2, end],
                    [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4],
                )
                arrow = FancyArrowPatch(
                    path=path,
                    arrowstyle="-|>",
                    # color='black',
                    color=cmap(norm(strength)),
                    lw=lw_scaler * line_strength,
                    mutation_scale=arrow_mutation,
                )
                ax.add_patch(arrow)

            # ---SAME_CELL_TYPE----
            else:
                # Within-category curved arrow
                cat_edge = category_edges_scaled.get(src_cat)
                off_vec = outward_normal(cat_edge[0], cat_edge[1], center)
                start = src_pos + off_vec * marker_radius
                end = tgt_pos + off_vec * marker_radius
                mid = (start + end) / 2
                twin_control = (
                    1.2
                    if np.linalg.norm(cat_edge[0] - start)
                    > np.linalg.norm(cat_edge[0] - end)
                    else 0.8
                )

                control = mid + 3.5 * marker_radius * off_vec * twin_control
                path = Path(
                    [start, control, end], [Path.MOVETO, Path.CURVE3, Path.CURVE3]
                )
                arrow = FancyArrowPatch(
                    path=path,
                    arrowstyle="-|>",
                    # color='black',
                    color=cmap(norm(strength)),
                    lw=lw_scaler * line_strength,
                    mutation_scale=arrow_mutation,
                )
                ax.add_patch(arrow)
        # ---- DRAW NODES ----
        for node, pos in scaled_positions.items():
            ax.plot(
                pos[0],
                pos[1],
                "o",
                color=category_colors[instance_to_category[node]],
                markersize=markersize,
                zorder=10,
            )
            ax.text(
                pos[0],
                pos[1],
                node,
                fontsize=fontsize,
                ha="center",
                va="center",
                color="white",
                fontweight="bold",
                zorder=11,
            )

        # ----COLORBAR----
        sm = cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])  # dummy array
        cbar = plt.colorbar(sm, ax=ax, fraction=0.03, pad=0.04)
        if not logd:
            cbar.set_label("lr-pairs", fontsize=fontsize)
        else:
            cbar.set_label(r"$log_2(lr-pairs)$", fontsize=fontsize)
        cbar.ax.tick_params(labelsize=fontsize - 2)

        if catnames:
        # -----EDGE_LABLES------
            for cat, (p1, p2) in category_edges_scaled.items():
                mid = (p1 + p2) / 2
                direction = p2 - p1
                ortho = np.array([direction[1], -direction[0]])
                ortho = ortho / np.linalg.norm(ortho)
                offset = 0.2 * scale
                label_pos = mid + ortho * offset
                ang = np.degrees(np.arctan2(direction[1], direction[0]))
                if cat != "endothelial":
                    ang -= 180
                ax.text(
                    label_pos[0],
                    label_pos[1],
                    cat,
                    fontsize= 2 * fontsize,
                    ha="center",
                    va="center",
                    color=category_colors[cat],
                    fontweight="bold",
                    rotation=ang,
                )
        if title!=None:
            fig.suptitle(title, fontsize=27)
        # ---- FINALIZE PLOT ----
        ax.set_aspect("equal")
        ax.axis("off")
        plt.tight_layout()
        return fig


@app.cell(hide_code=True)
def _():
    mo.md(r"""
    # testing
    """)
    return


@app.cell
def _():
    #temporary data till i get it working
    interactions = [
        ('AF1', 'AM3', 1.0),
        ('AM3', 'AF1', 0.5),
        ('AM2', 'AF4', 2.0),
        ('AF4', 'AF4', 2),
        ('AE2', 'AM5', 0.8),
        # ('AF3', 'AF5', 2.1),
        ('AM2', 'AM2', 2),
        ('AE1', 'AE1', 2),
        ('AM1', 'AM3', 3),
        ('AM3', 'AM1', 2)
    ]
    return (interactions,)


@app.cell
def _(interactions):
    trig(interactions, title="testing trig function", catnames=True)
    return


if __name__ == "__main__":
    app.run()
