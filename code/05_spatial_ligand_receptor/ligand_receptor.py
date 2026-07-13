import os
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as mcolors
import pandas as pd
import numpy as np
from scipy import sparse
from scipy.sparse import issparse
from scipy.stats import pearsonr
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests
from tqdm import trange

# ── Constants ─────────────────────────────────────────────────────────────────
PALETTE = {
    "endothelial": "#457B9D",
    "macrophage":  "#E63946",
    "fibroblast":  "#2A9D8F",
    "Other":       "#D3D3D3",
}
CT_COLORS = {
    'fibroblast':  'blue',
    'macrophage':  'green',
    'endothelial': 'red',
}


# ── Data loading ──────────────────────────────────────────────────────────────
def broad_category(ct):
    ct_lower = ct.lower()
    if "fibroblast" in ct_lower:
        return "fibroblast"
    elif "macrophage" in ct_lower:
        return "macrophage"
    elif "endothelial" in ct_lower:
        return "endothelial"
    else:
        return "Other"


def load_data():
    adata = ad.read_h5ad("../heart_visium/cdaf82fb-7d23-4234-b871-b65b3bce91f4.h5ad")
    print(f"Spots: {adata.n_obs}")
    print("Cell types present:", adata.obs["cell_type"].unique().tolist())
    adata.obs["broad_type"] = adata.obs["cell_type"].map(broad_category).astype("category")
    print(adata.obs["broad_type"].value_counts())
    print(len(adata.obs["sample"].unique()))

    return adata


def load_lr_pairs():
    lr_pairs = pd.read_csv("../heart_visium/lr_pairs_data.csv", index_col=0)
    return lr_pairs[lr_pairs['tiss_corl'] == 1]


def load_archetype_genes(top_n=30):
    archetype_genes = {}
    endo_enrich = pd.read_csv('../heart_visium/Endothelials_gene_enrichment.csv')
    for arch_num in endo_enrich['archetype #'].unique():
        subset = endo_enrich[endo_enrich['archetype #'] == arch_num].nlargest(top_n, 'Median Difference')
        archetype_genes[f'AE{arch_num}'] = subset['Feature Name'].tolist()

    macro_enrich = pd.read_csv('../heart_visium/Macrophages_gene_enrichment.csv')
    for arch_num in macro_enrich['archetype #'].unique():
        subset = macro_enrich[macro_enrich['archetype #'] == arch_num].nlargest(top_n, 'Median Difference')
        archetype_genes[f'AM{arch_num}'] = subset['Feature Name'].tolist()

    fibro_path = '../heart_visium/Fibroblasts_gene_enrichment.csv'
    if os.path.exists(fibro_path):
        fibro_enrich = pd.read_csv(fibro_path)
        for arch_num in fibro_enrich['archetype #'].unique():
            subset = fibro_enrich[fibro_enrich['archetype #'] == arch_num].nlargest(top_n, 'Median Difference')
            archetype_genes[f'AF{arch_num}'] = subset['Feature Name'].tolist()
    return archetype_genes


# ── Spatial cell-type overview ────────────────────────────────────────────────
def plot_all_samples(adata):
    for sample_id in adata.obs['sample'].unique():
        adata_sample = adata[adata.obs['sample'] == sample_id]
        sc.pl.spatial(
            adata_sample,
            color="broad_type",
            palette=PALETTE,
            size=1.5,
            spot_size=55,
            title=f"Sample: {sample_id}",
            img_key=None,
            show=False,
        )
        plt.savefig(f"../heart_visium/visium_plots/broad_cell_types_{sample_id}.png", dpi=150, bbox_inches="tight")
        plt.show()


# ── Neighbor index ────────────────────────────────────────────────────────────
MAX_NEIGHBOR_DIST = 290  # pixels; set to None to disable the cap

def build_neighbor_index(adata_sub, k=6, max_dist=MAX_NEIGHBOR_DIST):
    samples_arr   = adata_sub.obs['sample'].values
    n_spots       = adata_sub.n_obs
    neighbor_idx  = np.zeros((n_spots, k), dtype=np.int64)
    neighbor_mask = np.zeros((n_spots, k), dtype=bool)
    for sample_id in np.unique(samples_arr):
        g_idx    = np.where(samples_arr == sample_id)[0]
        coords_s = adata_sub.obsm['spatial'][g_idx]
        k_fit    = min(k + 1, len(g_idx))
        if k_fit < 2:
            neighbor_idx[g_idx] = g_idx[:, None]
            continue
        nn = NearestNeighbors(n_neighbors=k_fit, metric='euclidean')
        nn.fit(coords_s)
        dist_local, nbr_local = nn.kneighbors(coords_s)
        dist_local = dist_local[:, 1:]
        nbr_local  = nbr_local[:, 1:]
        if nbr_local.shape[1] < k:
            pad      = np.tile(np.arange(len(g_idx))[:, None], (1, k - nbr_local.shape[1]))
            pad_dist = np.full((len(g_idx), k - nbr_local.shape[1]), np.inf)
            nbr_local  = np.concatenate([nbr_local, pad], axis=1)
            dist_local = np.concatenate([dist_local, pad_dist], axis=1)
        neighbor_idx[g_idx]  = g_idx[nbr_local]
        neighbor_mask[g_idx] = (dist_local <= max_dist) if max_dist is not None else np.isfinite(dist_local)
    kept = neighbor_mask.sum() / max(neighbor_mask.size, 1)
    print(f"Neighbors kept after distance cap ({max_dist} px): {kept:.1%}")
    return neighbor_idx, neighbor_mask, samples_arr


def save_neighbors_to_adata(adata_sub, neighbor_idx, neighbor_mask,
                            save_path='../heart_visium/adata_sub_with_neighbors.h5ad'):
    adata_sub.obsm['neighbor_idx']  = neighbor_idx.astype(np.int64)
    adata_sub.obsm['neighbor_mask'] = neighbor_mask.astype(np.int8)
    adata_sub.uns['neighbor_info']  = {
        'k': int(neighbor_idx.shape[1]),
        'max_dist_px': MAX_NEIGHBOR_DIST if MAX_NEIGHBOR_DIST is not None else 0,
        'description': 'neighbor_idx[i,j] = global index of i-th spot j-th kNN neighbor; '
                       'neighbor_mask[i,j] = 1 if within max_dist',
    }
    adata_sub.write_h5ad(save_path)


def plot_celltype_pair_distribution(neighbor_idx, neighbor_mask, broad_labels,
                                    save_path='../heart_visium/visium_plots/celltype_neighbor_pairs.png'):
    types = sorted(np.unique(broad_labels).tolist())
    t2i   = {t: i for i, t in enumerate(types)}
    nbr_types    = np.array([t2i[t] for t in broad_labels])[neighbor_idx]
    center_types = np.array([t2i[t] for t in broad_labels])
    n = len(types)
    counts = np.zeros((n, n), dtype=np.int64)
    for a in range(n):
        rows = center_types == a
        if not rows.any():
            continue
        sub_nbr  = nbr_types[rows]
        sub_mask = neighbor_mask[rows]
        for b in range(n):
            counts[a, b] = int(((sub_nbr == b) & sub_mask).sum())

    fig, ax = plt.subplots(figsize=(1.4 * n + 2, 1.2 * n + 2))
    im = ax.imshow(counts, cmap='viridis')
    ax.set_xticks(range(n)); ax.set_xticklabels(types, rotation=40, ha='right')
    ax.set_yticks(range(n)); ax.set_yticklabels(types)
    ax.set_xlabel('Neighbor cell type')
    ax.set_ylabel('Center cell type')
    ax.set_title('Observed neighbor pairs by cell type\n(directed; each center -> its valid neighbors)')
    vmax = counts.max() if counts.max() > 0 else 1
    for i in range(n):
        for j in range(n):
            ax.text(j, i, f"{counts[i, j]:,}", ha='center', va='center',
                    color='white' if counts[i, j] < vmax * 0.6 else 'black',
                    fontsize=9)
    plt.colorbar(im, ax=ax, label='# neighbor pairs', shrink=0.7)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    df = pd.DataFrame(counts, index=types, columns=types)
    print("Observed neighbor pairs by cell type:")
    print(df)
    return df


def plot_neighborhood_composition(neighbor_idx, neighbor_mask, broad_labels,
                                  save_path='../heart_visium/visium_plots/neighborhood_composition.png'):
    short = {'macrophage': 'M', 'fibroblast': 'F', 'endothelial': 'E'}
    labels_arr = np.asarray([str(x) for x in broad_labels])

    print("Unique broad labels:", np.unique(labels_arr))
    print("Spots per label:", {l: int((labels_arr == l).sum()) for l in np.unique(labels_arr)})
    print("Mean valid neighbors per spot:", neighbor_mask.sum(axis=1).mean())




    print("Unique broad labels:", np.unique(labels_arr))
    nbr_labels = labels_arr[neighbor_idx]

    categories = ['M', 'F', 'E', 'M+E', 'M+F', 'F+E', 'M+F+E']
    counts = {c: 0 for c in categories}

    for i in range(len(labels_arr)):
        if labels_arr[i] not in short:
            continue
        types = {short[labels_arr[i]]}
        for j in range(neighbor_idx.shape[1]):
            if neighbor_mask[i, j] and nbr_labels[i, j] in short:
                types.add(short[nbr_labels[i, j]])
        if not types or len(types) > 3:
            continue
        key = '+'.join([t for t in ['M', 'F', 'E'] if t in types])
        if key in counts:
            counts[key] += 1

    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(categories, [counts[c] for c in categories],
                  color=['red', 'green', 'blue', 'purple', 'orange', 'teal', 'black'],
                  alpha=0.75, edgecolor='grey')
    for b, c in zip(bars, categories):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height(),
                f'{counts[c]:,}', ha='center', va='bottom', fontsize=9)
    ax.set_ylabel('# neighborhoods')
    ax.set_xlabel('Cell types present (center + neighbors)')
    ax.set_title('Neighborhood composition distribution\n(M=macrophage, F=fibroblast, E=endothelial)')
    ax.grid(alpha=0.2, axis='y')
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()
    print(f"Saved -> {save_path}")
    print("Counts:", counts)
    return counts


def plot_neighbor_count_distribution(neighbor_mask,
                                     save_path='../heart_visium/visium_plots/neighbor_count_distribution.png'):
    counts = neighbor_mask.sum(axis=1)
    k      = neighbor_mask.shape[1]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(counts, bins=np.arange(k + 2) - 0.5,
            color='steelblue', edgecolor='grey', alpha=0.85)
    ax.axvline(np.median(counts), color='red', linestyle='--', linewidth=1.5,
               label=f'median = {np.median(counts):.1f}')
    ax.axvline(np.mean(counts), color='black', linestyle='--', linewidth=1.5,
               label=f'mean = {np.mean(counts):.2f}')
    ax.set_xticks(range(k + 1))
    ax.set_xlabel(f'Valid neighbors per spot (k={k}, max_dist={MAX_NEIGHBOR_DIST} px)')
    ax.set_ylabel('# spots')
    ax.set_title(f'Neighborhood size distribution\n'
                 f'spots with 0 neighbors = {(counts == 0).sum():,} / {len(counts):,}')
    ax.legend()
    ax.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved -> {save_path}")


def plot_neighborhood_example(adata_sub, neighbor_idx, samples_arr, n_centers=5, seed=0):
    sample_id     = np.unique(samples_arr)[0]
    g_idx         = np.where(samples_arr == sample_id)[0]
    coords_all    = adata_sub.obsm['spatial'][g_idx]
    rng           = np.random.default_rng(seed)
    n_centers     = min(n_centers, len(g_idx))
    center_local  = rng.choice(len(g_idx), size=n_centers, replace=False)
    center_global = g_idx[center_local]

    neighbor_global = np.unique(neighbor_idx[center_global].ravel())
    neighbor_global = neighbor_global[~np.isin(neighbor_global, center_global)]

    center_set   = set(center_global)
    neighbor_set = set(neighbor_global)
    colors = ['red' if gi in center_set else 'blue' if gi in neighbor_set else 'lightgray'
              for gi in g_idx]

    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(coords_all[:, 0], coords_all[:, 1], c=colors, s=30, linewidths=0)
    ax.legend(handles=[
        Line2D([0], [0], marker='o', color='w', markerfacecolor='red',       markersize=8, label='Center'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue',      markersize=8, label='Neighbor'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='lightgray', markersize=8, label='Other'),
    ], loc='upper right')
    ax.set_title(f"Neighborhood example — sample {sample_id} ({n_centers} centers, k=6)")
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig("../heart_visium/visium_plots/neighborhood_example.png", dpi=150, bbox_inches="tight")
    plt.show()


# ── Expression matrix ─────────────────────────────────────────────────────────
def prepare_expression(adata):
    X            = adata.X.toarray() if issparse(adata.X) else np.array(adata.X)
    gene_to_idx  = {g: i for i, g in enumerate(adata.var['feature_name'])}
    X_bin        = X > 0
    broad_labels = adata.obs['broad_type'].values
    return X, X_bin, gene_to_idx, broad_labels


# ── Count observed LR neighbor pairs ─────────────────────────────────────────
def count_lr_pairs(lr_pairs, neighbor_idx, neighbor_mask, X_bin, gene_to_idx, broad_labels):
    records = []
    for _, row in lr_pairs.iterrows():
        lig    = row['enriched_gene_ligand']
        rec    = row['enriched_gene_receptor']
        lig_ct = row['cell_type_ligand']
        rec_ct = row['cell_type_receptor']
        if lig not in gene_to_idx or rec not in gene_to_idx:
            continue
        lig_vec   = X_bin[:, gene_to_idx[lig]] & (broad_labels == lig_ct)
        rec_vec   = X_bin[:, gene_to_idx[rec]] & (broad_labels == rec_ct)
        lig_spots = np.where(lig_vec)[0]
        count     = int(sum((rec_vec[neighbor_idx[i]] & neighbor_mask[i]).sum() for i in lig_spots))
        records.append({
            'ligand':             lig,
            'receptor':           rec,
            'archetype_ligand':   row['archetype_ligand'],
            'archetype_receptor': row['archetype_receptor'],
            'cell_type_ligand':   lig_ct,
            'cell_type_receptor': rec_ct,
            'pair_weight':        row['pair_weight'],
            'n_neighbor_pairs':   count,
        })
    lr_counts = (
        pd.DataFrame(records)
        .sort_values('n_neighbor_pairs', ascending=False)
        .reset_index(drop=True)
    )
    print(f"Pairs with spatial data: {len(lr_counts)}")
    return lr_counts


# ── Permutation test 1: neighbor-count ───────────────────────────────────────
def run_permutation_test(lr_counts, X_bin, gene_to_idx, broad_labels, neighbor_idx, neighbor_mask, samples_arr, n_perms=1000):
    n_spots  = len(broad_labels)
    n_pairs  = len(lr_counts)
    observed = lr_counts["n_neighbor_pairs"].values

    pair_lig_idx = np.array([gene_to_idx[g] for g in lr_counts["ligand"]])
    pair_rec_idx = np.array([gene_to_idx[g] for g in lr_counts["receptor"]])
    pair_lig_ct  = lr_counts["cell_type_ligand"].values
    pair_rec_ct  = lr_counts["cell_type_receptor"].values

    ct_mask  = {ct: (broad_labels == ct) for ct in np.unique(np.concatenate([pair_lig_ct, pair_rec_ct]))}
    lig_vecs = [(X_bin[:, pair_lig_idx[j]] & ct_mask[pair_lig_ct[j]]).astype(np.int32)
                for j in range(n_pairs)]

    perm_counts = np.zeros((n_perms, n_pairs), dtype=np.int32)
    rng         = np.random.default_rng(42)

    # Group pairs by receptor cell type; shuffle only that type's spots
    for rec_ct in np.unique(pair_rec_ct):
        group        = np.where(pair_rec_ct == rec_ct)[0]
        unique_ri    = np.unique(pair_rec_idx[group])
        sample_idx_list = []
        for sample_id in np.unique(samples_arr):
            mask = (samples_arr == sample_id) & (broad_labels == rec_ct)
            idx  = np.where(mask)[0]
            if len(idx) > 1:
                sample_idx_list.append(idx)
        for p in trange(n_perms, desc=f"Perm rec_ct={rec_ct}"):
            perm = np.arange(n_spots)
            for idx in sample_idx_list:
                perm[idx] = idx[rng.permutation(len(idx))]
            rec_nbr_sum = {}
            for ri in unique_ri:
                rec_vec_perm = (X_bin[:, ri] & ct_mask[rec_ct])[perm].astype(np.int32)
                rec_nbr_sum[ri] = (rec_vec_perm[neighbor_idx] * neighbor_mask).sum(axis=1)
            for j in group:
                perm_counts[p, j] = lig_vecs[j].dot(rec_nbr_sum[pair_rec_idx[j]])

    p_values = (perm_counts >= observed).sum(axis=0).astype(float)
    p_values = (p_values + 1) / (n_perms + 1)
    _, p_adj, _, _ = multipletests(p_values, method="fdr_bh")

    lr_counts = lr_counts.copy()
    lr_counts["p_value"] = p_values
    lr_counts["p_adj"]   = p_adj
    return lr_counts, perm_counts, observed


def plot_obs_vs_exp(lr_counts, perm_counts, observed, selected_pairs=None):
    expected = perm_counts.mean(axis=0)
    sig_mask = lr_counts['p_adj'] < 0.05

    fig, ax = plt.subplots(figsize=(9, 8))
    ax.scatter(expected[sig_mask], observed[sig_mask],
               c='red', s=30, alpha=0.85, zorder=2,
               label=f'p_adj < 0.05 (n={sig_mask.sum()})')
    ax.scatter(expected[~sig_mask], observed[~sig_mask],
               c='lightgrey', s=20, alpha=0.6, zorder=1, label='Not significant')

    lim = max(expected.max(), observed.max()) * 1.05
    ax.plot([1e-3, lim], [1e-3, lim], 'k--', linewidth=1, alpha=0.4, zorder=0)

    if selected_pairs is not None:
        for entry in selected_pairs:
            if len(entry) == 2:
                lig, rec = entry
                mask = (lr_counts['ligand'] == lig) & (lr_counts['receptor'] == rec)
            elif len(entry) == 4:
                lig, rec, al, ar = entry
                mask = (
                    (lr_counts['ligand']             == lig) &
                    (lr_counts['receptor']           == rec) &
                    (lr_counts['archetype_ligand']   == al)  &
                    (lr_counts['archetype_receptor'] == ar)
                )
            else:
                print(f'Invalid entry (expected 2 or 4 elements): {entry}')
                continue
            matches = np.where(mask)[0]
            if len(matches) == 0:
                print(f'Pair not found: {entry}')
                continue
            for i in matches:
                row   = lr_counts.iloc[i]
                label = (f"{row['ligand']}->{row['receptor']}\n"
                         f"({row['archetype_ligand']}/{row['archetype_receptor']})")
                ax.annotate(label, (expected[i], observed[i]),
                            fontsize=8, alpha=0.75, color='black',
                            xytext=(-55, 10), textcoords='offset points')
                ax.scatter(expected[i], observed[i],
                           s=80, facecolors='none', edgecolors='black',
                           linewidths=1.5, zorder=3)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim(left=max(expected[expected > 0].min() * 0.8, 1e-3))
    ax.set_ylim(bottom=max(observed[observed > 0].min() * 0.8, 1e-3))
    ax.set_xlabel('Expected neighbor pairs (permutation mean) (log)', fontsize=12)
    ax.set_ylabel('Observed neighbor pairs (log)', fontsize=12)
    ax.set_title('Observed vs Expected L-R co-localizing pairs', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.2, which='both')
    plt.tight_layout()
    plt.savefig("../heart_visium/lr_obs_vs_exp.pdf")


# ── Bubble plot ───────────────────────────────────────────────────────────────
def plot_top_significant_LR_pairs(sig_pairs, top_n=30,
                                   p_col='p_adj',
                                   size_col: str | None = 'pair_weight',
                                   save_path='../heart_visium/lr_bubble_plot.pdf'):
    sig_pairs = sig_pairs[sig_pairs[p_col] < 0.05]
    top_pairs = sig_pairs.drop_duplicates(subset=['ligand', 'receptor']).head(top_n)
    sig = sig_pairs.merge(top_pairs[['ligand', 'receptor']], on=['ligand', 'receptor']).copy()
    sig['pair_label'] = sig['ligand'] + '->' + sig['receptor']
    sig['ct_label']   = (sig['cell_type_ligand'].str.capitalize()   + ' (' + sig['archetype_ligand']   + ')'
                         + '->' +
                         sig['cell_type_receptor'].str.capitalize() + ' (' + sig['archetype_receptor'] + ')')

    pair_labels = sig['pair_label'].unique()[::-1]
    ct_labels   = sig['ct_label'].unique()
    pair_to_y   = {p: i for i, p in enumerate(pair_labels)}
    ct_to_x     = {c: i for i, c in enumerate(ct_labels)}

    if size_col is None:
        scaled_sizes = 150
    else:
        sizes    = sig[size_col].fillna(0).abs()
        max_size = sizes.max() or 1
        scaled_sizes = np.sqrt(sizes / max_size) * 400 + 30

    fig, ax = plt.subplots(figsize=(max(5, len(ct_labels) * 2.2), max(6, len(pair_labels) * 0.45)))
    scat = ax.scatter(
        [ct_to_x[c] for c in sig['ct_label']],
        [pair_to_y[p] for p in sig['pair_label']],
        s=scaled_sizes,
        c=sig[p_col], cmap='YlOrRd_r',
        edgecolors='grey', linewidths=0.5, alpha=0.9,
    )
    plt.colorbar(scat, ax=ax, label=p_col, shrink=0.6)
    ax.set_xticks(range(len(ct_labels)))
    ax.set_xticklabels(ct_labels, rotation=40, ha='right', fontsize=8)
    ax.set_yticks(range(len(pair_labels)))
    ax.set_yticklabels(pair_labels, fontsize=8)
    ax.set_xlabel('Archetypes')
    ax.set_ylabel('L-R pair')
    title_size = f'   size = {size_col}' if size_col is not None else ''
    ax.set_title(f'Top {top_n} significant L-R interactions\ncolour = {p_col}{title_size}')
    ax.grid(True, alpha=0.3, linestyle='--')
    if size_col is not None:
        for frac, lbl in [(0.25, round(max_size * 0.25, 2)),
                          (0.60, round(max_size * 0.60, 2)),
                          (1.00, round(max_size,        2))]:
            ax.scatter([], [], s=np.sqrt(frac) * 400 + 30, c='grey', alpha=0.6, label=str(lbl))
        ax.legend(title=size_col, fontsize=8,
                  handletextpad=1.5, borderpad=1.2, labelspacing=1.2)
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    plt.show()


# ── Archetype scoring ─────────────────────────────────────────────────────────
def score_visium_heart(adata, lr_counts, archetype_genes):
    print("  Scoring archetypes...")
    if 'feature_name' in adata.var.columns:
        symbol_to_varname = dict(zip(adata.var['feature_name'], adata.var_names))
        varname_to_idx    = {s: i for i, s in enumerate(adata.var['feature_name'])}
    else:
        symbol_to_varname = {v: v for v in adata.var_names}
        varname_to_idx    = {v: i for i, v in enumerate(adata.var_names)}

    varnames_set = set(adata.var_names)

    def resolve_genes(genes):
        out = []
        for g in genes:
            if g in symbol_to_varname:
                out.append(symbol_to_varname[g])
            elif g in varnames_set:
                out.append(g)
        return out

    arch_list = sorted(archetype_genes.items())
    for i, (arch, genes) in enumerate(arch_list):
        mapped = resolve_genes(genes)
        if len(mapped) < 3:
            print(f"    Skipping {arch} ({len(mapped)} genes found)")
            continue
        sc.tl.score_genes(adata, gene_list=mapped, score_name=f"score_{arch}")
        print(f"    Scored {arch} ({i+1}/{len(arch_list)}, {len(mapped)} genes)")

    arch_pairs = lr_counts[['archetype_ligand', 'archetype_receptor']].drop_duplicates()
    score_cols = [c for c in adata.obs.columns if c.startswith('score_')]
    shifted    = {col.replace('score_', ''): adata.obs[col] - adata.obs[col].min() + 1
                  for col in score_cols}
    for _, row in arch_pairs.iterrows():
        lig, rec = row['archetype_ligand'], row['archetype_receptor']
        if lig not in shifted or rec not in shifted:
            continue
        adata.obs[f"logratio_{lig}_{rec}"] = np.log2(shifted[lig] / shifted[rec])

    X_dense    = adata.X.toarray() if sparse.issparse(adata.X) else np.array(adata.X)
    gene_pairs = lr_counts[['ligand', 'receptor']].drop_duplicates()
    gene_ratios = {}
    for _, row in gene_pairs.iterrows():
        lig_gene, rec_gene = row['ligand'], row['receptor']
        if lig_gene not in varname_to_idx or rec_gene not in varname_to_idx:
            continue
        lig_expr = X_dense[:, varname_to_idx[lig_gene]]
        rec_expr = X_dense[:, varname_to_idx[rec_gene]]
        ratio    = np.log2((lig_expr + 1) / (rec_expr + 1))
        ratio[(lig_expr == 0) & (rec_expr == 0)] = np.nan
        gene_ratios[f"generatio_{lig_gene}_{rec_gene}"] = ratio

    adata.obs = pd.concat([adata.obs, pd.DataFrame(gene_ratios, index=adata.obs.index)], axis=1)
    print(f"  Done: {len(score_cols)} archetypes scored, {len(gene_ratios)} gene ratios computed")
    return adata


# ── Archetype vs gene log-ratio scatter ───────────────────────────────────────
def compute_pearson_r(adata, lr_counts):
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.duplicated()]
    if 'feature_name' in adata.var.columns:
        symbol_to_var = dict(zip(adata.var['feature_name'], adata.var_names))
    else:
        symbol_to_var = {v: v for v in adata.var_names}

    def gene_expr(g):
        var = symbol_to_var.get(g, g)
        if var not in adata.var_names:
            return None
        col = adata[:, var].X
        return col.toarray().ravel() if sparse.issparse(col) else np.asarray(col).ravel()

    rs = []
    for _, row in lr_counts.iterrows():
        arch_col = f"logratio_{row['archetype_ligand']}_{row['archetype_receptor']}"
        gene_col = f"generatio_{row['ligand']}_{row['receptor']}"
        if arch_col not in adata.obs or gene_col not in adata.obs:
            rs.append(np.nan); continue
        x = adata.obs[arch_col].values.astype(float)
        y = adata.obs[gene_col].values.astype(float)
        m = np.isfinite(x) & np.isfinite(y)
        lig_vec = gene_expr(row['ligand'])
        rec_vec = gene_expr(row['receptor'])
        if lig_vec is not None:
            m &= (lig_vec > 0)
        if rec_vec is not None:
            m &= (rec_vec > 0)
        rs.append(pearsonr(x[m], y[m])[0] if m.sum() > 2 else np.nan)
    lr_counts = lr_counts.copy()
    lr_counts['pearson_r'] = rs
    return lr_counts


def _neighbor_mean(arr, neighbor_idx, neighbor_mask):
    nbr_vals = arr[neighbor_idx]
    valid = neighbor_mask.astype(float)
    n_valid = valid.sum(axis=1)
    s = (nbr_vals * valid).sum(axis=1)
    out = np.divide(s, n_valid, out=np.zeros_like(s, dtype=float), where=n_valid > 0)
    return out





def plot_ligand_receptor_archetype_scatter(adata, lr_counts,
        save_path='../heart_visium/lr_ligand_receptor_archetype_scatter.pdf'):
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.duplicated()]
    if 'feature_name' in adata.var.columns:
        symbol_to_var = dict(zip(adata.var['feature_name'], adata.var_names))
    else:
        symbol_to_var = {v: v for v in adata.var_names}

    def gene_expr(g):
        var = symbol_to_var.get(g, g)
        if var not in adata.var_names:
            return None
        col = adata[:, var].X
        return col.toarray().ravel() if sparse.issparse(col) else np.asarray(col).ravel()

    ct_all = adata.obs['broad_type'].values
    with PdfPages(save_path) as pdf:
        for _, row in lr_counts.reset_index(drop=True).iterrows():
            lig, rec = row['ligand'], row['receptor']
            ax_arch  = row['archetype_ligand']
            ay_arch  = row['archetype_receptor']
            ct_lig   = broad_category(row['cell_type_ligand'])
            ct_rec   = broad_category(row['cell_type_receptor'])
            s_lig    = adata.obs.get(f"score_{ax_arch}")
            s_rec    = adata.obs.get(f"score_{ay_arch}")
            L        = gene_expr(lig)
            R        = gene_expr(rec)
            if s_lig is None or s_rec is None or L is None or R is None:
                continue

            fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

            mask_x = (ct_all == ct_lig) & (L > 0)
            xL = L[mask_x]
            yL = np.asarray(s_lig, dtype=float)[mask_x]
            if len(xL) >= 3 and xL.std() > 0 and yL.std() > 0:
                rL, pL = pearsonr(xL, yL)
                axes[0].scatter(xL, yL, s=6, alpha=0.5,
                                color=CT_COLORS.get(ct_lig, 'black'), rasterized=True)
                axes[0].set_title(f"{lig} in {ct_lig} cells (n={int(mask_x.sum())})\n"
                                  f"r = {rL:.3f}   p = {pL:.2e}", fontsize=10)
            else:
                axes[0].set_title(f"{lig} in {ct_lig} (n={int(mask_x.sum())}) — too few", fontsize=10)
            axes[0].set_xlabel(f"{lig} expression")
            axes[0].set_ylabel(f"score_{ax_arch}")
            axes[0].grid(alpha=0.2)

            mask_y = (ct_all == ct_rec) & (R > 0)
            xR = R[mask_y]
            yR = np.asarray(s_rec, dtype=float)[mask_y]
            if len(xR) >= 3 and xR.std() > 0 and yR.std() > 0:
                rR, pR = pearsonr(xR, yR)
                axes[1].scatter(xR, yR, s=6, alpha=0.5,
                                color=CT_COLORS.get(ct_rec, 'black'), rasterized=True)
                axes[1].set_title(f"{rec} in {ct_rec} cells (n={int(mask_y.sum())})\n"
                                  f"r = {rR:.3f}   p = {pR:.2e}", fontsize=10)
            else:
                axes[1].set_title(f"{rec} in {ct_rec} (n={int(mask_y.sum())}) — too few", fontsize=10)
            axes[1].set_xlabel(f"{rec} expression")
            axes[1].set_ylabel(f"score_{ay_arch}")
            axes[1].grid(alpha=0.2)

            fig.suptitle(f"{lig} ({ax_arch}) -> {rec} ({ay_arch})", fontsize=12)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
    print(f"Saved {save_path}")


def select_doubly_significant(lr_counts, alpha=0.05,
                              save_path='../heart_visium/lr_doubly_significant.csv'):
    mask = (lr_counts['p_adj'] < alpha) & (lr_counts['weighted_p_adj'] < alpha)
    out  = lr_counts.loc[mask].copy()
    out  = out.sort_values('pearson_r', ascending=False).reset_index(drop=True)
    cols = ['ligand', 'receptor', 'cell_type_ligand', 'cell_type_receptor',
            'archetype_ligand', 'archetype_receptor',
            'pearson_r', 'p_adj', 'weighted_p_adj',
            'n_neighbor_pairs', 'weighted_pair_score']
    cols = [c for c in cols if c in out.columns]
    print(f"Pairs significant in both p_adj and weighted_p_adj (alpha={alpha}): {len(out)}")
    with pd.option_context('display.max_rows', None, 'display.width', 200):
        print(out[cols].to_string(index=False))
    out.to_csv(save_path, index=False)
    print(f"Saved -> {save_path}")
    return out


def plot_arch_vs_gene_logratio(adata, lr_counts,
                               ligand=None, receptor=None,
                               archetype1=None, archetype2=None):
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.duplicated()]

    if ligand is not None:
        lr_counts = lr_counts[lr_counts['ligand'] == ligand]
    if receptor is not None:
        lr_counts = lr_counts[lr_counts['receptor'] == receptor]
    if archetype1 is not None:
        lr_counts = lr_counts[lr_counts['archetype_ligand'] == archetype1]
    if archetype2 is not None:
        lr_counts = lr_counts[lr_counts['archetype_receptor'] == archetype2]

    if 'feature_name' in adata.var.columns:
        symbol_to_var = dict(zip(adata.var['feature_name'], adata.var_names))
    else:
        symbol_to_var = {v: v for v in adata.var_names}

    def gene_expr(g):
        var = symbol_to_var.get(g, g)
        if var not in adata.var_names:
            return None
        col = adata[:, var].X
        return (col.toarray().ravel() if sparse.issparse(col) else np.asarray(col).ravel())

    pairs = lr_counts.reset_index(drop=True)
    correlations = []
    for _, row in pairs.iterrows():
        arch_col = f"logratio_{row['archetype_ligand']}_{row['archetype_receptor']}"
        gene_col = f"generatio_{row['ligand']}_{row['receptor']}"
        if arch_col not in adata.obs or gene_col not in adata.obs:
            correlations.append(np.nan)
            continue
        x    = adata.obs[arch_col].values.astype(float)
        y    = adata.obs[gene_col].values.astype(float)
        lig_vec = gene_expr(row['ligand'])
        rec_vec = gene_expr(row['receptor'])
        mask = np.isfinite(x) & np.isfinite(y)
        if lig_vec is not None:
            mask &= (lig_vec > 0)
        if rec_vec is not None:
            mask &= (rec_vec > 0)
        r, _ = pearsonr(x[mask], y[mask]) if mask.sum() > 2 else (np.nan, np.nan)
        correlations.append(r)

    pairs['pearson_r'] = correlations
    pairs = pairs.sort_values('pearson_r', ascending=False).reset_index(drop=True)

    suffix = ''
    for tag in (ligand, receptor, archetype1, archetype2):
        if tag is not None:
            suffix += f'_{tag}'
    with PdfPages(f'../heart_visium/lr_arch_vs_gene_logratio{suffix}.pdf') as pdf:
        for _, row in pairs.iterrows():
            arch_col = f"logratio_{row['archetype_ligand']}_{row['archetype_receptor']}"
            gene_col = f"generatio_{row['ligand']}_{row['receptor']}"
            if arch_col not in adata.obs or gene_col not in adata.obs:
                continue
            x    = adata.obs[arch_col].values.astype(float)
            y    = adata.obs[gene_col].values.astype(float)
            ct   = adata.obs['broad_type'].values
            lig_vec = gene_expr(row['ligand'])
            rec_vec = gene_expr(row['receptor'])
            mask = np.isfinite(x) & np.isfinite(y)
            if lig_vec is not None:
                mask &= (lig_vec > 0)
            if rec_vec is not None:
                mask &= (rec_vec > 0)
            x, y, ct = x[mask], y[mask], ct[mask]
            r, pval  = pearsonr(x, y) if len(x) > 2 else (np.nan, np.nan)
            pstr     = f"p = {pval:.2e}" if not np.isnan(pval) else ""
            relevant_cts = {row['cell_type_ligand'], row['cell_type_receptor']}
            fig, ax  = plt.subplots(figsize=(6, 5))
            other_mask = ~np.isin(ct, list(relevant_cts))
            if other_mask.any():
                ax.scatter(x[other_mask], y[other_mask], s=10, alpha=0.25,
                           color='lightgrey', edgecolors='none',
                           rasterized=False, zorder=1)
            leg_handles = []
            for ct_name in sorted(relevant_cts):
                color   = CT_COLORS.get(ct_name, 'black')
                ct_mask = ct == ct_name
                if ct_mask.any():
                    ax.scatter(x[ct_mask], y[ct_mask], s=14, alpha=0.7,
                               color=color, edgecolors='none',
                               rasterized=False, zorder=2)
                leg_handles.append(
                    plt.Line2D([0], [0], marker='o', color='w',
                               markerfacecolor=color, markersize=9, label=ct_name)
                )
            if other_mask.any():
                leg_handles.append(
                    plt.Line2D([0], [0], marker='o', color='w',
                               markerfacecolor='lightgrey', markersize=8, label='other')
                )
            ax.axhline(0, color='grey', linewidth=0.7, linestyle='--')
            ax.axvline(0, color='grey', linewidth=0.7, linestyle='--')
            ax.set_title(
                f"{row['ligand']} -> {row['receptor']}\n"
                f"Archetypes: {row['archetype_ligand']} / {row['archetype_receptor']}\n"
                f"Pearson r = {r:.3f}    {pstr}", fontsize=11)
            ax.set_xlabel(f"log2({row['archetype_ligand']} / {row['archetype_receptor']}) archetype score", fontsize=9)
            ax.set_ylabel(f"log2({row['ligand']} / {row['receptor']}) expression", fontsize=9)
            ax.legend(handles=leg_handles, fontsize=9, title='Cell type', title_fontsize=9,
                      framealpha=0.8, loc='best')
            ax.grid(alpha=0.2)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight', dpi=600)
            plt.close(fig)
    print(f"Saved lr_arch_vs_gene_logratio{suffix}.pdf")


# ── Spatial map for a specific LR pair ───────────────────────────────────────
def plot_spatial_lr_pair(adata, lr_counts, pair_ligand, pair_receptor, selected_arch_pairs=None):
    all_rows = lr_counts[
        (lr_counts['ligand']   == pair_ligand) &
        (lr_counts['receptor'] == pair_receptor)
    ]
    if selected_arch_pairs:
        arch_rows = [
            all_rows[
                (all_rows['archetype_ligand']   == al) &
                (all_rows['archetype_receptor'] == ar)
            ].iloc[0]
            for al, ar in selected_arch_pairs
        ]
    else:
        arch_rows = [row for _, row in all_rows.iterrows()]

    if 'feature_name' in adata.var.columns:
        gene_to_var = dict(zip(adata.var['feature_name'], adata.var_names))
    else:
        gene_to_var = {v: v for v in adata.var_names}
    var_list = list(adata.var_names)

    def get_expr(adata_obj, gene):
        varname = gene_to_var.get(gene, gene)
        if varname not in adata_obj.var_names:
            return np.zeros(adata_obj.n_obs)
        idx = var_list.index(varname)
        x   = adata_obj.X[:, idx]
        return x.toarray().flatten() if sparse.issparse(adata_obj.X) else np.asarray(x).flatten()

    ncols = 1 + len(arch_rows)
    with PdfPages(f'../heart_visium/{pair_ligand}_{pair_receptor}_spatial.pdf') as pdf:
        for sample_id in sorted(adata.obs['sample'].unique()):
            adata_s  = adata[adata.obs['sample'] == sample_id]
            coords   = adata_s.obsm['spatial']
            broad    = adata_s.obs['broad_type'].values
            lig_expr = get_expr(adata_s, pair_ligand)
            rec_expr = get_expr(adata_s, pair_receptor)
            all_cts  = {r['cell_type_ligand'] for r in arch_rows} | {r['cell_type_receptor'] for r in arch_rows}
            other    = ~np.isin(broad, list(all_cts))

            fig, axes = plt.subplots(1, ncols, figsize=(7 * ncols, 6))
            if ncols == 1:
                axes = [axes]

            ax = axes[0]
            ax.scatter(coords[other, 0], coords[other, 1],
                       c='lightgrey', s=2, alpha=0.2, rasterized=True, zorder=1)
            ct_gene_shown = set()
            for row in arch_rows:
                for ct, gene, cmap in [
                    (row['cell_type_ligand'],   pair_ligand,   'Reds'),
                    (row['cell_type_receptor'], pair_receptor, 'Blues'),
                ]:
                    if (ct, gene) in ct_gene_shown:
                        continue
                    ct_gene_shown.add((ct, gene))
                    expr = lig_expr if gene == pair_ligand else rec_expr
                    mask = broad == ct
                    if not mask.any():
                        continue
                    vmax = np.percentile(expr[mask], 95)
                    sc_p = ax.scatter(coords[mask, 0], coords[mask, 1],
                                      c=expr[mask], cmap=cmap, s=6,
                                      vmin=0, vmax=max(vmax, 1),
                                      rasterized=True, zorder=2)
                    cb = plt.colorbar(sc_p, ax=ax, shrink=0.4, pad=0.02)
                    cb.set_label(f'{gene} ({ct})', fontsize=8)
            ax.set_title(f'Gene expression\n{pair_ligand} -> {pair_receptor}', fontsize=10)
            ax.axis('off')

            for col_idx, row in enumerate(arch_rows, start=1):
                arch_lig = row['archetype_ligand']
                arch_rec = row['archetype_receptor']
                ct_lig   = row['cell_type_ligand']
                ct_rec   = row['cell_type_receptor']
                ax2 = axes[col_idx]
                ax2.scatter(coords[other, 0], coords[other, 1],
                            c='lightgrey', s=2, alpha=0.2, rasterized=True, zorder=1)
                for ct, arch, cmap, pad in [
                    (ct_lig, arch_lig, 'Reds',  0.02),
                    (ct_rec, arch_rec, 'Blues', 0.08),
                ]:
                    mask = broad == ct
                    col  = f'score_{arch}'
                    if col not in adata_s.obs.columns or not mask.any():
                        continue
                    scores = adata_s.obs[col].values[mask].astype(float)
                    valid  = np.isfinite(scores)
                    if not valid.any():
                        continue
                    sc_p = ax2.scatter(coords[mask][valid, 0], coords[mask][valid, 1],
                                       c=scores[valid], cmap=cmap, s=6,
                                       rasterized=True, zorder=2)
                    cb = plt.colorbar(sc_p, ax=ax2, shrink=0.4, pad=pad)
                    cb.set_label(f'{arch} score ({ct})', fontsize=8)
                leg = [plt.Line2D([0], [0], marker='o', color='w',
                                   markerfacecolor=c, markersize=7, label=lbl)
                       for c, lbl in [('lightgrey', 'other'),
                                      ('red',  f'{ct_lig} ({arch_lig})'),
                                      ('blue', f'{ct_rec} ({arch_rec})')]]
                ax2.legend(handles=leg, fontsize=7, title='Cell type / archetype',
                           title_fontsize=7, loc='upper left',
                           bbox_to_anchor=(1.02, 1.0), borderaxespad=0.,
                           framealpha=0.8)
                ax2.set_title(f'Archetypes\n{ct_lig}: {arch_lig}  |  {ct_rec}: {arch_rec}', fontsize=10)
                ax2.axis('off')

            fig.suptitle(f'Sample: {sample_id}  {pair_ligand} -> {pair_receptor}', fontsize=12)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
    print(f"Saved {pair_ligand}_{pair_receptor}_spatial.pdf")

# ── Neighbor adjacency from kNN index ────────────────────────────────────────
def build_neighbor_adjacency(neighbor_idx, neighbor_mask):
    n_sub, k = neighbor_idx.shape
    rows = np.repeat(np.arange(n_sub), k)[neighbor_mask.ravel()]
    cols = neighbor_idx.ravel()[neighbor_mask.ravel()]
    A = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.float32), (rows, cols)),
        shape=(n_sub, n_sub),
    )
    print(f"kNN neighbor pairs after cap: {A.nnz:,}  ({A.nnz / n_sub:.1f} per spot)")
    return A


def plot_archetype_score_distributions(adata, save_path='../heart_visium/archetype_score_distributions.pdf'):
    score_cols = sorted(c for c in adata.obs.columns if c.startswith('score_'))
    broad      = adata.obs['broad_type'].values
    with PdfPages(save_path) as pdf:
        for col in score_cols:
            arch   = col.replace('score_', '')
            scores = adata.obs[col].values.astype(float)
            valid  = np.isfinite(scores)
            fig, ax = plt.subplots(figsize=(6, 4))
            ax.hist(scores[valid], bins=60, color='steelblue',
                    edgecolor='grey', alpha=0.85)
            ax.axvline(np.median(scores[valid]), color='red', linestyle='--',
                       linewidth=1.5, label=f'median = {np.median(scores[valid]):.3f}')
            ax.axvline(np.mean(scores[valid]), color='black', linestyle='--',
                       linewidth=1.5, label=f'mean = {np.mean(scores[valid]):.3f}')
            ax.set_xlabel(f'{arch} score', fontsize=10)
            ax.set_ylabel('# spots', fontsize=10)
            ax.set_title(f'{arch}  (n={valid.sum():,} spots)', fontsize=11)
            ax.legend(fontsize=9)
            ax.grid(alpha=0.2)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
    print(f"Saved {len(score_cols)} pages -> {save_path}")


def build_score_lookup(adata):
    return {
        col.replace('score_', ''): np.clip(
            adata.obs[col].values.astype(np.float32), 0, None)
        for col in adata.obs.columns if col.startswith('score_')
    }


# ── Permutation test 2: distance-weighted ────────────────────────────────────
def run_weighted_permutation_test(lr_counts, X, gene_to_idx, score_lookup, A_nbr, samples_arr, broad_labels, n_perms=1000):
    n_sub   = X.shape[0]
    n_pairs = len(lr_counts)

    pair_keys_w = list(zip(lr_counts['ligand'],
                           lr_counts['receptor'],
                           lr_counts['archetype_ligand'],
                           lr_counts['archetype_receptor']))
    pair_lig_ct = lr_counts['cell_type_ligand'].values
    pair_rec_ct = lr_counts['cell_type_receptor'].values
    ct_mask     = {ct: (broad_labels == ct).astype(np.float32)
                   for ct in np.unique(np.concatenate([pair_lig_ct, pair_rec_ct]))}

    lig_w_mat = np.zeros((n_sub, n_pairs), dtype=np.float32)
    rec_w_mat = np.zeros((n_sub, n_pairs), dtype=np.float32)
    valid_w   = np.ones(n_pairs, dtype=bool)

    for k, (L, R, AL, AR) in enumerate(pair_keys_w):
        if L not in gene_to_idx or R not in gene_to_idx or AL not in score_lookup or AR not in score_lookup:
            valid_w[k] = False
            continue
        lig_w_mat[:, k] = (X[:, gene_to_idx[L]].astype(np.float32)
                           * score_lookup[AL] * ct_mask[pair_lig_ct[k]])
        rec_w_mat[:, k] = (X[:, gene_to_idx[R]].astype(np.float32)
                           * score_lookup[AR] * ct_mask[pair_rec_ct[k]])

    weighted_scores = []
    n_pairs_within  = []
    for k, (L, R, AL, AR) in enumerate(pair_keys_w):
        if not valid_w[k]:
            weighted_scores.append(np.nan)
            n_pairs_within.append(0)
            continue
        score   = float(lig_w_mat[:, k] @ A_nbr @ rec_w_mat[:, k])
        lig_bin = (X[:, gene_to_idx[L]] > 0).astype(np.float32)
        rec_bin = (X[:, gene_to_idx[R]] > 0).astype(np.float32)
        weighted_scores.append(score)
        n_pairs_within.append(int(lig_bin @ A_nbr @ rec_bin))

    lr_counts = lr_counts.copy()
    lr_counts['weighted_pair_score']     = weighted_scores
    lr_counts['n_neighbor_pairs_binary'] = n_pairs_within

    observed_w      = lr_counts['weighted_pair_score'].values.astype(np.float32)
    perm_scores_w   = np.zeros((n_perms, n_pairs), dtype=np.float32)
    pair_key_to_col = {k: i for i, k in enumerate(pair_keys_w)}
    rng_w           = np.random.default_rng(42)

    # Group pairs by receptor cell type; shuffle only that type's spots
    for rec_ct in np.unique(pair_rec_ct):
        group = np.where((pair_rec_ct == rec_ct) & valid_w)[0]
        if len(group) == 0:
            continue
        sample_idx_list = []
        for sample_id in np.unique(samples_arr):
            mask = (samples_arr == sample_id) & (broad_labels == rec_ct)
            idx  = np.where(mask)[0]
            if len(idx) > 1:
                sample_idx_list.append(idx)
        sub_lig = lig_w_mat[:, group]
        sub_rec = rec_w_mat[:, group]
        for p in trange(n_perms, desc=f"Weighted perm rec_ct={rec_ct}"):
            perm = np.arange(n_sub)
            for idx in sample_idx_list:
                perm[idx] = idx[rng_w.permutation(len(idx))]
            rec_perm = sub_rec[perm]
            Arec     = A_nbr @ rec_perm
            perm_scores_w[p, group] = (sub_lig * Arec).sum(axis=0)

    ge_counts  = (perm_scores_w >= observed_w[None, :]).sum(axis=0).astype(np.int32)
    p_values_w = (ge_counts + 1) / (n_perms + 1)
    p_values_w[~valid_w] = np.nan

    p_adj_w = np.full(n_pairs, np.nan)
    mask_w  = ~np.isnan(p_values_w)
    _, p_adj_w[mask_w], _, _ = multipletests(p_values_w[mask_w], method='fdr_bh')

    lr_counts['weighted_ge_count'] = ge_counts
    lr_counts['weighted_p_value']  = p_values_w
    lr_counts['weighted_p_adj']    = p_adj_w

    return lr_counts, perm_scores_w, observed_w, pair_key_to_col


def plot_perm_dists(lr_counts, perm_counts, observed, n_perms=1000,
                    save_path='../heart_visium/lr_permutation_dists.pdf'):
    rows_sorted = lr_counts.sort_values('p_value').reset_index(drop=True)
    # map (ligand, receptor, AL, AR) to column in perm_counts using original lr_counts order
    keys = list(zip(lr_counts['ligand'], lr_counts['receptor'],
                    lr_counts['archetype_ligand'], lr_counts['archetype_receptor']))
    key_to_col = {k: i for i, k in enumerate(keys)}
    with PdfPages(save_path) as pdf:
        for _, row in rows_sorted.iterrows():
            key = (row['ligand'], row['receptor'],
                   row['archetype_ligand'], row['archetype_receptor'])
            if key not in key_to_col:
                continue
            col  = key_to_col[key]
            null = perm_counts[:, col]
            obs  = observed[col]
            max_val = max(int(null.max()), int(obs))
            bins = np.arange(max_val + 2) - 0.5 if max_val < 80 else 50
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.hist(null, bins=bins, color='lightsteelblue', edgecolor='grey',
                    alpha=0.85, label=f'permutations (n={n_perms})')
            ax.axvline(obs, color='red', linewidth=2, label=f'observed = {obs}')
            ax.axvline(np.median(null), color='black', linestyle='--',
                       linewidth=1, label=f'null median = {np.median(null):.1f}')
            ax.set_xlabel('Neighbor-pair count', fontsize=10)
            ax.set_ylabel('# permutations', fontsize=10)
            ax.set_title(
                f"{row['ligand']} ({row['cell_type_ligand']}, {row['archetype_ligand']}) -> "
                f"{row['receptor']} ({row['cell_type_receptor']}, {row['archetype_receptor']})\n"
                f"p = {row['p_value']:.4f}    p_adj = {row['p_adj']:.4f}",
                fontsize=10)
            ax.legend(fontsize=8, loc='best')
            ax.grid(alpha=0.2)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
    print(f"Saved -> {save_path}")


def plot_weighted_perm_dists(lr_counts, perm_scores_w, pair_key_to_col, n_perms_w=1000):
    rows_sorted = lr_counts.sort_values('weighted_p_value').reset_index(drop=True)
    with PdfPages('../heart_visium/lr_weighted_permutation_dists.pdf') as pdf:
        for _, row in rows_sorted.iterrows():
            key = (row['ligand'], row['receptor'],
                   row['archetype_ligand'], row['archetype_receptor'])
            if key not in pair_key_to_col:
                continue
            col  = pair_key_to_col[key]
            null = perm_scores_w[:, col]
            obs  = row['weighted_pair_score']
            fig, ax = plt.subplots(figsize=(7, 4.5))
            ax.hist(null, bins=50, color='lightsteelblue', edgecolor='grey',
                    alpha=0.85, label=f'permutations (n={n_perms_w})')
            ax.axvline(obs, color='red', linewidth=2, label=f'observed = {obs:.2f}')
            ax.axvline(np.median(null), color='black', linestyle='--',
                       linewidth=1, label=f'null median = {np.median(null):.2f}')
            ax.set_xlabel('Distance-weighted L×R score', fontsize=10)
            ax.set_ylabel('# permutations', fontsize=10)
            ax.set_title(
                f"{row['ligand']} ({row['cell_type_ligand']}, {row['archetype_ligand']}) -> "
                f"{row['receptor']} ({row['cell_type_receptor']}, {row['archetype_receptor']})\n"
                f"ge_count = {int(row['weighted_ge_count'])}/{n_perms_w}    "
                f"p = {row['weighted_p_value']:.4f}    p_adj = {row['weighted_p_adj']:.4f}",
                fontsize=10)
            ax.legend(fontsize=8, loc='best')
            ax.grid(alpha=0.2)
            plt.tight_layout()
            pdf.savefig(fig, bbox_inches='tight')
            plt.close(fig)
    print("Saved lr_weighted_permutation_dists.pdf")


def plot_volcano_weighted(lr_counts, perm_scores_w, observed_w, alpha=0.05,
                          selected_pairs=None,
                          save_path='../heart_visium/lr_volcano_weighted.pdf'):
    df       = lr_counts.reset_index(drop=True)
    expected = perm_scores_w.mean(axis=0)
    eps      = np.nextafter(0, 1, dtype=np.float64)
    obs_safe = np.where(np.isfinite(observed_w), observed_w, 0.0)
    finite_ratio = (expected > 0) & (obs_safe > 0)
    log2fc   = np.full_like(obs_safe, np.nan, dtype=float)
    log2fc[finite_ratio] = np.log2(obs_safe[finite_ratio] / expected[finite_ratio])
    if np.isfinite(log2fc).any():
        left_sentinel  = float(np.nanmin(log2fc) - 1)
        right_sentinel = float(np.nanmax(log2fc) + 1)
    else:
        left_sentinel, right_sentinel = -1.0, 1.0
    log2fc[(obs_safe == 0) & (expected > 0)] = left_sentinel
    log2fc[(obs_safe > 0)  & (expected == 0)] = right_sentinel
    padj     = df['weighted_p_adj'].values.astype(float)
    neglogp  = -np.log10(np.clip(padj, eps, 1))

    valid    = np.isfinite(log2fc) & np.isfinite(neglogp)
    sig      = valid & (padj < alpha)

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(log2fc[valid & ~sig], neglogp[valid & ~sig],
               c='lightgrey', s=18, alpha=0.6, label='Not significant')
    ax.scatter(log2fc[sig], neglogp[sig],
               c='red', s=28, alpha=0.85, label=f'p_adj < {alpha} (n={sig.sum()})')
    ax.axhline(-np.log10(alpha), color='black', linestyle='--', linewidth=1, alpha=0.5)
    ax.axvline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)

    if selected_pairs is not None:
        for entry in selected_pairs:
            if len(entry) == 2:
                lig, rec = entry
                mask = (df['ligand'] == lig) & (df['receptor'] == rec)
            elif len(entry) == 4:
                lig, rec, al, ar = entry
                mask = ((df['ligand']==lig) & (df['receptor']==rec)
                        & (df['archetype_ligand']==al) & (df['archetype_receptor']==ar))
            else:
                continue
            for i in np.where(mask & valid)[0]:
                row = df.iloc[i]
                ax.scatter(log2fc[i], neglogp[i], s=90, facecolors='none',
                           edgecolors='black', linewidths=1.5, zorder=3)
                ax.annotate(f"{row['ligand']}->{row['receptor']}\n"
                            f"({row['archetype_ligand']}/{row['archetype_receptor']})",
                            (log2fc[i], neglogp[i]),
                            fontsize=8, xytext=(6, 6), textcoords='offset points')

    ax.set_xlabel('log2(observed / expected)  [weighted]', fontsize=12)
    ax.set_ylabel('-log10(weighted_p_adj)', fontsize=12)
    ax.set_title('Volcano - distance-weighted L×R score', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.2)
    plt.tight_layout()
    plt.savefig(save_path)
    print(f"Saved -> {save_path}")


def plot_obs_vs_exp_weighted(lr_counts, perm_scores_w, observed_w, selected_pairs=None):
    expected_w = perm_scores_w.mean(axis=0)
    df         = lr_counts.reset_index(drop=True)
    sig_mask   = df['weighted_p_adj'] < 0.05

    fig, ax = plt.subplots(figsize=(9, 8))
    ax.scatter(expected_w[sig_mask], observed_w[sig_mask],
               c='red', s=30, alpha=0.85, zorder=2,
               label=f'p_adj < 0.05 (n={sig_mask.sum()})')
    ax.scatter(expected_w[~sig_mask], observed_w[~sig_mask],
               c='lightgrey', s=20, alpha=0.6, zorder=1, label='Not significant')

    pos = (expected_w > 0) & (observed_w > 0)
    if pos.any():
        lim = max(expected_w[pos].max(), observed_w[pos].max()) * 1.05
        ax.plot([1e-3, lim], [1e-3, lim], 'k--', linewidth=1, alpha=0.4, zorder=0)

    if selected_pairs is not None:
        for entry in selected_pairs:
            if len(entry) == 2:
                lig, rec = entry
                mask = (df['ligand'] == lig) & (df['receptor'] == rec)
            elif len(entry) == 4:
                lig, rec, al, ar = entry
                mask = (
                    (df['ligand']             == lig) &
                    (df['receptor']           == rec) &
                    (df['archetype_ligand']   == al)  &
                    (df['archetype_receptor'] == ar)
                )
            else:
                print(f'Invalid entry (expected 2 or 4 elements): {entry}')
                continue
            matches = np.where(mask)[0]
            if len(matches) == 0:
                print(f'Pair not found: {entry}')
                continue
            for i in matches:
                row   = df.iloc[i]
                label = (f"{row['ligand']}->{row['receptor']}\n"
                         f"({row['archetype_ligand']}/{row['archetype_receptor']})")
                ax.annotate(label, (expected_w[i], observed_w[i]),
                            fontsize=8, alpha=0.75, color='black',
                            xytext=(-55, 10), textcoords='offset points')
                ax.scatter(expected_w[i], observed_w[i],
                           s=80, facecolors='none', edgecolors='black',
                           linewidths=1.5, zorder=3)

    ax.set_xscale('log')
    ax.set_yscale('log')
    if pos.any():
        ax.set_xlim(left=expected_w[pos].min() * 0.8)
        ax.set_ylim(bottom=observed_w[pos].min() * 0.8)
    ax.set_xlabel('Expected weighted score (permutation mean, log)', fontsize=12)
    ax.set_ylabel('Observed weighted score (log)', fontsize=12)
    ax.set_title('Observed vs Expected  -  distance-weighted L×R score', fontsize=13)
    ax.legend(fontsize=9)
    ax.grid(alpha=0.2, which='both')
    plt.tight_layout()
    plt.savefig("../heart_visium/lr_obs_vs_exp_weighted.pdf")


def plot_bubble_plots(lr_counts):
    lr_counts = lr_counts[(lr_counts['p_adj'] < 0.05)]
    # lr_cross = lr_counts[lr_counts['cell_type_ligand'] != lr_counts['cell_type_receptor']]
    plot_top_significant_LR_pairs(
        lr_counts.sort_values(by=['p_adj', 'pearson_r'], ascending=[True, False]),
        p_col='p_adj',
        size_col='pearson_r',
        save_path='../heart_visium/lr_bubble_plot.pdf',
    )
    # plot_top_significant_LR_pairs(
    #     lr_counts.dropna(subset=['weighted_p_adj'])
    #             .sort_values(by=['weighted_p_adj', 'weighted_pair_score'], ascending=[True, False]),
    #     p_col='weighted_p_adj',
    #     size_col=None,
    #     save_path='../heart_visium/lr_bubble_plot_weighted.pdf',
    # )

def van(lr_counts):
    from matplotlib_venn import venn2

    sig1_idx = set(lr_counts.index[lr_counts['p_adj'] < 0.05])
    sig2_idx = set(lr_counts.index[lr_counts['weighted_p_adj'] < 0.05])

    fig, ax = plt.subplots(figsize=(6, 5))
    v = venn2(
        [sig1_idx, sig2_idx],
        set_labels=(
            f"6-NN binary\n(p_adj < 0.05, n={len(sig1_idx)})",
            f"Distance-weighted\n(w_p_adj < 0.05, n={len(sig2_idx)})",
        ),
        ax=ax,
    )
    ax.set_title("Overlap between two permutation tests", fontsize=13)
    plt.tight_layout()
    plt.savefig("../heart_visium/lr_method_venn.pdf", bbox_inches="tight")
    plt.show()


if __name__ == '__main__':
    lr_counts = pd.read_csv('../heart_visium/lr_counts.csv')
    plot_bubble_plots(lr_counts)

    # ── Load data (required by everything below) ───────────────────────────────
    adata = load_data()
    lr_pairs = load_lr_pairs()
    archetype_genes = load_archetype_genes()

    # plot_all_samples(adata)

    # ── Normalize once; everything downstream uses adata.X ─────────────────────
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # ── Build spatial structures & expression matrix ───────────────────────────
    neighbor_idx, neighbor_mask, samples_arr = build_neighbor_index(adata)
    plot_neighbor_count_distribution(neighbor_mask)
    # save_neighbors_to_adata(adata, neighbor_idx, neighbor_mask)
    X, X_bin, gene_to_idx, broad_labels = prepare_expression(adata)
    # plot_neighborhood_composition(neighbor_idx, neighbor_mask, broad_labels)

    # plot_neighborhood_example(adata, neighbor_idx, samples_arr)

    # ── Count LR pairs ─────────────────────────────────────────────────────────
    lr_counts = count_lr_pairs(lr_pairs, neighbor_idx, neighbor_mask, X_bin, gene_to_idx, broad_labels)
    # plot_celltype_pair_distribution(neighbor_idx, neighbor_mask, broad_labels)

    # ── Archetype scoring ──────────────────────────────────────────────────────
    adata = score_visium_heart(adata, lr_counts, archetype_genes)
    # plot_ligand_receptor_archetype_scatter(adata, lr_counts)
    #plot_arch_vs_gene_logratio(adata, lr_counts)
    # plot_archetype_score_distributions(adata)
    lr_counts = compute_pearson_r(adata, lr_counts)
    plot_arch_vs_gene_logratio(adata, lr_counts, ligand='APP', receptor='RPSA', archetype1='AE1', archetype2='AF4')
    plot_spatial_lr_pair(adata, lr_counts,
                         pair_ligand='APP', pair_receptor='RPSA',
                         selected_arch_pairs=[('AE1', 'AF4')])
    # ── Permutation test 1: neighbor-count ────────────────────────────────────
    lr_counts, perm_counts, observed = run_permutation_test(
        lr_counts, X_bin, gene_to_idx, broad_labels, neighbor_idx, neighbor_mask, samples_arr, n_perms=10000)
    lr_significant = lr_counts[lr_counts['p_adj'] < 0.05]
    print(f"Significant pairs (p_adj < 0.05): {len(lr_significant)}")
    plot_obs_vs_exp(lr_counts, perm_counts, observed,
                   selected_pairs=[('APP', 'RPSA', 'AE1', 'AF4')])
    # plot_perm_dists(lr_counts, perm_counts, observed)





    # ── Permutation test 2: distance-weighted ─────────────────────────────────
    # A_nbr = build_neighbor_adjacency(neighbor_idx, neighbor_mask)
    # score_lookup = build_score_lookup(adata)
    # lr_counts, perm_scores_w, observed_w, pair_key_to_col = run_weighted_permutation_test(
    #     lr_counts, X, gene_to_idx, score_lookup, A_nbr, samples_arr, broad_labels)
    # plot_volcano_weighted(lr_counts, perm_scores_w, observed_w)
    # plot_obs_vs_exp_weighted(lr_counts, perm_scores_w, observed_w)

    # plot_weighted_perm_dists(lr_counts, perm_scores_w, pair_key_to_col)
    # mask = (lr_counts['ligand'] == 'LGALS1') & (lr_counts['receptor'] == 'PTPRC') & \
    #        (lr_counts['archetype_ligand'] == 'AE6') & (lr_counts['archetype_receptor'] == 'AM4')

    # plot_obs_vs_exp_weighted(lr_counts, perm_scores_w, observed_w,
    #                          selected_pairs=[('LGALS1', 'PTPRC', 'AE6', 'AM4')])

    lr_counts.to_csv("../heart_visium/lr_counts.csv")
    # pairs_sig = select_doubly_significant(lr_counts)
    # ── Bubble plots ───────────────────────────────────────────────────────────
    plot_bubble_plots(lr_counts)
    # van(lr_counts)

