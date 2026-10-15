"""Embedding plots of scpagwas cell scores and significant cells."""

import matplotlib.pyplot as plt
import numpy as np


def _finish(fig, save_path):
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=150)
    return fig


def plot_umap_scores(adata, key="scPagwas.gPAS.score", obsm="X_umap", save_path=None):
    """Scatter cells in an embedding colored by a per-cell score."""
    emb = np.asarray(adata.obsm[obsm])
    values = np.asarray(adata.obs[key], dtype=float)
    fig, ax = plt.subplots(figsize=(5, 5))
    sc = ax.scatter(emb[:, 0], emb[:, 1], c=values, s=6, cmap="viridis")
    fig.colorbar(sc, ax=ax, label=key)
    ax.set_title(key)
    ax.set_xticks([])
    ax.set_yticks([])
    return _finish(fig, save_path)


def plot_significant_cells(adata, p_thre=0.05, obsm="X_umap", save_path=None,
                           key="Random_Correct_BG_adjp"):
    """Highlight cells whose adjusted p-value is below ``p_thre``."""
    emb = np.asarray(adata.obsm[obsm])
    sig = np.asarray(adata.obs[key], dtype=float) <= p_thre
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.scatter(emb[~sig, 0], emb[~sig, 1], c="#d9d9d9", s=6, label="not significant")
    ax.scatter(emb[sig, 0], emb[sig, 1], c="#ea5455", s=6,
               label=f"adj p <= {p_thre}")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.legend(loc="best", fontsize=8)
    ax.set_title("significant cells")
    return _finish(fig, save_path)
