"""Violin plots of top gene expression across cell types."""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def plot_top_gene_violin(adata, genes, cell_type_key="cell_type", save_path=None):
    """Grid of violins: expression of ``genes`` per cell type."""
    genes = list(genes)
    X = adata[:, genes].X
    X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
    long = pd.DataFrame({
        "cell_type": np.tile(adata.obs[cell_type_key].astype(str).to_numpy(), len(genes)),
        "gene": np.repeat(genes, adata.n_obs),
        "expression": X.T.reshape(-1),
    })
    g = sns.FacetGrid(long, col="gene", sharey=False, height=3)
    g.map_dataframe(sns.violinplot, x="cell_type", y="expression", cut=0)
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=90)
    fig = g.figure
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=150)
    return fig
