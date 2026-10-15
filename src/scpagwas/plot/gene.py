"""Gene-level PCC scatter plot."""

import matplotlib.pyplot as plt
import numpy as np


def plot_gene_pcc(gene_pcc, topn=10, save_path=None):
    """Scatter of ranked gene PCC values, labeling the top ``topn`` genes."""
    pcc = gene_pcc.sort_values(ascending=False)
    values = pcc.to_numpy()
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(np.arange(len(values)), values, c=values, cmap="coolwarm", s=8)
    for i in range(min(topn, len(values))):
        ax.annotate(str(pcc.index[i]), (i, values[i]), fontsize=6)
    ax.axhline(0, color="#888888", linestyle="--", linewidth=0.8)
    ax.set_xlabel("gene rank")
    ax.set_ylabel("PCC")
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=150)
    return fig
