"""Bar plots: bootstrap p-values and per-celltype positive-cell fractions."""

import matplotlib.pyplot as plt
import numpy as np


def _finish(fig, save_path):
    fig.tight_layout()
    if save_path is not None:
        fig.savefig(save_path, dpi=150)
    return fig


def plot_bootstrap_bar(bootstrap_df, save_path=None):
    """Bar plot of ``-log2(bp_value)`` per celltype."""
    df = bootstrap_df.drop(index="Intercept", errors="ignore")
    df = df.sort_values("bp_value", ascending=False)
    logp = -np.log2(df["bp_value"].to_numpy())
    fig, ax = plt.subplots(figsize=(6, max(3, 0.4 * len(df))))
    ax.barh(np.arange(len(df)), logp)
    ax.axvline(4.321, color="#990000", linestyle="--")
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df.index)
    ax.set_xlabel("-log2(p)")
    return _finish(fig, save_path)


def plot_positive_cells_bar(adata, cell_type_key="cell_type", p_thre=0.05,
                            save_path=None, p_key="Random_Correct_BG_adjp"):
    """Stacked bar of significant vs non-significant fraction per celltype."""
    ct = adata.obs[cell_type_key].astype(str)
    pos = np.asarray(adata.obs[p_key], dtype=float) <= p_thre
    frame = (
        adata.obs.assign(_ct=ct, _pos=pos)
        .groupby("_ct")["_pos"]
        .mean()
        .sort_values(ascending=False)
    )
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(frame.index.astype(str), frame.to_numpy())
    ax.tick_params(axis="x", rotation=90)
    ax.set_ylabel(f"fraction adj p <= {p_thre}")
    return _finish(fig, save_path)
