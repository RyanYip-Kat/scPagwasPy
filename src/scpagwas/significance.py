"""Cell-level significance: expression-matched module scores and an empirical
background correction, plus merging cell p-values into celltype p-values.

This mirrors scPagwas ``AddModuleScore`` + ``Get_CorrectBg_p`` +
``Merge_celltype_p``.
"""

import numpy as np
import pandas as pd
from scipy.stats import norm


def _bins(mean_expr, n_bins):
    order = np.argsort(mean_expr)
    ranks = np.empty_like(order)
    ranks[order] = np.arange(len(mean_expr))
    return (ranks * n_bins // max(1, len(mean_expr))).clip(0, n_bins - 1)


def module_score(values, genes, feature_genes, n_bins=24, ctrl=100, seed=0):
    """Mean feature-gene expression minus expression-matched control genes.

    ``values`` is genes x cells.
    """
    values = np.asarray(values, dtype=float)
    bins = _bins(values.mean(axis=1), n_bins)
    rng = np.random.default_rng(seed)
    idx = {g: i for i, g in enumerate(genes)}
    feat = [idx[g] for g in feature_genes if g in idx]
    ctrl_idx = []
    for fi in feat:
        pool = np.where(bins == bins[fi])[0]
        pool = pool[pool != fi]
        if pool.size:
            ctrl_idx.append(rng.choice(pool, size=min(ctrl, pool.size), replace=False))
    ctrl_idx = np.concatenate(ctrl_idx) if ctrl_idx else np.array([], dtype=int)
    base = values[ctrl_idx].mean(axis=0) if ctrl_idx.size else 0.0
    return values[feat].mean(axis=0) - base


def _empirical_p(obs, null):
    null = np.sort(null)
    pos = np.searchsorted(null, obs, side="right")
    return (null.size - pos + 1) / (null.size + 1)


def empirical_background_p(scores, null_matrix, v_ratio):
    """One-sided empirical p-values against a variance-matched null."""
    obs = np.asarray(scores, dtype=float)
    obs = obs - obs.mean()
    null = np.asarray(null_matrix, dtype=float) / np.sqrt(np.asarray(v_ratio))[None, :]
    null = null - null.mean(axis=0, keepdims=True)
    pooled = _empirical_p(obs, null.ravel())
    order = np.argsort(-pooled)
    q = np.minimum.accumulate(pooled[order] * len(pooled) / np.arange(1, len(pooled) + 1))
    adj = np.empty_like(q)
    adj[order] = np.minimum(q, 1.0)
    return pd.DataFrame({
        "raw_score": obs,
        "pooled_p": pooled,
        "adj_p": adj,
        "nlog10_pooled_p": -np.log10(pooled),
        "pooled_z": np.clip(-norm.ppf(pooled), -10, 10),
    })


def background_correction(sc, feature_genes, cfg):
    """Module score for ``feature_genes`` corrected by random gene-set nulls."""
    dense = np.asarray(sc.data_mat.todense())
    rng = np.random.default_rng(cfg.seed)
    n_top = len(feature_genes)
    null = np.zeros((len(sc.cells), cfg.iters_singlecell))
    ratios = np.ones(cfg.iters_singlecell)
    var_all = dense.var(axis=1)
    tvar = var_all[[sc.genes.get_loc(g) for g in feature_genes if g in sc.genes]].sum()
    for i in range(cfg.iters_singlecell):
        ctrl = rng.choice(len(sc.genes), size=n_top, replace=False)
        null[:, i] = module_score(dense, sc.genes, list(sc.genes[ctrl]),
                                  seed=int(rng.integers(1e6)))
        ratios[i] = var_all[ctrl].sum() / (tvar + 1e-12)
    out = empirical_background_p(module_score(dense, sc.genes, feature_genes), null, ratios)
    out.index = sc.cells
    return out


def merge_celltype_p(pvalues, cell_type):
    """Stouffer-Z merge of per-cell p-values within each cell type."""
    p = pd.Series(np.asarray(pvalues, dtype=float))
    ct = pd.Series(np.asarray(cell_type))
    z = -np.sqrt(2) * norm.ppf(p / 2)
    agg = pd.DataFrame({"celltype": ct, "z": z}).groupby("celltype")["z"].mean()
    return pd.DataFrame({
        "celltype": agg.index,
        "pvalue": norm.cdf(-agg.abs().to_numpy()),
    }).reset_index(drop=True)
