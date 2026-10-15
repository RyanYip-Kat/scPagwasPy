"""Regressions linking GWAS effect size to pathway activity x expression.

The model mirrors scPagwas: ``beta^2 ~ offset(se^2) + design`` where the
design columns are pathway-activity-weighted, LD-weighted gene expression.
"""

import numpy as np
import pandas as pd
from scipy.stats import norm


def singlecell_regression(link, pca_cell_mat, data_mat, genes, pathway_genes, cells):
    """Per-cell univariate coefficients for one pathway (length = n_cells)."""
    gl = [g for g in pathway_genes if g in set(genes)]
    gene_pos = {g: i for i, g in enumerate(genes)}
    local = {g: i for i, g in enumerate(gl)}
    sub = np.asarray(data_mat[[gene_pos[g] for g in gl], :], dtype=float)
    colsum = sub.sum(axis=0)
    frac = np.divide(sub, colsum, out=np.zeros_like(sub), where=colsum != 0)
    pca = pca_cell_mat.loc[link.pathway].to_numpy()
    design = frac * pca[None, :]

    rows = np.array([local[g] for g in link.snps["label"].to_numpy()])
    X = link.ldsq @ design[rows, :]
    y = link.y - link.snps["se"].to_numpy() ** 2
    Xc = X - X.mean(axis=0)
    denom = (Xc ** 2).sum(axis=0)
    coef = (Xc * (y - y.mean())[:, None]).sum(axis=0) / np.where(denom == 0, np.nan, denom)
    return np.nan_to_num(coef)


def build_celltype_design(link, pca_cell_df, merge_expr, pathway_genes):
    """Return ``(x, y, noise)`` for the celltype-level regression."""
    genes = [g for g in pathway_genes if g in merge_expr.index]
    sub = merge_expr.loc[genes]
    sub = sub.apply(lambda c: (c - c.min()) / (c.max() - c.min() + 1e-12))
    gene_expr = sub.reindex(link.snps["label"].to_numpy()).to_numpy()
    pca = pca_cell_df.loc[link.pathway].to_numpy()
    x = link.ldsq @ (gene_expr * pca[None, :])
    return x, link.y, link.snps["se"].to_numpy() ** 2


def _ols_offset(X, y, noise):
    A = np.column_stack([np.ones(len(y)), X])
    coef, *_ = np.linalg.lstsq(A, y - noise, rcond=None)
    return coef


def celltype_regression(designs, columns=None):
    """Global regression across pathway blocks; one coefficient per celltype."""
    X = np.vstack([d[0] for d in designs])
    y = np.concatenate([d[1] for d in designs])
    noise = np.concatenate([d[2] for d in designs])
    coef = _ols_offset(X, y, noise)
    index = ["Intercept"] + list(columns) if columns is not None else None
    return pd.Series(coef, index=index)


def bootstrap_celltype(designs, celltypes, base_params, cfg):
    """Bootstrap over pathway blocks; ``base_params`` is the full vector."""
    rng = np.random.default_rng(cfg.seed)
    n = len(designs)
    draws = []
    for _ in range(cfg.iters_celltype):
        pick = rng.choice(n, size=max(1, int(n * 0.5)), replace=False)
        X = np.vstack([designs[i][0] for i in pick])
        y = np.concatenate([designs[i][1] for i in pick])
        noise = np.concatenate([designs[i][2] for i in pick])
        draws.append(_ols_offset(X, y, noise))
    draws = np.vstack(draws)
    names = ["Intercept"] + list(celltypes)
    est = draws.mean(axis=0)
    err = draws.std(axis=0, ddof=1)
    bt = est / np.where(err == 0, np.nan, err)
    out = pd.DataFrame({
        "annotation": names,
        "bootstrap_estimate": est,
        "bootstrap_error": err,
        "bt_value": bt,
        "bp_value": norm.sf(bt),
        "bias_corrected_estimate": 2 * np.asarray(base_params) - est,
        "CI_lo": np.quantile(draws, 0.025, axis=0),
        "CI_hi": np.quantile(draws, 0.975, axis=0),
    })
    return out.set_index("annotation")
