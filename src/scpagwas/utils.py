"""Numerical helpers used across scpagwas."""

import numpy as np


def winsorize(values, lower=0.01, upper=0.99):
    """Clip values to the [lower, upper] quantile range; NaN/inf become finite."""
    v = np.asarray(values, dtype=float).copy()
    v[~np.isfinite(v)] = 0.0
    lo, hi = np.quantile(v, lower), np.quantile(v, upper)
    return np.clip(v, lo, hi)


def storey_qvalue(pvals, lambdas=None):
    """Storey q-values.

    Returns ``(qvalues, pi0)`` with NaN preserved for non-finite inputs.
    """
    p = np.asarray(pvals, dtype=float)
    finite = np.isfinite(p)
    out = np.full(p.shape, np.nan)
    pf = p[finite]
    if pf.size == 0:
        return out, 1.0
    if lambdas is None:
        lambdas = np.arange(0.05, 0.96, 0.05)
    pi0 = np.array([np.mean(pf >= lam) / (1 - lam) for lam in lambdas])
    pi0 = min(float(np.interp(lambdas[-1], lambdas, pi0)), 1.0)
    pi0 = max(pi0, 1e-8)
    m = pf.size
    order = np.argsort(-pf)
    ranked = pf[order]
    q = pi0 * ranked * m / np.arange(1, m + 1)
    q = np.minimum(np.minimum.accumulate(q), 1.0)
    inv = np.empty(m, dtype=int)
    inv[order] = np.arange(m)
    out[finite] = q[inv]
    return out, float(pi0)


def marchenko_pastur_edge(n_cells, n_genes, scale=1.0):
    """Upper support edge of the Marchenko-Pastur law for aspect ratio p/n."""
    return scale * (1.0 + np.sqrt(n_genes / n_cells)) ** 2
