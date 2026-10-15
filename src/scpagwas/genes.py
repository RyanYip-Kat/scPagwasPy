"""Gene-level heritability correlation with the per-cell gPAS score."""

import numpy as np
import pandas as pd


def gene_pcc(gpas, data_mat, genes, cells):
    """Pearson correlation of each gene across cells with the gPAS score."""
    g = np.asarray(gpas.reindex(cells).to_numpy(), dtype=float)
    g = (g - g.mean()) / (g.std() + 1e-12)
    X = np.asarray(data_mat, dtype=float)
    Xc = X - X.mean(axis=1, keepdims=True)
    Xs = Xc / (Xc.std(axis=1, keepdims=True) + 1e-12)
    return pd.Series(Xs @ g / len(g), index=genes)
