"""Single-cell input preparation: cell filtering and pseudobulk aggregation."""

from dataclasses import dataclass

import numpy as np
import pandas as pd
import scipy.sparse as sp


@dataclass
class SingleCellData:
    """Prepared single-cell matrices.

    ``data_mat`` is genes x cells (sparse CSR); ``merge_expr`` is
    genes x celltypes pseudobulk means.
    """

    genes: pd.Index
    cells: pd.Index
    cell_type: pd.Series
    data_mat: sp.csr_matrix
    merge_expr: pd.DataFrame
    celltypes: list


def prepare_singlecell(adata, pathways, cfg, layer=None):
    """Filter small cell types and build the gene x cell matrix + pseudobulk."""
    X = adata.layers[layer] if layer else adata.X
    X = sp.csr_matrix(X).T.tocsr()
    genes = pd.Index(adata.var_names)
    if genes.has_duplicates:
        keep = ~genes.duplicated()
        X = X[keep, :]
        genes = genes[keep]
    cells = pd.Index(adata.obs_names)
    ct = pd.Series(adata.obs["cell_type"].astype(str).values, index=cells)

    counts = ct.value_counts()
    keep_types = sorted(counts.index[counts > cfg.min_clustercells])
    keep_cells = cells[ct.isin(keep_types)]
    col = [cells.get_loc(c) for c in keep_cells]
    data_mat = X[:, col].tocsr()
    cells = pd.Index(keep_cells)
    ct = ct.loc[cells]

    overlap = set(genes) & {g for gs in pathways.values() for g in gs}
    if len(overlap) < cfg.min_shared_genes:
        raise ValueError(
            f"Too few overlapping pathway genes "
            f"({len(overlap)} < {cfg.min_shared_genes})"
        )

    merge = pd.DataFrame(
        {t: np.asarray(data_mat[:, (ct == t).values].mean(axis=1)).ravel()
         for t in keep_types},
        index=genes,
    )
    return SingleCellData(genes, cells, ct, data_mat, merge, keep_types)
