"""Per-celltype pathway activity scoring via a single principal component.

Mirrors ``Pathway_pcascore_run`` / ``PathwayPCAtest``: for each pathway and
cell type, take PC1 of the centered expression submatrix, standardize the
observed variance by the Marchenko-Pastur upper edge, and keep the PC1 score
as the per-cell pathway activity.
"""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .utils import marchenko_pastur_edge


@dataclass
class PathwayScore:
    pca_cell_mat: pd.DataFrame  # pathways x cells
    pca_cell_df: pd.DataFrame   # pathways x celltypes
    pathways: dict


def pca1_scores(matrix_cells_by_genes):
    """Return (PC1 scores per row, loadings) with the sign made interpretable."""
    m = np.asarray(matrix_cells_by_genes, dtype=float)
    mc = m - m.mean(axis=0)
    u, s, vt = np.linalg.svd(mc, full_matrices=False)
    loadings = vt[0]
    scores = u[:, 0] * s[0]
    ref = (np.abs(loadings) * m.mean(axis=0)).sum()
    if np.sign(np.corrcoef(scores, mc @ np.abs(loadings))[0, 1]) != np.sign(ref):
        loadings, scores = -loadings, -scores
    return scores, loadings


def _pathway_val(mat_cells_by_genes):
    # Degenerate pathway (no variation after centering): drop it, as R does.
    if np.ptp(mat_cells_by_genes, axis=0).sum() == 0:
        return np.nan, None
    try:
        scores, _ = pca1_scores(mat_cells_by_genes)
    except Exception:
        return np.nan, None
    n_cells, n_genes = mat_cells_by_genes.shape
    return float(np.var(scores)) / marchenko_pastur_edge(n_cells, n_genes), scores


def score_pathways(sc, pathways, cfg):
    """Score every pathway in every cell type and per cell."""
    genes = set(sc.genes)
    pw = {k: list(dict.fromkeys(g for g in v if g in genes))
          for k, v in pathways.items()}
    pw = {k: v for k, v in pw.items()
          if cfg.min_pathway_size <= len(v) <= cfg.max_pathway_size}
    if not pw:
        raise ValueError("No pathways passed the size filter")

    dense = np.asarray(sc.data_mat.todense())  # genes x cells
    cell_scores, type_scores = {}, {}
    for i, (name, gl) in enumerate(pw.items()):
        idx = [sc.genes.get_loc(g) for g in gl]
        rows = dense[idx, :]
        row = np.full(len(sc.cells), np.nan)
        vals = []
        for t in sc.celltypes:
            mask = (sc.cell_type == t).values
            sub = rows[:, mask].T
            if sub.shape[0] < 3 or sub.shape[1] < 2:
                vals.append(np.nan)
                continue
            val, cell_vals = _pathway_val(sub)
            vals.append(val)
            if cell_vals is not None:
                row[mask] = cell_vals
        cell_scores[name] = row
        type_scores[name] = vals
        if i % 50 == 0:
            print(f"  pathway {i + 1}/{len(pw)}", flush=True)

    pca_cell_mat = pd.DataFrame(cell_scores, index=sc.cells).T
    pca_cell_df = pd.DataFrame(type_scores, index=sc.celltypes).T
    pca_cell_df.columns = sc.celltypes
    keep = pca_cell_df.notna().all(axis=1)
    pca_cell_df = pca_cell_df.loc[keep]
    pca_cell_mat = pca_cell_mat.loc[keep].fillna(0.0)
    pw = {k: pw[k] for k in pca_cell_mat.index}
    return PathwayScore(pca_cell_mat, pca_cell_df, pw)
