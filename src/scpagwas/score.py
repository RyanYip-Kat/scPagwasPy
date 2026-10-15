"""Per-cell gPAS scoring: pathway coefficients x activity x expression, summed."""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .utils import winsorize


@dataclass
class CellScores:
    gpas: pd.Series
    pathway_cell: pd.DataFrame
    rank_pvalue: pd.DataFrame


def compute_gpas(pathway_cell_lm, pca_cell_mat, data_mat, genes, pathway_list, cfg):
    """Compute the summed per-cell gPAS score.

    ``pathway_cell_lm`` is pathways x cells; ``data_mat`` is genes x cells.
    """
    names = [p for p in pathway_cell_lm.index if p in pca_cell_mat.index]
    gpos = {g: i for i, g in enumerate(genes)}
    expr = {}
    for p in names:
        gl = [g for g in pathway_list[p] if g in gpos]
        if gl:
            expr[p] = data_mat[[gpos[g] for g in gl]].mean(axis=0)
    expr = pd.DataFrame(expr, index=pca_cell_mat.columns)
    pathway_cell = (pathway_cell_lm.loc[names].T * pca_cell_mat.loc[names].T * expr).T
    gpas = pathway_cell.sum(axis=0)
    if cfg.remove_outlier:
        gpas = pd.Series(winsorize(gpas.to_numpy()), index=gpas.index)
    return CellScores(gpas, pathway_cell, pd.DataFrame(index=pathway_cell.index))
