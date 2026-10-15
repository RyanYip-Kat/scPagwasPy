import numpy as np
import pandas as pd
import scipy.sparse as sp

from scpagwas.config import PagwasConfig
from scpagwas.pathway_score import pca1_scores, score_pathways
from scpagwas.preprocess import SingleCellData


def _sc():
    rng = np.random.default_rng(1)
    genes = pd.Index([f"G{i}" for i in range(40)])
    cells = pd.Index([f"c{i}" for i in range(40)])
    mat = sp.csr_matrix(rng.random((40, 40)))
    ct = pd.Series(["A"] * 20 + ["B"] * 20, index=cells)
    merge = pd.DataFrame(rng.random((40, 2)), index=genes, columns=["A", "B"])
    return SingleCellData(genes, cells, ct, mat, merge, ["A", "B"])


def test_pca1_shapes():
    scores, loadings = pca1_scores(np.random.default_rng(2).random((30, 8)))
    assert scores.shape == (30,)
    assert loadings.shape == (8,)


def test_score_pathways_shapes():
    pw = {"p1": [f"G{i}" for i in range(20)], "p2": [f"G{i}" for i in range(20, 40)]}
    cfg = PagwasConfig(min_pathway_size=2, max_pathway_size=100)
    ps = score_pathways(_sc(), pw, cfg)
    assert set(ps.pca_cell_mat.index) == {"p1", "p2"}
    assert ps.pca_cell_mat.shape[1] == 40
    assert ps.pca_cell_df.shape == (2, 2)
    assert (ps.pca_cell_df.values >= 0).all()


def test_score_pathways_skips_degenerate_pathway():
    rng = np.random.default_rng(3)
    genes = pd.Index([f"G{i}" for i in range(40)])
    cells = pd.Index([f"c{i}" for i in range(40)])
    mat = np.zeros((40, 40))
    mat[20:] = rng.random((20, 40))  # only the second pathway's genes vary
    ct = pd.Series(["A"] * 20 + ["B"] * 20, index=cells)
    merge = pd.DataFrame(rng.random((40, 2)), index=genes, columns=["A", "B"])
    sc = SingleCellData(genes, cells, ct, sp.csr_matrix(mat), merge, ["A", "B"])
    pw = {"p1": list(genes[:20]), "p2": list(genes[20:])}
    cfg = PagwasConfig(min_pathway_size=2, max_pathway_size=100)
    ps = score_pathways(sc, pw, cfg)
    assert "p1" not in ps.pca_cell_mat.index
    assert "p2" in ps.pca_cell_mat.index


def test_score_pathways_dedupes_pathway_genes():
    pw = {"p1": ["G0", "G0", "G1", "G2", "G3"],
          "p2": [f"G{i}" for i in range(4, 24)]}
    cfg = PagwasConfig(min_pathway_size=2, max_pathway_size=100)
    ps = score_pathways(_sc(), pw, cfg)
    assert len(ps.pathways["p1"]) == len(set(ps.pathways["p1"]))
