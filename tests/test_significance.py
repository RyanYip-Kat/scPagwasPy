import numpy as np
import pandas as pd
import scipy.sparse as sp

from scpagwas.config import PagwasConfig
from scpagwas.preprocess import SingleCellData
from scpagwas.significance import (background_correction, empirical_background_p,
                                   merge_celltype_p, module_score)


def test_module_score_shape():
    genes = pd.Index([f"g{i}" for i in range(200)])
    values = np.random.default_rng(0).random((200, 50))
    s = module_score(values, genes, list(genes[:20]), n_bins=10, ctrl=10, seed=1)
    assert s.shape == (50,)


def test_module_score_ranks_feature_genes():
    genes = pd.Index([f"g{i}" for i in range(100)])
    values = np.random.default_rng(0).random((100, 30))
    values[:5] += 5.0  # feature genes are strongly expressed
    s = module_score(values, genes, list(genes[:5]), n_bins=10, ctrl=10, seed=1)
    assert s.mean() > 1.0


def test_empirical_background_p_columns():
    obs = np.array([1.0, 2.0, 3.0])
    null = np.random.default_rng(0).random((3, 100))
    out = empirical_background_p(obs, null, np.ones(null.shape[1]))
    expected = {"raw_score", "pooled_p", "adj_p", "nlog10_pooled_p", "pooled_z"}
    assert expected <= set(out.columns)
    assert (out["pooled_p"] > 0).all() and (out["pooled_p"] <= 1).all()
    assert (out["pooled_z"].abs() <= 10).all()


def test_merge_celltype_p():
    p = pd.Series([0.01, 0.02, 0.5, 0.6])
    ct = pd.Series(["A", "A", "B", "B"])
    out = merge_celltype_p(p, ct)
    assert set(out["celltype"]) == {"A", "B"}
    assert out.loc[out["celltype"] == "A", "pvalue"].iloc[0] < 0.05


def test_background_correction_columns():
    rng = np.random.default_rng(0)
    genes = pd.Index([f"g{i}" for i in range(40)])
    cells = pd.Index([f"c{i}" for i in range(50)])
    mat = sp.csr_matrix(rng.random((40, 50)))
    ct = pd.Series(["A"] * 25 + ["B"] * 25, index=cells)
    sc = SingleCellData(genes, cells, ct, mat, pd.DataFrame(index=genes), ["A", "B"])
    cfg = PagwasConfig(iters_singlecell=5, n_topgenes=5, seed=1)
    out = background_correction(sc, list(genes[:5]), cfg)
    assert len(out) == 50
    assert "pooled_p" in out.columns
